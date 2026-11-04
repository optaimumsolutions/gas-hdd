#!/usr/bin/env python3
"""
hdd.py - US gas-weighted heating degree days (HDD) from NOAA's GEFS ensemble.

What this script does, in order:
  1. Pick the newest GEFS run that is COMPLETE on NOAA's public S3 bucket
     (all 31 members x all 65 lead times), or the run given with --run.
  2. Download 2-metre temperature (6-hourly value, 6-h max and 6-h min) for
     every member and lead time. herbie pulls just those three fields out of
     each GRIB file (byte-range requests), not the whole file.
  3. Cut each field down to the Lower 48 and stack everything into arrays
     T, TX, TN [member, lead_time, lat, lon] in Kelvin.
  4. Turn the daily (Tmax + Tmin) / 2 at 16 big cities into HDDs, weight them by
     population, and summarise across the 31 members (median, p10, p90).
     Then, for each region in regions.py, weight that region's cities by gas
     consumption and convert HDD to Bcf/d of res/comm demand, member by member.
  5. Compare the 15-day totals with the previous run and the run 24 h earlier.
  6. Archive each airport's daily Tmax/Tmin for all 31 members in Neon
     (gefs.forecast_station; see PRD-calibration.md). With --no-docs, stop here.
  7. Draw the charts and maps (national, then per region).
  8. Write docs/latest.json, the PNGs, and the per-member CSVs in history/.

Nothing touches disk (other than the GRIB cache) until every number and
picture has been computed, so a failure part-way through never leaves a
half-updated dashboard.

Usage:
  python hdd.py                      # newest complete run
  python hdd.py --run 2026-11-03T00  # a specific run (UTC cycle hour)
  python hdd.py --run 2025-01-15T00 --no-docs   # backfill: archive only
"""

import argparse
import io
import json
import os
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw to files; no display needed
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from herbie import Herbie

import db
from regions import REGIONS, calibrated_band, demand_bcfd, implied_withdrawal, scenario_10pct

sys.stdout.reconfigure(encoding="utf-8")  # lets the degree sign print on Windows too

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
HISTORY = ROOT / "history"
CACHE = ROOT / "data"  # downloaded GRIB pieces land here (git-ignored)

# ----------------------------------------------------------------- locations
# HDD cities: name, latitude, longitude, rough population weight (millions),
# and the airport ASOS station (ICAO) that observe.py scores the forecast
# against. The weights are normalised to sum to 1 in weighted_hdd().
CITIES = [
    ("New York",      40.71,  -74.01, 8.3, "LGA"),
    ("Los Angeles",   34.05, -118.24, 3.9, "LAX"),
    ("Chicago",       41.88,  -87.63, 2.7, "ORD"),
    ("Houston",       29.76,  -95.37, 2.3, "IAH"),
    ("Philadelphia",  39.95,  -75.17, 1.6, "PHL"),
    ("Dallas",        32.78,  -96.80, 1.3, "DFW"),
    ("Columbus",      39.96,  -83.00, 0.9, "CMH"),
    ("Boston",        42.36,  -71.06, 0.7, "BOS"),
    ("Detroit",       42.33,  -83.05, 0.6, "DTW"),
    ("Denver",        39.74, -104.99, 0.7, "DEN"),
    ("Washington DC", 38.91,  -77.04, 0.7, "DCA"),
    ("Minneapolis",   44.98,  -93.27, 0.4, "MSP"),
    ("Pittsburgh",    40.44,  -79.99, 0.3, "PIT"),
    ("Kansas City",   39.10,  -94.58, 0.5, "MCI"),
    ("St. Louis",     38.63,  -90.20, 0.3, "STL"),
    ("Atlanta",       33.75,  -84.39, 0.5, "ATL"),
]

# Every airport the forecast is archived at: ICAO -> (lat, lon). The regional
# tables carry the airport's own coordinates, the national table the city
# centre; the airport wins because that is where the thermometer is.
STATIONS = {c[5]: (c[2], c[3]) for r in REGIONS for c in r["cities"]}
for _c in CITIES:
    STATIONS.setdefault(_c[4], (_c[1], _c[2]))
STATIONS = dict(sorted(STATIONS.items()))

# Gas-producing basins: name, map label, latitude, longitude, freeze-off
# threshold in F. A basin is flagged when >= 30% of members take the period's
# minimum temperature below its threshold.
BASINS = [
    ("Permian / Midland TX",        "PERM", 32.00, -102.08, 20),
    ("Haynesville / Shreveport LA", "HAYN", 32.53,  -93.75, 20),
    ("Anadarko / Oklahoma City",    "ANAD", 35.47,  -97.52, 20),
    ("Eagle Ford / San Antonio",    "EFRD", 29.42,  -98.49, 20),
    ("Appalachia / Pittsburgh",     "APPA", 40.44,  -79.99,  0),
    ("DJ / Greeley CO",             "DJ",   40.42, -104.71,  0),
    ("Bakken / Williston ND",       "BAKK", 48.15, -103.62,  0),
]

# ---------------------------------------------------------------- GEFS setup
MEMBERS = list(range(31))          # 0 = control run, 1..30 = perturbed members
FXX = list(range(0, 385, 6))       # lead times in hours: 0, 6, ..., 384 (65 of them)
# The GRIB fields we want: 2 m temperature at the valid time, and its max and
# min over the 6 hours ending at the valid time (absent at f000). One regex, so
# each file is still one byte-range download.
FIELDS = [":TMP:2 m above ground:", ":TMAX:2 m above ground:", ":TMIN:2 m above ground:"]
SEARCH = ":(?:TMP|TMAX|TMIN):2 m above ground:"
LAT_MIN, LAT_MAX = 24.0, 50.0      # Lower-48 bounding box
LON_MIN, LON_MAX = -125.0, -66.0
THREADS = int(os.environ.get("THREADS", 8))   # parallel S3 requests (more = more "Slow Down"); backfill sets 4
FREEZE_F = 32.0                    # the map's freeze line
HDD_BASE_F = 65.0                  # HDD = max(0, 65 - daily mean temperature)
PERIODS = [("d1-5", 0, 5), ("d6-10", 5, 10), ("d11-15", 10, 15)]  # map panels


class step:
    """Context manager that prints a step's name and timing, and re-raises any
    error with the step name in front so the Actions log says what broke."""

    def __init__(self, name):
        self.name = name

    def __enter__(self):
        print(f"\n[{self.name}]", flush=True)
        self.t0 = time.time()

    def __exit__(self, exc_type, exc, tb):
        if exc is not None:
            raise RuntimeError(f"step '{self.name}' failed: {exc}") from exc
        print(f"[{self.name}] done in {time.time() - self.t0:.0f}s", flush=True)


def retry(call):
    """Run call() and, if it raises, wait and try again (up to 5 attempts).
    NOAA's S3 bucket sometimes answers "503 Slow Down" or drops a connection;
    the second try almost always works."""
    for attempt in range(5):
        try:
            return call()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt + random.random())


# ------------------------------------------------------------ 1. pick a run
def find_files(run, pairs):
    """Make one Herbie object per (member, lead time) pair. Creating it checks
    that the GRIB file and its .idx index exist on S3. Returns the dict of
    objects and the list of pairs that are missing."""

    def make(mf):
        m, f = mf
        H = retry(lambda: Herbie(run, model="gefs", product="atmos.5", member=m, fxx=f,
                                 priority="aws", save_dir=CACHE, verbose=False))
        return mf, H

    with ThreadPoolExecutor(THREADS) as pool:
        files = dict(pool.map(make, pairs))
    missing = [mf for mf, H in files.items() if H.grib is None or H.idx is None]
    if missing:
        # S3 sometimes answers "503 Slow Down", and herbie then reports the file
        # as absent. Wait a moment and look again before believing it is missing.
        time.sleep(10)
        with ThreadPoolExecutor(THREADS) as pool:
            files.update(dict(pool.map(make, missing)))
        missing = [mf for mf in missing if files[mf].grib is None or files[mf].idx is None]
    return files, missing


def probe(run):
    """Check whether a run is complete. GEFS writes lead times in order, so the
    last lead time is checked first for every member (31 files, seconds); only
    if all of those exist are all 2015 files checked (minutes)."""
    files, missing = find_files(run, [(m, FXX[-1]) for m in MEMBERS])
    if missing:
        return files, missing
    return find_files(run, [(m, f) for m in MEMBERS for f in FXX])


def pick_run(explicit):
    """Return (run datetime, Herbie objects, note). With --run, that run must be
    complete. Otherwise try the newest cycle and step back 6 h at a time."""
    if explicit is not None:
        files, missing = probe(explicit)
        if missing:
            m, f = missing[0]
            raise RuntimeError(f"run {explicit:%Y-%m-%d %H}Z is missing {len(missing)} "
                               f"of {len(files)} files checked on S3 (first: member {m}, f{f:03d})")
        return explicit, files, ""

    now = datetime.now(timezone.utc).replace(tzinfo=None)
    latest = now.replace(minute=0, second=0, microsecond=0)
    latest = latest.replace(hour=latest.hour - latest.hour % 6)  # 00/06/12/18Z
    note = ""
    for k in range(4):  # the newest cycle and the three before it
        run = latest - timedelta(hours=6 * k)
        files, missing = probe(run)
        if not missing:
            print(f"  using run {run:%Y-%m-%d %H}Z (all {len(files)} files present)")
            return run, files, note
        print(f"  run {run:%Y-%m-%d %H}Z: {len(missing)} of {len(files)} files checked not on S3 yet")
        note += f"{run:%Y-%m-%d %H}Z was incomplete at {now:%H:%M} UTC. "
    raise RuntimeError("no complete GEFS run in the last 24 hours")


# ------------------------------------------------------------- 2. download
def download_all(files):
    """Fetch the three 2 m temperature fields for every file, THREADS at a time."""
    t0 = time.time()
    with ThreadPoolExecutor(THREADS) as pool:
        list(pool.map(lambda H: retry(lambda: H.download(SEARCH)), files.values()))
    print(f"  downloaded {len(files)} files x {len(FIELDS)} fields in {time.time() - t0:.0f}s")


# ------------------------------------------------------ 3. build the array
def load_cube(files):
    """Open every downloaded file and stack the fields into T, TX, TN
    [member, lead, lat, lon] (Kelvin), keeping only grid cells inside the
    Lower-48 box. TX/TN are the 6-h max/min ending at each lead; where a file
    has no TMAX/TMIN (always at f000) the 6-hourly TMP stands in for both."""
    T = None
    for mi, m in enumerate(MEMBERS):
        for fi, f in enumerate(FXX):
            ds = files[(m, f)].xarray(SEARCH, remove_grib=False)
            if T is None:  # first file: work out which rows/columns to keep
                lat_all = ds.latitude.values                        # 90 .. -90
                lon_all = (ds.longitude.values + 180) % 360 - 180   # 0..360 -> -180..180
                ilat = np.where((lat_all >= LAT_MIN) & (lat_all <= LAT_MAX))[0]
                ilon = np.where((lon_all >= LON_MIN) & (lon_all <= LON_MAX))[0]
                lat, lon = lat_all[ilat], lon_all[ilon]
                T = np.full((len(MEMBERS), len(FXX), len(ilat), len(ilon)), np.nan, np.float32)
                TX, TN = T.copy(), T.copy()
            # Guard against a mixed-up file: the GRIB header must say this member and lead.
            if int(ds.number) != m or int(ds.step.values / np.timedelta64(1, "h")) != f:
                raise RuntimeError(f"file for member {m} f{f:03d} holds member "
                                   f"{int(ds.number)} step {ds.step.values}")
            T[mi, fi] = ds.t2m.values[np.ix_(ilat, ilon)]
            if "tmax" in ds and "tmin" in ds:
                TX[mi, fi] = ds.tmax.values[np.ix_(ilat, ilon)]
                TN[mi, fi] = ds.tmin.values[np.ix_(ilat, ilon)]
            ds.close()
        print(f"  member {m:2d} loaded", flush=True)
    if np.isnan(T).any():
        raise RuntimeError("temperature array has gaps")
    for A in (T, TX, TN):
        if np.nanmin(A) < 180 or np.nanmax(A) > 340:  # Kelvin sanity check
            raise RuntimeError(f"temperatures out of range: {np.nanmin(A):.0f}..{np.nanmax(A):.0f} K")
    missing = np.isnan(TX[:, :, 0, 0]) | np.isnan(TN[:, :, 0, 0])   # [member, lead]
    for fi in np.where(missing.any(axis=0))[0]:
        if FXX[fi] == 0:
            continue  # f000 never has TMAX/TMIN, and no HDD day uses it
        members = [MEMBERS[mi] for mi in np.where(missing[:, fi])[0]]
        print(f"  WARNING: no TMAX/TMIN at f{FXX[fi]:03d} for members {members}; using 6-h TMP")
    TX = np.where(np.isnan(TX), T, TX)
    TN = np.where(np.isnan(TN), T, TN)
    return T, TX, TN, lat, lon


# ------------------------------------------------------------ 4. the maths
def k_to_f(kelvin):
    return (kelvin - 273.15) * 9 / 5 + 32


def nearest_index(lat, lon, point_lat, point_lon):
    """Row and column of the 0.5-degree grid cell closest to a point."""
    return int(np.abs(lat - point_lat).argmin()), int(np.abs(lon - point_lon).argmin())


def full_days(run):
    """Group the 65 valid times by UTC calendar day and keep the first 15 days
    that have all four 6-hourly values (00, 06, 12, 18 UTC). Returns the list
    of dates and, for each date, the lead-time indices that fall in it."""
    groups = {}
    for fi, f in enumerate(FXX):
        day = (run + timedelta(hours=f)).date()
        groups.setdefault(day, []).append(fi)
    days = [d for d in sorted(groups) if len(groups[d]) == 4][:15]
    if len(days) != 15:
        raise RuntimeError(f"expected 15 full days, found {len(days)}")
    if groups[days[-1]][-1] + 1 >= len(FXX):   # city_hdd needs the next day's 00Z TMAX/TMIN
        raise RuntimeError("day 15 has no 24Z lead for TMAX/TMIN")
    return days, [groups[d] for d in days]


def hdd_from_mean(daily_mean_f):
    """HDD = max(0, 65 - daily mean F). Shared by the forecast and observe.py,
    so forecast and observed HDDs can never drift apart in definition."""
    return np.maximum(0.0, HDD_BASE_F - daily_mean_f)


def weighted_hdd(daily_mean_f, weights):
    """City HDDs combined into one index: sum over cities of weight x HDD, with
    the weights normalised to sum to 1. daily_mean_f[..., city] in F.
    Shared by the forecast (city_hdd) and observe.py."""
    weights = np.asarray(weights, dtype=float)
    return hdd_from_mean(daily_mean_f) @ (weights / weights.sum())


def daily_extremes(TX, TN, lat, lon, day_steps, points):
    """Daily Tmax and Tmin in F at each point: arrays [member, day, point].
    points = list of (latitude, longitude). TMAX/TMIN cover the 6 hours ENDING
    at their lead, so a UTC day 00-24Z is the four leads valid 06, 12, 18 and
    24Z: the day's own leads shifted by one."""
    rows, cols = zip(*[nearest_index(lat, lon, p[0], p[1]) for p in points])
    tmax_f = k_to_f(TX[:, :, list(rows), list(cols)])       # [member, lead, point]
    tmin_f = k_to_f(TN[:, :, list(rows), list(cols)])
    tmax = np.zeros((len(MEMBERS), len(day_steps), len(points)))
    tmin = tmax.copy()
    for di, steps in enumerate(day_steps):
        window = [s + 1 for s in steps]                     # leads valid 06, 12, 18, 24Z
        tmax[:, di] = tmax_f[:, window].max(axis=1)
        tmin[:, di] = tmin_f[:, window].min(axis=1)
    return tmax, tmin


def city_hdd(TX, TN, lat, lon, day_steps, cities):
    """hdd[member, day] = sum over cities of weight x max(0, 65 - daily mean F).
    cities = list of (latitude, longitude, weight).
    Daily mean = (Tmax + Tmin) / 2, the EIA/NOAA definition."""
    tmax, tmin = daily_extremes(TX, TN, lat, lon, day_steps, [(c[0], c[1]) for c in cities])
    hdd = np.zeros((len(MEMBERS), len(day_steps)))
    for di in range(len(day_steps)):
        hdd[:, di] = weighted_hdd((tmax[:, di] + tmin[:, di]) / 2, [c[2] for c in cities])
    return hdd


def period_fields(T, lat, lon, day_steps, first, last):
    """For one map panel (days first+1..last):
    - ensemble-mean period-average temperature (F) at every grid cell,
    - 'uncertain' mask: 30-70% of members put the period average below 32 F,
    - per basin, the fraction of members whose period minimum is below its threshold."""
    steps = [s for d in day_steps[first:last] for s in d]
    temp_f = k_to_f(T[:, steps])                 # [member, step, lat, lon]
    period_avg = temp_f.mean(axis=1)             # [member, lat, lon]
    mean_f = period_avg.mean(axis=0)             # ensemble mean
    frac_cold = (period_avg < FREEZE_F).mean(axis=0)
    uncertain = (frac_cold >= 0.3) & (frac_cold <= 0.7)
    risk = {}
    for name, _, blat, blon, threshold in BASINS:
        r, c = nearest_index(lat, lon, blat, blon)
        period_min = temp_f[:, :, r, c].min(axis=1)   # [member]
        risk[name] = float((period_min < threshold).mean())
    return mean_f, uncertain, risk


def region_series(region, TX, TN, lat, lon, day_steps):
    """One region's per-member numbers. Everything is converted member by member
    and only then summarised, never on the percentiles:
      hdd[m, d]    consumption-weighted HDD
      bcfd[m, d]   res/comm demand, Bcf/d
      total_bcf[m] 15-day sum of bcfd, Bcf
      cal          calibrated band (see regions.calibrated_band)"""
    hdd = city_hdd(TX, TN, lat, lon, day_steps, [(c[2], c[3], c[4]) for c in region["cities"]])
    bcfd = demand_bcfd(region, hdd)
    return {"hdd": hdd, "bcfd": bcfd, "total_bcf": bcfd.sum(axis=1), "cal": calibrated_band(region, bcfd)}


def region_extent(region):
    """Map box [lon_min, lon_max, lat_min, lat_max]: the region's cities and
    basins plus 3 degrees on every side."""
    pts = [(c[2], c[3]) for c in region["cities"]] + [(b[2], b[3]) for b in BASINS if b[0] in region["basins"]]
    lats, lons = zip(*pts)
    return [min(lons) - 3, max(lons) + 3, min(lats) - 3, max(lats) + 3]


# ------------------------------------------------- 5. compare with history
def history_runs(prefix="hdd"):
    """Map run datetime -> CSV path for every history/<prefix>_YYYYMMDD_HHz.csv."""
    return {datetime.strptime(p.stem, f"{prefix}_%Y%m%d_%Hz"): p for p in HISTORY.glob(f"{prefix}_*.csv")}


def compare(values, days, other_run, csv_path):
    """Change in 15-day total versus another run, measured over the days both
    runs forecast (so a day dropping off the front of the window does not
    count as a 'change'). Totals are the median of each member's sum.
    values[member, day] is HDD for the national index, Bcf/d for a region."""
    other = pd.read_csv(csv_path, index_col="date")
    shared = [d for d in days if d.isoformat() in other.index]
    if not shared:
        return None
    now_total = np.median(values[:, [days.index(d) for d in shared]].sum(axis=1))
    other_total = np.median(other.loc[[d.isoformat() for d in shared]].to_numpy().sum(axis=0))
    return {"run": f"{other_run:%Y-%m-%d %H}Z", "delta": round(float(now_total - other_total), 1),
            "shared_days": len(shared), "other_total": round(float(other_total), 1)}


# ------------------------------------------------------------- 6. pictures
def png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def make_chart(run, days, hdd, p10, median, p90, total, title="US gas-weighted HDD"):
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    ax.fill_between(days, p10, p90, color="tab:blue", alpha=0.2, label="10th-90th percentile")
    ax.plot(days, median, color="tab:blue", lw=2.2, marker="o", ms=3.5, label="Ensemble median")
    ax.set_ylabel("Gas-weighted HDD (°F-days)")
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left")
    ax.set_title(f"{title} - GEFS {run:%Y-%m-%d %H}Z, 31 members\n"
                 f"15-day total: median {total['median']:.0f}  "
                 f"(p10 {total['p10']:.0f} / p90 {total['p90']:.0f})")
    return png_bytes(fig)


def make_bcf_chart(region, run, days, bcfd, cal, total):
    """Daily res/comm demand: median line, raw p10-p90 shaded, calibrated
    p10-p90 as a dashed outline."""
    p10, median, p90 = np.percentile(bcfd, [10, 50, 90], axis=0)
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    ax.fill_between(days, p10, p90, color="tab:orange", alpha=0.25, label="Raw p10-p90 (GEFS spread)")
    ax.plot(days, cal["daily_lo"], color="tab:red", lw=1.3, ls="--", label="Calibrated p10-p90")
    ax.plot(days, cal["daily_hi"], color="tab:red", lw=1.3, ls="--")
    ax.plot(days, median, color="tab:orange", lw=2.2, marker="o", ms=3.5, label="Ensemble median")
    ax.set_ylabel("Res/comm gas demand (Bcf/d)")
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left")
    ax.set_title(f"{region['name']} res/comm gas demand - GEFS {run:%Y-%m-%d %H}Z, 31 members\n"
                 f"15-day total: median {total['median']:.0f} Bcf  "
                 f"(calibrated {total['cal_p10']:.0f}-{total['cal_p90']:.0f}; "
                 f"raw {total['raw_p10']:.0f}-{total['raw_p90']:.0f})")
    return png_bytes(fig)


def make_map(run, first, last, days, mean_f, uncertain, risk, lat, lon,
             extent=(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX), points=None, basins=BASINS, name=""):
    """points = city (lat, lon) dots; defaults to the national CITIES."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    plt.rcParams["hatch.color"] = "#333333"
    plt.rcParams["hatch.linewidth"] = 0.6
    pc = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9, 5.8), dpi=130)
    ax = plt.axes(projection=pc)
    ax.set_extent(list(extent), crs=pc)
    levels = np.arange(-30, 101, 5)
    filled = ax.contourf(lon, lat, mean_f, levels=levels, cmap="RdYlBu_r", extend="both", transform=pc)
    ax.contour(lon, lat, mean_f, levels=[FREEZE_F], colors="black", linewidths=2, transform=pc)
    ax.contourf(lon, lat, uncertain.astype(float), levels=[0.5, 1.5], colors="none", hatches=["///"], transform=pc)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.4, edgecolor="#555555")
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.6, edgecolor="#333333")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, edgecolor="#333333")
    for clat, clon in points or [(c[1], c[2]) for c in CITIES]:
        ax.plot(clon, clat, "o", ms=4, color="black", transform=pc)
    for bname, short, blat, blon, threshold in basins:
        flagged = risk[bname] >= 0.3
        ax.plot(blon, blat, "^", ms=12, color="red" if flagged else "limegreen",
                markeredgecolor="black", transform=pc, zorder=5)
        ax.text(blon + 0.5, blat + 0.5, f"{short} {risk[bname]:.0%}", fontsize=7.5, weight="bold",
                transform=pc, zorder=6, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    fig.colorbar(filled, ax=ax, orientation="horizontal", pad=0.03, shrink=0.8,
                 label="Ensemble-mean 2 m temperature (°F)")
    ax.set_title(f"{name}GEFS {run:%Y-%m-%d %H}Z - days {first + 1}-{last} "
                 f"({days[first]:%b %d} to {days[last - 1]:%b %d})\n"
                 "black line = 32°F; hatched = <70% of members agree on below/above 32°F; "
                 "triangle red = freeze-off risk (% of members below basin threshold)", fontsize=9)
    return png_bytes(fig)


# ------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", type=lambda s: datetime.strptime(s, "%Y-%m-%dT%H"),
                        help="GEFS run to process, e.g. 2026-11-03T00 (UTC)")
    parser.add_argument("--no-docs", action="store_true",
                        help="archive the run in Neon only: no charts, no docs/, no history/ (backfill)")
    args = parser.parse_args()

    with step("1/8 find a complete GEFS run"):
        run, files, note = pick_run(args.run)

    with step("2/8 download 2 m temperature (31 members x 65 lead times)"):
        download_all(files)

    with step("3/8 read the GRIB files into one array"):
        T, TX, TN, lat, lon = load_cube(files)
        print(f"  array shape {T.shape} (member, lead, lat, lon)")

    with step("4/8 compute HDDs and map fields"):
        days, day_steps = full_days(run)
        hdd = city_hdd(TX, TN, lat, lon, day_steps, [(c[1], c[2], c[3]) for c in CITIES])  # [member, day]
        p10, median, p90 = np.percentile(hdd, [10, 50, 90], axis=0)
        member_totals = hdd.sum(axis=1)
        total = {k: round(float(v), 1) for k, v in
                 zip(["p10", "median", "p90"], np.percentile(member_totals, [10, 50, 90]))}
        st_tmax, st_tmin = daily_extremes(TX, TN, lat, lon, day_steps, list(STATIONS.values()))
        maps = {label: period_fields(T, lat, lon, day_steps, a, b) for label, a, b in PERIODS}
        regional = {r["key"]: region_series(r, TX, TN, lat, lon, day_steps) for r in REGIONS}
        for r in REGIONS:
            s = regional[r["key"]]
            raw = np.percentile(s["total_bcf"], [10, 50, 90])
            s["total"] = {"median": round(float(raw[1]), 1), "raw_p10": round(float(raw[0]), 1),
                          "raw_p90": round(float(raw[2]), 1),
                          "cal_p10": round(float(s["cal"]["total_lo"]), 1),
                          "cal_p90": round(float(s["cal"]["total_hi"]), 1)}
            h = np.percentile(s["hdd"].sum(axis=1), [10, 50, 90])
            s["total_hdd"] = {"median": round(float(h[1]), 1), "p10": round(float(h[0]), 1),
                              "p90": round(float(h[2]), 1)}
            s["hdd_median"] = np.median(s["hdd"], axis=0)
            s["scenario_10pct"] = scenario_10pct(r, s["hdd_median"])
            s["implied_withdrawal"] = implied_withdrawal(r, s["bcfd"], days)

    with step("5/8 compare with earlier runs"):
        runs = history_runs()
        earlier = [r for r in runs if r < run]
        prev_run = max(earlier) if earlier else None
        ago_run = run - timedelta(hours=24) if (run - timedelta(hours=24)) in runs else None
        vs_prev = compare(hdd, days, prev_run, runs[prev_run]) if prev_run else None
        vs_24h = compare(hdd, days, ago_run, runs[ago_run]) if ago_run else None
        for r in REGIONS:  # same comparison, in median 15-day Bcf
            s, rruns = regional[r["key"]], history_runs(f"{r['key']}_bcf")
            rprev = max([x for x in rruns if x < run], default=None)
            rago = run - timedelta(hours=24) if (run - timedelta(hours=24)) in rruns else None
            for key, other in [("vs_prev", rprev), ("vs_24h", rago)]:
                cmp = compare(s["bcfd"], days, other, rruns[other]) if other else None
                if cmp:
                    cmp["delta_bcf"] = cmp.pop("delta")
                s[key] = cmp

    with step("6/8 archive per-airport Tmax/Tmin in Neon"):
        # A failure here never blocks the published forecast (the next run
        # re-ingests this one), except in backfill mode where the archive is the
        # only output.
        try:
            with retry(db.connect) as conn:
                n = db.write_forecast(conn, run.replace(tzinfo=timezone.utc), days, STATIONS, st_tmax, st_tmin)
            print(f"  wrote {n} rows ({len(STATIONS)} airports x {len(days)} days) to gefs.forecast_station")
        except Exception as err:
            if args.no_docs:
                raise
            print(f"  WARNING: Neon archive skipped: {err}")
        if args.no_docs:
            print("  --no-docs: done")
            return

    with step("7/8 draw chart and maps"):
        pngs = {"hdd_forecast.png": make_chart(run, days, hdd, p10, median, p90, total)}
        for label, a, b in PERIODS:
            mean_f, uncertain, risk = maps[label]
            pngs[f"map_{label}.png"] = make_map(run, a, b, days, mean_f, uncertain, risk, lat, lon)
        for r in REGIONS:
            k, s = r["key"], regional[r["key"]]
            hp10, hmed, hp90 = np.percentile(s["hdd"], [10, 50, 90], axis=0)
            pngs[f"{k}_hdd.png"] = make_chart(run, days, s["hdd"], hp10, hmed, hp90, s["total_hdd"],
                                              title=f"{r['name']} gas-consumption-weighted HDD")
            pngs[f"{k}_bcf.png"] = make_bcf_chart(r, run, days, s["bcfd"], s["cal"], s["total"])
            basins = [b for b in BASINS if b[0] in r["basins"]]
            points = [(c[2], c[3]) for c in r["cities"]]
            for label, a, b in PERIODS:
                mean_f, uncertain, risk = maps[label]
                pngs[f"{k}_map_{label}.png"] = make_map(run, a, b, days, mean_f, uncertain, risk, lat, lon,
                                                        extent=region_extent(r), points=points,
                                                        basins=basins, name=f"{r['name']} - ")

    # ---- summary table to stdout (shows up in the Actions log)
    table = pd.DataFrame({"date": [d.isoformat() for d in days], "median": median, "p10": p10, "p90": p90})
    table.index = range(1, 16)
    table.index.name = "day"
    print(f"\nGEFS run {run:%Y-%m-%d %H}Z   {note}".rstrip())
    print(table.round(1).to_string())
    print(f"15-day total: median {total['median']}  (p10 {total['p10']} / p90 {total['p90']})")
    for name, cmp in [("previous run", vs_prev), ("run 24h earlier", vs_24h)]:
        if cmp is None:
            print(f"vs {name}: no prior run")
        else:
            print(f"vs {name} ({cmp['run']}): {cmp['delta']:+.1f} over {cmp['shared_days']} shared days")
    for r in REGIONS:
        s, t = regional[r["key"]], regional[r["key"]]["total"]
        print(f"\n{r['name']}: 15-day res/comm demand {t['median']} Bcf "
              f"(calibrated {t['cal_p10']}-{t['cal_p90']}, raw {t['raw_p10']}-{t['raw_p90']}); "
              f"HDD {s['total_hdd']['median']}; "
              f"sum of daily medians {np.median(s['bcfd'], axis=0).sum():.1f} Bcf")
        sc = s["scenario_10pct"]
        print(f"  10% colder: {sc['colder_bcf']:+.1f} Bcf, 10% warmer: {sc['warmer_bcf']:+.1f} Bcf")
        for w in s["implied_withdrawal"]:
            print(f"  implied withdrawal, week ending {w['week_ending']} ({w['days_in_forecast']}/7 days): "
                  f"{w['low_bcf']}-{w['high_bcf']} Bcf")
        for name, key in [("previous run", "vs_prev"), ("run 24h earlier", "vs_24h")]:
            cmp = s[key]
            print(f"  vs {name}: " + ("no prior run" if cmp is None else
                  f"{cmp['delta_bcf']:+.1f} Bcf over {cmp['shared_days']} shared days ({cmp['run']})"))

    with step("8/8 write outputs"):
        DOCS.mkdir(exist_ok=True)
        HISTORY.mkdir(exist_ok=True)
        run_id = f"{run:%Y%m%d_%H}z"
        latest = {
            "run": f"{run:%Y-%m-%d %H}Z",
            "run_id": run_id,
            "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "note": note.strip(),
            "days": [{"day": i + 1, "date": d.isoformat(), "median": round(float(median[i]), 1),
                      "p10": round(float(p10[i]), 1), "p90": round(float(p90[i]), 1)}
                     for i, d in enumerate(days)],
            "total": total,
            "vs_prev": vs_prev,
            "vs_24h": vs_24h,
            "basin_freeze_risk": {label: {k: round(v, 2) for k, v in maps[label][2].items()}
                                  for label, _, _ in PERIODS},
            "images": {"chart": "hdd_forecast.png", "maps": [f"map_{label}.png" for label, _, _ in PERIODS]},
            "regions": {},
        }
        csvs = {f"hdd_{run_id}.csv": hdd}
        for r in REGIONS:
            k, s = r["key"], regional[r["key"]]
            hp10, hmed, hp90 = np.percentile(s["hdd"], [10, 50, 90], axis=0)
            bp10, bmed, bp90 = np.percentile(s["bcfd"], [10, 50, 90], axis=0)
            latest["regions"][k] = {
                "name": r["name"], "states": r["states"], "eia_storage_region": r["eia_storage_region"],
                "days": [{"day": i + 1, "date": d.isoformat(),
                          "hdd_median": round(float(hmed[i]), 1), "hdd_p10": round(float(hp10[i]), 1),
                          "hdd_p90": round(float(hp90[i]), 1),
                          "bcfd_median": round(float(bmed[i]), 2), "bcfd_p10": round(float(bp10[i]), 2),
                          "bcfd_p90": round(float(bp90[i]), 2)}
                         for i, d in enumerate(days)],
                "total_bcf": s["total"],
                "total_hdd": s["total_hdd"],
                "vs_prev": s["vs_prev"],
                "vs_24h": s["vs_24h"],
                "scenario_10pct": s["scenario_10pct"],
                "implied_withdrawal": s["implied_withdrawal"],
                "model": {"slope": r["slope_bcfd_per_hdd"], "intercept": r["intercept_bcfd"],
                          "floor": r["floor_bcfd"], "cold_kink_hdd": r["cold_kink_hdd"],
                          "cold_kink_extra_frac": r["cold_kink_extra_frac"],
                          "storage_ratio": list(r["storage_ratio"]),
                          "slope_sigma_frac": r["slope_sigma_frac"],
                          "spread_inflation": r["spread_inflation"], "refit": r["refit"],
                          "source": "EIA STEO Supplement, Weather Sensitivity in Natural Gas Markets, Oct 2014"},
                "basin_freeze_risk": {label: {b: round(maps[label][2][b], 2) for b in r["basins"]}
                                      for label, _, _ in PERIODS},
                "images": {"hdd_chart": f"{k}_hdd.png", "bcf_chart": f"{k}_bcf.png",
                           "maps": [f"{k}_map_{label}.png" for label, _, _ in PERIODS]},
            }
            csvs[f"{k}_hdd_{run_id}.csv"] = s["hdd"]
            csvs[f"{k}_bcf_{run_id}.csv"] = s["bcfd"]
        for name, data in pngs.items():
            (DOCS / name).write_bytes(data)
        for name, values in csvs.items():
            csv = pd.DataFrame(values.T, index=[d.isoformat() for d in days],
                               columns=[f"m{m:02d}" for m in MEMBERS]).rename_axis("date")
            csv.to_csv(HISTORY / name, float_format="%.3f")
        (DOCS / "latest.json").write_text(json.dumps(latest, indent=1), encoding="utf-8")
        print(f"  wrote docs/latest.json, {len(pngs)} PNGs, {len(csvs)} CSVs in history/")


if __name__ == "__main__":
    main()
