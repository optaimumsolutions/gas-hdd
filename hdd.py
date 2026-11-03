#!/usr/bin/env python3
"""
hdd.py - US gas-weighted heating degree days (HDD) from NOAA's GEFS ensemble.

What this script does, in order:
  1. Pick the newest GEFS run that is COMPLETE on NOAA's public S3 bucket
     (all 31 members x all 65 lead times), or the run given with --run.
  2. Download 2-metre temperature for every member and lead time. herbie
     pulls just that one field out of each GRIB file (a byte-range request),
     not the whole file.
  3. Cut each field down to the Lower 48 and stack everything into one array
     T[member, lead_time, lat, lon] in Kelvin.
  4. Turn the temperatures at 16 big cities into daily HDDs, weight them by
     population, and summarise across the 31 members (median, p10, p90).
  5. Compare the 15-day total with the previous run and the run 24 h earlier.
  6. Draw the chart and the three maps.
  7. Write docs/latest.json, the PNGs, and the per-member CSV in history/.

Nothing touches disk (other than the GRIB cache) until every number and
picture has been computed, so a failure part-way through never leaves a
half-updated dashboard.

Usage:
  python hdd.py                      # newest complete run
  python hdd.py --run 2026-11-03T00  # a specific run (UTC cycle hour)
"""

import argparse
import io
import json
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

sys.stdout.reconfigure(encoding="utf-8")  # lets the degree sign print on Windows too

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
HISTORY = ROOT / "history"
CACHE = ROOT / "data"  # downloaded GRIB pieces land here (git-ignored)

# ----------------------------------------------------------------- locations
# HDD cities: name, latitude, longitude, rough population weight (millions).
# The weights are normalised to sum to 1 in city_hdd().
CITIES = [
    ("New York",      40.71,  -74.01, 8.3),
    ("Los Angeles",   34.05, -118.24, 3.9),
    ("Chicago",       41.88,  -87.63, 2.7),
    ("Houston",       29.76,  -95.37, 2.3),
    ("Philadelphia",  39.95,  -75.17, 1.6),
    ("Dallas",        32.78,  -96.80, 1.3),
    ("Columbus",      39.96,  -83.00, 0.9),
    ("Boston",        42.36,  -71.06, 0.7),
    ("Detroit",       42.33,  -83.05, 0.6),
    ("Denver",        39.74, -104.99, 0.7),
    ("Washington DC", 38.91,  -77.04, 0.7),
    ("Minneapolis",   44.98,  -93.27, 0.4),
    ("Pittsburgh",    40.44,  -79.99, 0.3),
    ("Kansas City",   39.10,  -94.58, 0.5),
    ("St. Louis",     38.63,  -90.20, 0.3),
    ("Atlanta",       33.75,  -84.39, 0.5),
]

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
SEARCH = ":TMP:2 m above ground:"  # the one GRIB field we want (2 m temperature)
LAT_MIN, LAT_MAX = 24.0, 50.0      # Lower-48 bounding box
LON_MIN, LON_MAX = -125.0, -66.0
THREADS = 16                       # parallel HTTP requests to S3
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


# ------------------------------------------------------------ 1. pick a run
def find_files(run):
    """Make one Herbie object per (member, lead time). Creating it checks that
    the GRIB file and its .idx index exist on S3. Returns the dict of objects
    and a list of the (member, lead) pairs that are missing."""

    def make(mf):
        m, f = mf
        H = Herbie(run, model="gefs", product="atmos.5", member=m, fxx=f,
                   priority="aws", save_dir=CACHE, verbose=False)
        return mf, H

    pairs = [(m, f) for m in MEMBERS for f in FXX]
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


def pick_run(explicit):
    """Return (run datetime, Herbie objects, note). With --run, that run must be
    complete. Otherwise try the newest cycle and step back 6 h at a time."""
    if explicit is not None:
        files, missing = find_files(explicit)
        if missing:
            m, f = missing[0]
            raise RuntimeError(f"run {explicit:%Y-%m-%d %H}Z is missing {len(missing)} "
                               f"of {len(files)} files on S3 (first: member {m}, f{f:03d})")
        return explicit, files, ""

    now = datetime.now(timezone.utc).replace(tzinfo=None)
    latest = now.replace(minute=0, second=0, microsecond=0)
    latest = latest.replace(hour=latest.hour - latest.hour % 6)  # 00/06/12/18Z
    note = ""
    for k in range(4):  # the newest cycle and the three before it
        run = latest - timedelta(hours=6 * k)
        files, missing = find_files(run)
        if not missing:
            print(f"  using run {run:%Y-%m-%d %H}Z (all {len(files)} files present)")
            return run, files, note
        print(f"  run {run:%Y-%m-%d %H}Z: {len(missing)} of {len(files)} files not on S3 yet")
        note += f"{run:%Y-%m-%d %H}Z was incomplete at {now:%H:%M} UTC. "
    raise RuntimeError("no complete GEFS run in the last 24 hours")


# ------------------------------------------------------------- 2. download
def download_all(files):
    """Fetch the 2 m temperature field for every file, 16 at a time. S3 sometimes
    answers "503 Slow Down"; back off and retry up to 5 times before giving up."""

    def fetch(H):
        for attempt in range(5):
            try:
                return H.download(SEARCH)
            except Exception:  # herbie wraps the HTTP error in a RuntimeError
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt + random.random())

    with ThreadPoolExecutor(THREADS) as pool:
        list(pool.map(fetch, files.values()))


# ------------------------------------------------------ 3. build the array
def load_cube(files):
    """Open every downloaded field and stack them into T[member, lead, lat, lon]
    (Kelvin), keeping only grid cells inside the Lower-48 box."""
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
            # Guard against a mixed-up file: the GRIB header must say this member and lead.
            if int(ds.number) != m or int(ds.step.values / np.timedelta64(1, "h")) != f:
                raise RuntimeError(f"file for member {m} f{f:03d} holds member "
                                   f"{int(ds.number)} step {ds.step.values}")
            T[mi, fi] = ds.t2m.values[np.ix_(ilat, ilon)]
            ds.close()
        print(f"  member {m:2d} loaded", flush=True)
    if np.isnan(T).any():
        raise RuntimeError("temperature array has gaps")
    if T.min() < 180 or T.max() > 340:  # Kelvin sanity check
        raise RuntimeError(f"temperatures out of range: {T.min():.0f}..{T.max():.0f} K")
    return T, lat, lon


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
    return days, [groups[d] for d in days]


def city_hdd(T, lat, lon, day_steps):
    """hdd[member, day] = sum over cities of weight x max(0, 65 - daily mean F).
    Daily mean = average of the four 6-hourly values in that UTC day."""
    weights = np.array([c[3] for c in CITIES])
    weights = weights / weights.sum()                       # normalise to 1
    rows, cols = zip(*[nearest_index(lat, lon, c[1], c[2]) for c in CITIES])
    temp_f = k_to_f(T[:, :, list(rows), list(cols)])        # [member, lead, city]
    hdd = np.zeros((len(MEMBERS), len(day_steps)))
    for di, steps in enumerate(day_steps):
        daily_mean = temp_f[:, steps, :].mean(axis=1)       # [member, city]
        per_city = np.maximum(0.0, HDD_BASE_F - daily_mean)  # [member, city]
        hdd[:, di] = per_city @ weights                     # weighted sum
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


# ------------------------------------------------- 5. compare with history
def history_runs():
    """Map run datetime -> CSV path for every file in history/."""
    return {datetime.strptime(p.stem, "hdd_%Y%m%d_%Hz"): p for p in HISTORY.glob("hdd_*.csv")}


def compare(hdd, days, other_run, csv_path):
    """Change in 15-day total versus another run, measured over the days both
    runs forecast (so a day dropping off the front of the window does not
    count as a 'change'). Totals are the median of each member's sum."""
    other = pd.read_csv(csv_path, index_col="date")
    shared = [d for d in days if d.isoformat() in other.index]
    if not shared:
        return None
    now_total = np.median(hdd[:, [days.index(d) for d in shared]].sum(axis=1))
    other_total = np.median(other.loc[[d.isoformat() for d in shared]].to_numpy().sum(axis=0))
    return {"run": f"{other_run:%Y-%m-%d %H}Z", "delta": round(float(now_total - other_total), 1),
            "shared_days": len(shared), "other_total": round(float(other_total), 1)}


# ------------------------------------------------------------- 6. pictures
def png_bytes(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def make_chart(run, days, hdd, p10, median, p90, total):
    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    for m in range(hdd.shape[0]):
        ax.plot(days, hdd[m], color="gray", alpha=0.25, lw=0.7)
    ax.fill_between(days, p10, p90, color="tab:blue", alpha=0.2, label="10th-90th percentile")
    ax.plot(days, median, color="tab:blue", lw=2.2, marker="o", ms=3.5, label="Ensemble median")
    ax.set_ylabel("Gas-weighted HDD (°F-days)")
    ax.set_ylim(bottom=0)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left")
    ax.set_title(f"US gas-weighted HDD - GEFS {run:%Y-%m-%d %H}Z, 31 members\n"
                 f"15-day total: median {total['median']:.0f}  "
                 f"(p10 {total['p10']:.0f} / p90 {total['p90']:.0f})")
    return png_bytes(fig)


def make_map(run, first, last, days, mean_f, uncertain, risk, lat, lon):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    plt.rcParams["hatch.color"] = "#333333"
    plt.rcParams["hatch.linewidth"] = 0.6
    pc = ccrs.PlateCarree()
    fig = plt.figure(figsize=(9, 5.8), dpi=130)
    ax = plt.axes(projection=pc)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=pc)
    levels = np.arange(-30, 101, 5)
    filled = ax.contourf(lon, lat, mean_f, levels=levels, cmap="RdYlBu_r", extend="both", transform=pc)
    ax.contour(lon, lat, mean_f, levels=[FREEZE_F], colors="black", linewidths=2, transform=pc)
    ax.contourf(lon, lat, uncertain.astype(float), levels=[0.5, 1.5], colors="none", hatches=["///"], transform=pc)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.4, edgecolor="#555555")
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.6, edgecolor="#333333")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6, edgecolor="#333333")
    for name, clat, clon, _ in CITIES:
        ax.plot(clon, clat, "o", ms=4, color="black", transform=pc)
    for name, short, blat, blon, threshold in BASINS:
        flagged = risk[name] >= 0.3
        ax.plot(blon, blat, "^", ms=12, color="red" if flagged else "limegreen",
                markeredgecolor="black", transform=pc, zorder=5)
        ax.text(blon + 0.5, blat + 0.5, f"{short} {risk[name]:.0%}", fontsize=7.5, weight="bold",
                transform=pc, zorder=6, bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    fig.colorbar(filled, ax=ax, orientation="horizontal", pad=0.03, shrink=0.8,
                 label="Ensemble-mean 2 m temperature (°F)")
    ax.set_title(f"GEFS {run:%Y-%m-%d %H}Z - days {first + 1}-{last} "
                 f"({days[first]:%b %d} to {days[last - 1]:%b %d})\n"
                 "black line = 32°F; hatched = <70% of members agree on below/above 32°F; "
                 "triangle red = freeze-off risk (% of members below basin threshold)", fontsize=9)
    return png_bytes(fig)


# ------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", type=lambda s: datetime.strptime(s, "%Y-%m-%dT%H"),
                        help="GEFS run to process, e.g. 2026-11-03T00 (UTC)")
    args = parser.parse_args()

    with step("1/7 find a complete GEFS run"):
        run, files, note = pick_run(args.run)

    with step("2/7 download 2 m temperature (31 members x 65 lead times)"):
        download_all(files)

    with step("3/7 read the GRIB files into one array"):
        T, lat, lon = load_cube(files)
        print(f"  array shape {T.shape} (member, lead, lat, lon)")

    with step("4/7 compute HDDs and map fields"):
        days, day_steps = full_days(run)
        hdd = city_hdd(T, lat, lon, day_steps)                      # [member, day]
        p10, median, p90 = np.percentile(hdd, [10, 50, 90], axis=0)
        member_totals = hdd.sum(axis=1)
        total = {k: round(float(v), 1) for k, v in
                 zip(["p10", "median", "p90"], np.percentile(member_totals, [10, 50, 90]))}
        maps = {label: period_fields(T, lat, lon, day_steps, a, b) for label, a, b in PERIODS}

    with step("5/7 compare with earlier runs"):
        runs = history_runs()
        earlier = [r for r in runs if r < run]
        prev_run = max(earlier) if earlier else None
        ago_run = run - timedelta(hours=24) if (run - timedelta(hours=24)) in runs else None
        vs_prev = compare(hdd, days, prev_run, runs[prev_run]) if prev_run else None
        vs_24h = compare(hdd, days, ago_run, runs[ago_run]) if ago_run else None

    with step("6/7 draw chart and maps"):
        pngs = {"hdd_forecast.png": make_chart(run, days, hdd, p10, median, p90, total)}
        for label, a, b in PERIODS:
            mean_f, uncertain, risk = maps[label]
            pngs[f"map_{label}.png"] = make_map(run, a, b, days, mean_f, uncertain, risk, lat, lon)

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

    with step("7/7 write outputs"):
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
        }
        for name, data in pngs.items():
            (DOCS / name).write_bytes(data)
        csv = pd.DataFrame(hdd.T, index=[d.isoformat() for d in days],
                           columns=[f"m{m:02d}" for m in MEMBERS]).rename_axis("date")
        csv.to_csv(HISTORY / f"hdd_{run_id}.csv", float_format="%.3f")
        (DOCS / "latest.json").write_text(json.dumps(latest, indent=1), encoding="utf-8")
        print(f"  wrote docs/latest.json, {len(pngs)} PNGs, history/hdd_{run_id}.csv")


if __name__ == "__main__":
    main()
