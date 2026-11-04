#!/usr/bin/env python3
"""
observe.py - observed daily weather for the HDD cities, appended to
history/observed.csv so score.py can grade every past forecast.

Truth source: Iowa Environmental Mesonet (IEM) hourly ASOS temperatures for
each city's airport (free, no key), pulled in UTC. For every station and UTC
calendar day the max and min of the hourly values give
daily mean = (Tmax + Tmin) / 2, the definition the forecast uses, and the
shared functions in hdd.py / regions.py turn those means into the national
HDD, each region's HDD and Bcf/d. Nothing is reimplemented here.

Why hourly and not IEM's daily summaries: the daily product is a local
calendar day, the forecast is a UTC day; the two differ by 5-8 hours, which
matters on frontal days. Hourly data can be cut at 00Z exactly.

Reference, not truth: NOAA CPC's observed gas-weighted HDD for the East and
West North Central census divisions, consumption-weighted into one Midwest
number and stored alongside as cpc_midwest_hdd. It checks the city weights.

Rules for history/observed.csv:
  - fetch the last --days days (default 20); add only dates not already present
  - never overwrite a row: a refetch that differs by more than 0.5 HDD is
    logged and the first value kept (a blank cpc_midwest_hdd cell is filled in
    once CPC publishes that date; that is the one exception)
  - a UTC day is used only once it is over and every station reported at
    least MIN_HOURS hourly values for it; otherwise it is skipped and logged,
    and tried again next run

If IEM is unreachable the script prints why and exits 1; the forecast
workflow treats that as a warning, never as a failed forecast.

Every complete station-day also goes to Neon (gefs.observed, insert-only; see
PRD-calibration.md M3). --from/--to backfill a date range in 31-day chunks.
--climatology rebuilds gefs.climatology: per airport, 10 years of IEM hourly
data -> daily mean -> mean and sd per day of year, 15-day smoothed.

Usage: python observe.py [--days 20 | --from 2024-11-01 --to 2025-04-15] [--climatology]
"""

import argparse
import io
import sys
import time
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

import db
import hdd
from regions import REGIONS, demand_bcfd

ROOT = Path(__file__).resolve().parent
HISTORY = ROOT / "history"
OBSERVED = HISTORY / "observed.csv"
IEM_URL = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
CPC_URL = "https://ftp.cpc.ncep.noaa.gov/htdocs/degree_days/weighted/daily_data/{year}/UtilityGas.Heating.txt"
MIN_HOURS = 18             # hourly reports a station needs in a UTC day for the day to count
REVISION_TOLERANCE = 0.5   # HDD; a refetch that differs by more than this is logged, not applied
SOURCE = "iem_asos_hourly_utc"
CHUNK_DAYS = 31            # one IEM request covers at most this many days (all stations)
CLIMO_YEARS = 10           # climatology window
PAUSE = 2                  # seconds between IEM requests
# CPC census divisions that overlap the storage Midwest, keyed by CPC's row label
CPC_DIVISIONS = {"3": ["IL", "IN", "MI", "OH", "WI"],              # East North Central
                 "4": ["IA", "KS", "MN", "MO", "NE", "ND", "SD"]}  # West North Central

# (city, station, weight) lists straight from the forecast's own tables
NATIONAL = [(c[0], c[4], c[3]) for c in hdd.CITIES]
REGIONAL = {r["key"]: [(c[0], c[5], c[4]) for c in r["cities"]] for r in REGIONS}
STATIONS = sorted({s for _, s, _ in NATIONAL} | {s for cities in REGIONAL.values() for _, s, _ in cities})


# ------------------------------------------------------------------ fetching
def fetch_iem(stations, start, end):
    """Hourly temperatures (F) in UTC for all stations in one request, days
    start..end. Retried with backoff like the GEFS downloads; IEM answers
    "Too many requests" with a text body, which counts as a failure."""
    params = [("station", s) for s in stations] + [
        ("data", "tmpf"), ("tz", "Etc/UTC"), ("format", "onlycomma"), ("latlon", "no"),
        ("missing", "M"), ("trace", "T"), ("direct", "no"), ("report_type", "3"),
        ("year1", start.year), ("month1", start.month), ("day1", start.day),
        ("year2", end.year), ("month2", end.month), ("day2", end.day)]

    def get():
        r = requests.get(IEM_URL, params=params, timeout=180)
        r.raise_for_status()
        if not r.text.startswith("station,valid,tmpf"):
            raise RuntimeError(f"IEM answered: {r.text[:80]!r}")
        return r.text

    frame = pd.read_csv(io.StringIO(hdd.retry(get)), na_values=["M"])
    frame["valid"] = pd.to_datetime(frame["valid"])
    return frame.dropna(subset=["tmpf"])


def fetch_iem_range(stations, start, end):
    """fetch_iem over start..end inclusive, CHUNK_DAYS at a time with a pause
    between requests (IEM asks for that on long pulls). IEM's end day is
    exclusive, so each request asks for one day past the chunk."""
    frames = []
    chunk_start = start
    while chunk_start <= end:
        chunk_end = min(chunk_start + timedelta(days=CHUNK_DAYS - 1), end)
        frames.append(fetch_iem(stations, chunk_start, chunk_end + timedelta(days=1)))
        print(f"  IEM {chunk_start} to {chunk_end}: {len(frames[-1])} hourly rows", flush=True)
        chunk_start = chunk_end + timedelta(days=1)
        if chunk_start <= end:
            time.sleep(PAUSE)
    return pd.concat(frames, ignore_index=True)


def daily_extremes(hourly):
    """Per station and UTC date: max, min and how many hourly values there were."""
    hourly = hourly.assign(date=hourly["valid"].dt.strftime("%Y-%m-%d"))
    return hourly.groupby(["station", "date"])["tmpf"].agg(tmax="max", tmin="min", hours="count").reset_index()


def complete_days(extremes, stations, last_day):
    """{date: {station: daily mean F}} for dates <= last_day where every station
    has at least MIN_HOURS values, plus a list of (date, reason) for the rest."""
    means, dropped = {}, []
    for date, g in extremes.groupby("date"):
        if date > last_day:
            continue
        g = g.set_index("station")
        short = [s for s in stations if s not in g.index or g.loc[s, "hours"] < MIN_HOURS]
        if short:
            dropped.append((date, f"{len(short)} station(s) with fewer than {MIN_HOURS} hourly reports: {' '.join(short)}"))
            continue
        means[date] = {s: float(g.loc[s, "tmax"] + g.loc[s, "tmin"]) / 2 for s in stations}
    return means, dropped


def station_day_rows(extremes, last_day):
    """Rows for gefs.observed: every station-day <= last_day with at least
    MIN_HOURS hourly reports, judged per station (one short station does not
    discard the others' days there, unlike the index rows)."""
    ok = extremes[(extremes["date"] <= last_day) & (extremes["hours"] >= MIN_HOURS)]
    return [(r.station, date.fromisoformat(r.date), float(r.tmax), float(r.tmin), int(r.hours), SOURCE)
            for r in ok.itertuples(index=False)]


# -------------------------------------------------------------- climatology
def climatology_rows(station, daily):
    """daily = DataFrame[date, tmax, tmin, hours] for one station over
    CLIMO_YEARS years. Mean and sd of the daily mean per day of year (Feb 29
    folded into day 60), each smoothed with a 15-day circular window."""
    d = daily[daily["hours"] >= MIN_HOURS].copy()
    when = pd.to_datetime(d["date"])
    d["doy"] = when.dt.dayofyear - (when.dt.is_leap_year & (when.dt.dayofyear > 60)).astype(int)
    d["tmean"] = (d["tmax"] + d["tmin"]) / 2
    by = d.groupby("doy")["tmean"].agg(["mean", "std", "count"]).reindex(range(1, 366)).interpolate()
    years = int(round(by["count"].mean()))
    smooth = {}
    for col in ("mean", "std"):
        x = by[col].to_numpy()
        pad = np.concatenate([x[-7:], x, x[:7]])             # circular: Dec 25 neighbours Jan 1
        smooth[col] = np.convolve(pad, np.ones(15) / 15, mode="valid")
    return [(station, doy, float(smooth["mean"][i]), float(smooth["std"][i]), years)
            for i, doy in enumerate(range(1, 366))]


def build_climatology(conn, stations, today):
    """One IEM pull per station, CLIMO_YEARS years back; writes gefs.climatology."""
    start = date(today.year - CLIMO_YEARS, today.month, today.day)
    for i, s in enumerate(stations):
        daily = daily_extremes(fetch_iem([s], start, today))   # end exclusive: through yesterday
        rows = climatology_rows(s, daily)
        db.write_climatology(conn, rows)
        jan, jul = rows[14], rows[195]
        print(f"  {s}: {len(daily)} days, {rows[0][4]} yr; Jan 15 {jan[2]:.1f}+-{jan[3]:.1f} F, "
              f"Jul 15 {jul[2]:.1f}+-{jul[3]:.1f} F", flush=True)
        if i < len(stations) - 1:
            time.sleep(PAUSE)


def cpc_weights():
    """Shares of East and West North Central inside the Midwest: the Midwest
    city weights (res+comm gas consumption) summed by census division."""
    midwest = next((r for r in REGIONS if r["key"] == "midwest"), None)
    if midwest is None:
        return {}
    w = {d: sum(c[4] for c in midwest["cities"] if c[1] in states) for d, states in CPC_DIVISIONS.items()}
    total = sum(w.values())
    return {d: v / total for d, v in w.items()}


def fetch_cpc(years):
    """{date: cpc_midwest_hdd} from CPC's daily gas-weighted HDD file(s).
    Layout: three header lines, then 'Region|YYYYMMDD|...' and one row per
    census division ('3' = East North Central, '4' = West North Central)."""
    weights = cpc_weights()
    out = {}
    for year in sorted(set(years)):
        def get():
            r = requests.get(CPC_URL.format(year=year), timeout=60)
            r.raise_for_status()
            return r.text

        lines = hdd.retry(get).splitlines()
        header = next(line for line in lines if line.startswith("Region|"))
        dates = [datetime.strptime(d, "%Y%m%d").date().isoformat() for d in header.split("|")[1:]]
        rows = {line.split("|")[0]: line.split("|")[1:] for line in lines if line.split("|")[0] in weights}
        for i, date in enumerate(dates):
            try:
                out[date] = round(sum(weights[d] * float(rows[d][i]) for d in weights), 2)
            except (KeyError, ValueError, IndexError):
                pass  # that date has no value yet
    return out


# ------------------------------------------------------------------ the row
def observed_row(means):
    """One observed.csv row from {station: daily mean F}. Per-station HDD, the
    national index, and each region's HDD and Bcf/d come from the same
    functions the forecast calls (hdd.hdd_from_mean, hdd.weighted_hdd,
    regions.demand_bcfd)."""
    row = {f"hdd_{s}": round(float(hdd.hdd_from_mean(np.array(means[s]))), 3) for s in STATIONS}
    national = np.array([means[s] for _, s, _ in NATIONAL])
    row["us_hdd"] = round(float(hdd.weighted_hdd(national, [w for _, _, w in NATIONAL])), 3)
    for r in REGIONS:
        cities = REGIONAL[r["key"]]
        w = hdd.weighted_hdd(np.array([means[s] for _, s, _ in cities]), [wt for _, _, wt in cities])
        row[f"{r['key']}_hdd"] = round(float(w), 3)
        row[f"{r['key']}_bcfd"] = round(float(demand_bcfd(r, w)), 3)
    return row


def merge(existing, new, log=print):
    """Append rows for dates not already present. Existing rows are never
    overwritten: a refetch that differs by more than REVISION_TOLERANCE HDD is
    logged and the first value kept. The one exception is a blank
    cpc_midwest_hdd cell, filled when CPC publishes that date."""
    if existing is None or existing.empty:
        return new.copy(), list(new.index)
    merged = existing.copy()
    added = []
    for date, row in new.iterrows():
        if date not in merged.index:
            merged = pd.concat([merged, new.loc[[date]]])
            added.append(date)
            continue
        for col in [c for c in merged.columns if c.endswith("_hdd") and not c.startswith("cpc")]:
            old, cur = merged.loc[date, col], row.get(col)
            if pd.notna(old) and pd.notna(cur) and abs(float(old) - float(cur)) > REVISION_TOLERANCE:
                log(f"  revision ignored: {date} {col} stored {float(old):.2f}, refetched {float(cur):.2f} (kept the first)")
        if ("cpc_midwest_hdd" in merged.columns and pd.isna(merged.loc[date, "cpc_midwest_hdd"])
                and pd.notna(row.get("cpc_midwest_hdd"))):
            merged.loc[date, "cpc_midwest_hdd"] = row["cpc_midwest_hdd"]
            log(f"  filled cpc_midwest_hdd for {date}")
    return merged.sort_index(), added


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--days", type=int, default=20, help="how many past days to fetch (default 20)")
    parser.add_argument("--from", dest="start", type=date.fromisoformat, help="backfill: first UTC day")
    parser.add_argument("--to", dest="end", type=date.fromisoformat, help="backfill: last UTC day")
    parser.add_argument("--climatology", action="store_true", help="rebuild gefs.climatology (10 years per airport)")
    args = parser.parse_args()
    now = datetime.now(timezone.utc)
    today = now.date()
    yesterday = today - timedelta(days=1)                # the newest complete UTC day
    start = args.start or today - timedelta(days=args.days)
    end = min(args.end or yesterday, yesterday)
    last_day = end.isoformat()

    if args.climatology:
        print(f"[observe] climatology: {len(STATIONS)} stations, {CLIMO_YEARS} years each")
        with db.connect() as conn:
            build_climatology(conn, STATIONS, today)

    print(f"[observe] IEM hourly temperatures, {len(STATIONS)} stations, {start} to {last_day} (UTC days)")
    try:
        hourly = fetch_iem_range(STATIONS, start, end)
    except Exception as err:
        print(f"  IEM unavailable, skipping observations this run: {err}")
        sys.exit(1)
    extremes = daily_extremes(hourly)
    means, dropped = complete_days(extremes, STATIONS, last_day)
    for date_, why in dropped:
        print(f"  skipped {date_}: {why}")

    try:
        with db.connect() as conn:
            new = db.write_observed(conn, station_day_rows(extremes, last_day))
        print(f"  Neon gefs.observed: {new} new station-day(s)")
    except Exception as err:
        print(f"  WARNING: Neon write skipped: {err}")

    try:
        cpc = fetch_cpc(range(start.year, today.year + 1))
    except Exception as err:
        print(f"  CPC unavailable ({err}); cpc_midwest_hdd left blank for new dates")
        cpc = {}

    rows = {}
    for day in sorted(means):
        row = observed_row(means[day])
        row["cpc_midwest_hdd"] = cpc.get(day)
        row["source"] = SOURCE
        row["fetched_utc"] = now.strftime("%Y-%m-%d %H:%M")
        rows[day] = row
    new = pd.DataFrame.from_dict(rows, orient="index").rename_axis("date")
    existing = pd.read_csv(OBSERVED, index_col="date") if OBSERVED.exists() else None
    merged, added = merge(existing, new)

    HISTORY.mkdir(exist_ok=True)
    merged.index.name = "date"
    merged.to_csv(OBSERVED, float_format="%.3f")
    print(f"  added {len(added)} date(s): {', '.join(added) if added else 'none'}; "
          f"observed.csv now holds {len(merged)} dates")
    show = [c for c in ["us_hdd"] + [f"{r['key']}_hdd" for r in REGIONS]
            + [f"{r['key']}_bcfd" for r in REGIONS] + ["cpc_midwest_hdd"] if c in merged.columns]
    print(merged[show].tail(10).to_string())


if __name__ == "__main__":
    main()
