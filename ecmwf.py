#!/usr/bin/env python3
"""
ecmwf.py - archive ECMWF's open-data ensemble (IFS ENS, 51 members, 0.25°) at
the 22 HDD airports, the way hdd.py archives GEFS (PRD-calibration.md, M8).

  1. Pick the newest complete 00Z or 12Z run on the public bucket (the 360 h
     file of the run must exist), or the run given with --run.
  2. For every 6-hourly step 0..360 h, fetch only the 2 m temperature messages
     (51 members, ~100 MB per step via byte ranges), read the 22 airport grid
     cells, and delete the step's file: the whole run never sits on disk.
  3. Daily Tmax / Tmin per airport and member = max / min of the 6-hourly 2 m
     temperatures valid 00..24Z of each UTC day (the open-data ensemble has no
     6-hour max/min fields; the per-airport EMOS fit absorbs the difference
     from GEFS's definition). First 15 full days.
  4. Write gefs.forecast_station with model = 'ecmwf' (upsert).

No charts, no docs/: ECMWF is the second line on hdd.py's charts, read back
from Neon. Scoring and calibration: score.py --model ecmwf,
calibrate.py --model ecmwf.

Usage:
  python ecmwf.py                        # newest complete 00Z/12Z run
  python ecmwf.py --run 2026-11-01T00    # one run (backfill)
"""

import argparse
import os
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import numpy as np
from herbie import Herbie

import db
import hdd
from hdd import STATIONS, k_to_f, retry, step

sys.stdout.reconfigure(encoding="utf-8")
MODEL = "ecmwf"
CACHE = hdd.CACHE / "ecmwf"
STEPS = list(range(0, 361, 6))             # 6-hourly, 0..360 h (61 steps, 15 days)
CYCLES = (0, 12)                           # the open-data ensemble runs to 360 h at 00Z and 12Z only
THREADS = int(os.environ.get("THREADS", 3))   # the bucket answers "SlowDown" quickly; keep this low
SEARCH = ":2t:"


def herbie(run, fxx):
    return retry(lambda: Herbie(run, model="ecmwf", product="enfo", resolution="0p25", fxx=fxx,
                                priority="aws", save_dir=CACHE, verbose=False))


def complete(run):
    """A run is complete when its last step file and index exist."""
    H = herbie(run, STEPS[-1])
    return H.grib is not None and H.idx is not None


def pick_run(explicit):
    if explicit is not None:
        if explicit.hour not in CYCLES:
            raise RuntimeError(f"ECMWF ensemble runs to 360 h only at {CYCLES}Z; got {explicit:%H}Z")
        if not complete(explicit):
            raise RuntimeError(f"run {explicit:%Y-%m-%d %H}Z is not (yet) complete on the bucket")
        return explicit
    now = datetime.now(timezone.utc).replace(tzinfo=None, minute=0, second=0, microsecond=0)
    latest = now.replace(hour=12 if now.hour >= 12 else 0)
    for k in range(4):
        run = latest - timedelta(hours=12 * k)
        if complete(run):
            print(f"  using run {run:%Y-%m-%d %H}Z")
            return run
        print(f"  run {run:%Y-%m-%d %H}Z not complete yet")
    raise RuntimeError("no complete ECMWF run in the last 48 hours")


def station_cells(ds):
    """Row/column of each airport on this grid, in STATIONS order."""
    lat = ds.latitude.values
    lon = (ds.longitude.values + 180) % 360 - 180
    return [hdd.nearest_index(lat, lon, la, lo) for la, lo in STATIONS.values()]


def fetch_step(run, fxx):
    """2 m temperature at the 22 airports for all members at one step:
    array [member, station] in K. The step's GRIB subset is deleted after."""
    H = herbie(run, fxx)
    retry(lambda: H.download(SEARCH))
    ds = H.xarray(SEARCH, remove_grib=True)
    t2m = ds["t2m"].values                                         # [member, lat, lon]
    if t2m.ndim == 2:                                               # a single member would drop the axis
        t2m = t2m[None]
    cells = station_cells(ds)
    out = np.array([[t2m[m, r, c] for r, c in cells] for m in range(t2m.shape[0])], dtype=np.float32)
    ds.close()
    return fxx, out


def full_days(run, steps):
    """[(date, [steps valid in 00..24Z of that date])] for the first 15 days
    that have both boundaries (00Z and 24Z) inside the run."""
    by_day = {}
    for s in steps:
        valid = run + timedelta(hours=s)
        by_day.setdefault(valid.date(), []).append(s)
        if valid.hour == 0 and s > 0:                               # 00Z also closes the previous day
            by_day.setdefault((valid - timedelta(days=1)).date(), []).append(s)
    days = [d for d in sorted(by_day) if (run + timedelta(hours=min(by_day[d]))).hour == 0
            and (run + timedelta(hours=max(by_day[d]))) == datetime.combine(d, datetime.min.time()) + timedelta(days=1)]
    days = days[:15]
    if len(days) < 14:                                              # 00Z runs give 15 days, 12Z runs 14
        raise RuntimeError(f"expected 14 or 15 full days, found {len(days)}")
    return [(d, sorted(by_day[d])) for d in days]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", type=lambda s: datetime.strptime(s, "%Y-%m-%dT%H"), help="run, e.g. 2026-11-01T00")
    args = parser.parse_args()

    with step("1/4 find a complete ECMWF run"):
        run = pick_run(args.run)
        days = full_days(run, STEPS)
        needed = sorted({s for _, ss in days for s in ss})

    with step(f"2/4 fetch 2 m temperature at {len(STATIONS)} airports, {len(needed)} steps x 51 members"):
        t0 = time.time()
        with ThreadPoolExecutor(THREADS) as pool:
            temps = dict(pool.map(lambda f: fetch_step(run, f), needed))   # step -> [member, station] K
        shutil.rmtree(CACHE, ignore_errors=True)
        members = next(iter(temps.values())).shape[0]
        print(f"  {len(needed)} steps in {time.time() - t0:.0f}s, {members} members")

    with step("3/4 daily Tmax / Tmin per airport"):
        tmax = np.zeros((members, len(days), len(STATIONS)))
        tmin = np.zeros_like(tmax)
        for di, (_, ss) in enumerate(days):
            block = np.stack([temps[s] for s in ss])                   # [step, member, station]
            tmax[:, di] = k_to_f(block.max(axis=0))
            tmin[:, di] = k_to_f(block.min(axis=0))
        if np.nanmin(tmin) < -80 or np.nanmax(tmax) > 130:
            raise RuntimeError("temperatures out of range")
        dates = [d for d, _ in days]
        print(f"  days {dates[0]} .. {dates[-1]}")

    with step("4/4 archive in Neon"):
        with retry(db.connect) as conn:
            n = db.write_forecast(conn, run.replace(tzinfo=timezone.utc), dates, STATIONS, tmax, tmin, model=MODEL)
        print(f"  wrote {n} rows (model {MODEL}, {members} members) to gefs.forecast_station")


if __name__ == "__main__":
    main()
