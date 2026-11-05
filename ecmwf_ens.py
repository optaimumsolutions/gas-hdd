#!/usr/bin/env python3
"""
ecmwf_ens.py - archive ECMWF's open-data ensemble (IFS ENS, 50 perturbed members, 0.25°) at
the 22 HDD airports, the way hdd.py archives GEFS (PRD-calibration.md, M8).

  1. Pick the newest complete 00Z or 12Z run on the public bucket (the 360 h
     file of the run must exist), or the run given with --run.
  2. For every 6-hourly step 0..360 h, fetch only the 2 m temperature messages
     (50 members, ~16 MB per step) with ECMWF's own open-data client, which
     coalesces the byte ranges into one request per step and honours the
     servers' Retry-After (the raw S3 mirror answers 503 SlowDown to naive
     range requests). Sources tried in order: ecmwf, azure, aws. The 22 airport
     grid cells are read and the step's file deleted: a run never sits on disk.
  3. Daily Tmax / Tmin per airport and member = max / min of the 6-hourly 2 m
     temperatures valid 00..24Z of each UTC day (the open-data ensemble has no
     6-hour max/min fields; the per-airport EMOS fit absorbs the difference
     from GEFS's definition). First 15 full days.
  4. Write gefs.forecast_station with model = 'ecmwf' (upsert).

No charts, no docs/: ECMWF is the second line on hdd.py's charts, read back
from Neon. Scoring and calibration: score.py --model ecmwf,
calibrate.py --model ecmwf.

Usage:
  python ecmwf_ens.py                        # newest complete 00Z/12Z run
  python ecmwf_ens.py --run 2026-11-01T00    # one run (backfill)
"""

import argparse
import os
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import cfgrib
import numpy as np
from ecmwf.opendata import Client

import db
import hdd
from hdd import STATIONS, k_to_f, retry, step

sys.stdout.reconfigure(encoding="utf-8")
MODEL = "ecmwf"
CACHE = hdd.CACHE / "ecmwf"
STEPS = list(range(0, 361, 6))             # 6-hourly, 0..360 h (61 steps, 15 days)
CYCLES = (0, 12)                           # the open-data ensemble runs to 360 h at 00Z and 12Z only
THREADS = int(os.environ.get("THREADS", 2))   # the servers rate-limit; keep this low
SOURCES = os.environ.get("ECMWF_SOURCES", "ecmwf,azure,aws").split(",")


def client(source):
    return Client(source=source, model="ifs", resol="0p25", maximum_retries=12, retry_after=20)


def complete(run):
    """A run is complete when its last step can be listed (the client's
    latest() / the index for step 360 exist)."""
    for source in SOURCES:
        try:
            c = client(source)
            latest = c.latest(stream="enfo", type="cf", param="2t", step=STEPS[-1])
            return run <= latest.replace(tzinfo=None)
        except Exception as err:                                    # noqa: BLE001
            print(f"  {source}: cannot check completeness ({str(err)[:80]})")
    raise RuntimeError("no ECMWF source reachable")


def pick_run(explicit):
    if explicit is not None:
        if explicit.hour not in CYCLES:
            raise RuntimeError(f"ECMWF ensemble runs to 360 h only at {CYCLES}Z; got {explicit:%H}Z")
        return explicit                                             # the fetch itself fails if it is missing
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
    array [member, station] in K, perturbed members 1..50 (the open-data
    stream carries no control member for 2 m temperature).
    The step's GRIB subset is deleted after."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{run:%Y%m%d%H}_{fxx:03d}.grib2"
    last = None
    for source in SOURCES:
        try:
            client(source).retrieve(date=run.strftime("%Y%m%d"), time=run.hour, step=fxx, stream="enfo",
                                    type="pf", param="2t", target=str(path))
            last = None
            break
        except Exception as err:                                    # noqa: BLE001
            last = err
            print(f"  step {fxx:3d}: {source} failed ({str(err)[:80]}), trying the next source", flush=True)
    if last is not None:
        raise RuntimeError(f"step {fxx} could not be fetched from any source: {last}")
    blocks = []
    for ds in cfgrib.open_datasets(str(path)):                      # normally one dataset (pf)
        t2m = ds["t2m"].values
        if t2m.ndim == 2:
            t2m = t2m[None]
        cells = station_cells(ds)
        number = np.atleast_1d(ds["number"].values) if "number" in ds.coords else np.array([0])
        pts = np.array([[t2m[m, r, c] for r, c in cells] for m in range(t2m.shape[0])], dtype=np.float32)
        blocks.append((number, pts))
        ds.close()
    numbers = np.concatenate([n for n, _ in blocks])
    pts = np.concatenate([p for _, p in blocks])
    out = pts[np.argsort(numbers)]                                  # control (0) first, then 1..50
    path.unlink(missing_ok=True)
    for extra in CACHE.glob(path.name + "*.idx"):
        extra.unlink(missing_ok=True)
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

    with step(f"2/4 fetch 2 m temperature at {len(STATIONS)} airports, {len(needed)} steps x 50 members"):
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
