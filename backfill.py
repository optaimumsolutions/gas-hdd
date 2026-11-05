#!/usr/bin/env python3
"""
backfill.py - archive past GEFS 00Z runs in Neon (PRD-calibration.md, M2).

Three small commands, used by .github/workflows/backfill.yml:

  plan    --from 2024-11-01 --to 2025-03-31 --shards 20
          Print a JSON list of shards. Each shard is a space-separated string
          of runs ("2024-11-01T00 2024-11-02T00 ...") not yet in Neon.
  run     2024-11-01T00 2024-11-02T00 ...
          Archive each run with `hdd.py --run <run> --no-docs`, deleting the
          GRIB cache after each one. Keeps going after a failure; exits 1 at
          the end if any run failed.
  missing --from 2024-11-01 --to 2025-03-31
          Print the runs in the range that are still not in Neon (exit 1 if any).

Only 00Z runs are backfilled (one per day). Nothing is written to docs/ or
history/, and nothing is committed. --model ecmwf does the same for the ECMWF
ensemble (ecmwf_ens.py --run).
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import db

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "data"


def runs_between(start, end):
    """Every 00Z run from start to end inclusive, as naive UTC datetimes."""
    d = start
    while d <= end:
        yield datetime(d.year, d.month, d.day)
        d += timedelta(days=1)


def archived_runs(conn, model):
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT run FROM gefs.forecast_station WHERE model = %s", (model,))
        return {r[0].astimezone(timezone.utc).replace(tzinfo=None) for r in cur.fetchall()}


def missing_runs(start, end, model):
    with db.connect() as conn:
        done = archived_runs(conn, model)
    return [r for r in runs_between(start, end) if r not in done]


def cmd_plan(args):
    todo = missing_runs(args.start, args.end, args.model)
    shards = [todo[i::args.shards] for i in range(args.shards)]
    shards = [" ".join(f"{r:%Y-%m-%dT%H}" for r in s) for s in shards if s]
    print(json.dumps(shards))
    print(f"{len(todo)} runs to do in {len(shards)} shards", file=sys.stderr)


def cmd_run(args):
    failed = []
    for i, run in enumerate(args.runs, 1):
        print(f"\n===== run {i}/{len(args.runs)}: {run} =====", flush=True)
        if args.model == "ecmwf":
            cmd = [sys.executable, str(ROOT / "ecmwf_ens.py"), "--run", run]
        else:
            cmd = [sys.executable, str(ROOT / "hdd.py"), "--run", run, "--no-docs"]
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"===== {run} FAILED (exit {result.returncode})", flush=True)
            failed.append(run)
        shutil.rmtree(CACHE, ignore_errors=True)   # ~300 MB per run; the runner's disk is finite
    print(f"\n{len(args.runs) - len(failed)} of {len(args.runs)} runs archived"
          + (f"; failed: {' '.join(failed)}" if failed else ""))
    sys.exit(1 if failed else 0)


def cmd_missing(args):
    todo = missing_runs(args.start, args.end, args.model)
    for r in todo:
        print(f"{r:%Y-%m-%dT%H}")
    print(f"{len(todo)} run(s) missing between {args.start} and {args.end}", file=sys.stderr)
    sys.exit(1 if todo else 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", default="gefs", choices=["gefs", "ecmwf"])
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name, fn in [("plan", cmd_plan), ("missing", cmd_missing)]:
        p = sub.add_parser(name)
        p.add_argument("--from", dest="start", type=date.fromisoformat, required=True)
        p.add_argument("--to", dest="end", type=date.fromisoformat, required=True)
        if name == "plan":
            p.add_argument("--shards", type=int, default=20)
        p.set_defaults(fn=fn)
    p = sub.add_parser("run")
    p.add_argument("runs", nargs="+", help="runs as YYYY-MM-DDTHH")
    p.set_defaults(fn=cmd_run)
    args = parser.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
