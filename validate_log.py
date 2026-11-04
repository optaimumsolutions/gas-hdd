#!/usr/bin/env python3
"""Log one Thursday check: forecast implied Midwest withdrawal vs EIA's report.

  python validate_log.py --week-ending 2026-11-06 --eia-midwest -42

--eia-midwest is EIA's reported Midwest net change in Bcf, typed in by hand from
the Thursday storage report (negative = withdrawal). The forecast is the archived
run (Neon) nearest to 00Z on the Saturday the storage week starts.
"""

import argparse
import csv
from datetime import datetime, timedelta
from pathlib import Path

import db
from hdd import index_members, region_cities
from regions import MIDWEST, demand_bcfd, implied_withdrawal

HISTORY = Path(__file__).resolve().parent / "history"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--week-ending", required=True, type=lambda s: datetime.strptime(s, "%Y-%m-%d"))
parser.add_argument("--eia-midwest", required=True, type=float, help="EIA net change, Bcf (withdrawal < 0)")
args = parser.parse_args()

saturday = args.week_ending - timedelta(days=6)
with db.connect() as conn:
    run = min(db.list_runs(conn), key=lambda r: abs(r - saturday))
    loaded = db.load_forecasts(conn, run=run)[run]
D = demand_bcfd(MIDWEST, index_members(loaded["stations"], region_cities(MIDWEST)))   # [member, day]
week = [w for w in implied_withdrawal(MIDWEST, D, loaded["dates"])
        if w["week_ending"] == f"{args.week_ending:%Y-%m-%d}"]
if not week:
    raise SystemExit(f"run {run:%Y-%m-%d %H}Z does not cover the week ending {args.week_ending:%Y-%m-%d}")
w = week[0]
row = {"week_ending": w["week_ending"], "run_used": f"{run:%Y-%m-%d %H}Z", "days_in_forecast": w["days_in_forecast"],
       "implied_low": w["low_bcf"], "implied_high": w["high_bcf"], "eia_reported": args.eia_midwest,
       "inside_range": w["low_bcf"] <= -args.eia_midwest <= w["high_bcf"]}
out = HISTORY / "validation.csv"
new = not out.exists()
with out.open("a", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(row))
    if new:
        writer.writeheader()
    writer.writerow(row)
print(row)
