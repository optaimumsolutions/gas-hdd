#!/usr/bin/env python3
"""Log one Thursday check: forecast implied Midwest withdrawal vs EIA's report.

  python validate_log.py --week-ending 2026-11-06 --eia-midwest -42

--eia-midwest is EIA's reported Midwest net change in Bcf, typed in by hand from
the Thursday storage report (negative = withdrawal). The forecast is the stored
run nearest to 00Z on the Saturday the storage week starts. No network.
"""

import argparse
import csv
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd

from regions import MIDWEST, implied_withdrawal

HISTORY = Path(__file__).resolve().parent / "history"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--week-ending", required=True, type=lambda s: datetime.strptime(s, "%Y-%m-%d"))
parser.add_argument("--eia-midwest", required=True, type=float, help="EIA net change, Bcf (withdrawal < 0)")
args = parser.parse_args()

saturday = args.week_ending - timedelta(days=6)
runs = {datetime.strptime(p.stem, "midwest_bcf_%Y%m%d_%Hz"): p for p in HISTORY.glob("midwest_bcf_*.csv")}
run = min(runs, key=lambda r: abs(r - saturday))
D = pd.read_csv(runs[run], index_col="date")
days = [datetime.strptime(d, "%Y-%m-%d").date() for d in D.index]
week = [w for w in implied_withdrawal(MIDWEST, D.to_numpy().T, days)
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
