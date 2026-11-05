#!/usr/bin/env python3
"""
fit_demand.py - what an HDD is worth in gas, fitted from data (PRD-calibration.md, M11).

Consumption: EIA monthly residential + commercial natural gas consumption per
state (free API, DEMO_KEY), summed over each region's states and turned into
Bcf/d, against the month's mean observed daily HDD for that region (the same
airport thermometers and weights the forecast uses, from gefs.observed, which
holds 10 years once observe.py --climatology has run).
Storage: EIA weekly Midwest working gas -> weekly change, against the storage
week's summed observed Midwest HDD.

Fits (ordinary least squares on heating months / winter weeks of the last
FIT_WINTERS winters; every winter in the archive shown as the range):
    Bcf/d  = intercept + slope * HDD/day            (per region)
    weekly storage change (Bcf) = c + d * weekly HDD (Midwest)

Writes docs/demand.json (the two reference tables per region and the fit
numbers) and prints the regions.py constants to paste in by hand. No Neon
writes. Usage: python fit_demand.py
"""

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import db
import hdd
from regions import REGIONS

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
CACHE = ROOT / "data" / "eia"
EIA = "https://api.eia.gov/v2"
EIA_KEY = "DEMO_KEY"                       # a handful of requests; no account needed
FIT_WINTERS = 5                            # slope from the last 5 winters, range from all
MIN_DAYS = 25                              # a month needs this many observed days
HEATING_MIN_HDD = 5.0                      # months / weeks colder than this go into the fit
LEVELS = list(range(0, 55, 5))             # HDD/day rows of table A
DEVIATIONS = [-15, -10, -5, -2, 0, 2, 5, 10, 15]   # HDD/day rows of table B
STORAGE_SERIES = {"midwest": "NW2_EPG0_SWO_R33_BCF"}
STATES = {"us": ["US"]}
for _r in REGIONS:
    STATES[_r["key"]] = _r["states"]


# ----------------------------------------------------------------- EIA data
def eia(path, params, series):
    """Rows for the given EIA series, cached per series in data/eia/ (a week).
    Missing ones are fetched five at a time; DEMO_KEY rate limits (429) are
    waited out."""
    CACHE.mkdir(parents=True, exist_ok=True)
    rows, missing = [], []
    for sid in series:
        f = CACHE / f"{sid}.json"
        if f.exists() and (time.time() - f.stat().st_mtime) < 7 * 86400:
            rows += json.loads(f.read_text(encoding="utf-8"))
        else:
            missing.append(sid)
    for i in range(0, len(missing), 5):
        chunk = missing[i:i + 5]
        url = (f"{EIA}/{path}/data/?" + urllib.parse.urlencode(dict(params, api_key=EIA_KEY))
               + "".join(f"&facets[series][]={v}" for v in chunk))
        for wait in (0, 15, 45, 90):
            time.sleep(wait)
            try:
                data = json.load(urllib.request.urlopen(url, timeout=60))["response"]["data"]
                break
            except urllib.error.HTTPError as err:
                if err.code != 429 or wait == 90:
                    raise
                print("  EIA rate limit, waiting", flush=True)
        for sid in chunk:
            (CACHE / f"{sid}.json").write_text(json.dumps([r for r in data if r["series"] == sid]), encoding="utf-8")
        rows += data
        time.sleep(2)
    return rows


def monthly_bcfd(states, start="2014-01"):
    """Series[month 'YYYY-MM'] of residential + commercial Bcf/d over states."""
    series = [f"N30{k}{st}2" for st in states for k in ("10", "20")]
    rows = eia("natural-gas/cons/sum", {"frequency": "monthly", "data[0]": "value", "start": start, "length": 5000}, series)
    f = pd.DataFrame(rows)
    f["value"] = pd.to_numeric(f["value"], errors="coerce")
    mmcf = f.groupby("period")["value"].sum(min_count=len(series))      # every series present or NaN
    days = pd.Series({p: pd.Period(p).days_in_month for p in mmcf.index})
    return (mmcf / days / 1000).dropna().rename("bcfd")


def weekly_storage(series, start="2014-01-01"):
    """Series[week-ending Friday] of weekly change in working gas, Bcf (negative = withdrawal)."""
    rows = eia("natural-gas/stor/wkly", {"frequency": "weekly", "data[0]": "value", "start": start, "length": 5000}, [series])
    f = pd.DataFrame(rows)
    s = pd.Series(pd.to_numeric(f["value"]).to_numpy(), index=pd.to_datetime(f["period"])).sort_index()
    return s.diff().dropna().rename("change_bcf")


# ------------------------------------------------------------ observed HDD
def observed_hdd(conn):
    """DataFrame[date -> us, <region>...] daily observed weighted HDD, only days
    where every station of that index reported."""
    means = db.load_observed_means(conn)
    cities = {"us": hdd.national_cities()}
    for r in REGIONS:
        cities[r["key"]] = hdd.region_cities(r)
    rows = {}
    for day, m in means.items():
        row = {}
        for key, cl in cities.items():
            if all(st in m for st, _ in cl):
                row[key] = float(hdd.weighted_hdd(np.array([m[st] for st, _ in cl]), [w for _, w in cl]))
        if row:
            rows[pd.Timestamp(day)] = row
    return pd.DataFrame.from_dict(rows, orient="index").sort_index()


def winter_of(ts):
    """Winter label: Nov-Dec belong to the winter starting that year, Jan-Apr to the one before."""
    return ts.year if ts.month >= 7 else ts.year - 1


def ols(x, y):
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    return float(a), float(b), float(resid.std(ddof=2)), float(1 - resid.var() / y.var())


# ------------------------------------------------------------------- fits
def fit_region(key, daily, bcfd):
    """Monthly consumption fit for one index + the two reference tables."""
    d = daily[key].dropna()
    month = d.groupby(d.index.to_period("M")).agg(["mean", "count"])
    month = month[month["count"] >= MIN_DAYS]["mean"].rename("hdd")
    month.index = month.index.astype(str)
    both = pd.concat([month, bcfd], axis=1, join="inner").dropna()
    both["winter"] = [winter_of(pd.Period(p).to_timestamp()) for p in both.index]
    heating = both[both["hdd"] >= HEATING_MIN_HDD]
    winters = sorted(heating["winter"].unique())
    recent = heating[heating["winter"].isin(winters[-FIT_WINTERS:])]
    a, b, sd, r2 = ols(recent["hdd"].to_numpy(), recent["bcfd"].to_numpy())
    a_all, b_all, _, _ = ols(heating["hdd"].to_numpy(), heating["bcfd"].to_numpy())
    floor = float(both[both["hdd"] < 1]["bcfd"].median()) if (both["hdd"] < 1).any() else None

    level_rows = []
    for lv in LEVELS:
        near = heating[(heating["hdd"] - lv).abs() <= 2.5] if lv >= HEATING_MIN_HDD else both[(both["hdd"] - lv).abs() <= 2.5]
        level_rows.append({"hdd_per_day": lv, "fitted_bcfd": round(max(a + b * lv, floor or 0), 2),
                           "seen_n_months": int(len(near)),
                           "seen_min": round(float(near["bcfd"].min()), 2) if len(near) else None,
                           "seen_max": round(float(near["bcfd"].max()), 2) if len(near) else None,
                           "seen_months": [f"{p} ({h:.0f} HDD, {v:.1f})" for p, h, v in
                                           near.sort_values("hdd")[["hdd", "bcfd"]].tail(4).itertuples()]})
    return {
        "slope_bcfd_per_hdd": round(b, 4), "intercept_bcfd": round(a, 3), "resid_sd_bcfd": round(sd, 3), "r2": round(r2, 3),
        "floor_bcfd": None if floor is None else round(floor, 2),
        "fit_winters": [int(w) for w in winters[-FIT_WINTERS:]], "fit_months": int(len(recent)),
        "all_winters": {"from": int(winters[0]), "to": int(winters[-1]), "slope": round(b_all, 4), "intercept": round(a_all, 3),
                        "months": int(len(heating))},
        "months": [{"month": p, "hdd": round(h, 1), "bcfd": round(v, 2), "winter": int(w)}
                   for p, h, v, w in both[["hdd", "bcfd", "winter"]].itertuples()],
        "level_table": level_rows,
    }


def fit_storage(daily, change):
    """Weekly Midwest storage change vs the storage week's summed HDD (weeks end Friday)."""
    d = daily["midwest"].dropna()
    wk = d.resample("W-FRI").agg(["sum", "count"])
    wk = wk[wk["count"] == 7]["sum"].rename("hdd_week")
    both = pd.concat([wk, change], axis=1, join="inner").dropna()
    both["winter"] = [winter_of(t) for t in both.index]
    cold = both[both["hdd_week"] >= 7 * HEATING_MIN_HDD]
    winters = sorted(cold["winter"].unique())
    recent = cold[cold["winter"].isin(winters[-FIT_WINTERS:])]
    c, dd, sd, r2 = ols(recent["hdd_week"].to_numpy(), recent["change_bcf"].to_numpy())
    return {"bcf_per_weekly_hdd": round(dd, 4), "intercept_bcf": round(c, 2), "resid_sd_bcf": round(sd, 2), "r2": round(r2, 3),
            "fit_winters": [int(w) for w in winters[-FIT_WINTERS:]], "fit_weeks": int(len(recent)),
            "weeks": [{"week_ending": t.date().isoformat(), "hdd_week": round(h, 1), "change_bcf": round(v, 1)}
                      for t, h, v in both[["hdd_week", "change_bcf"]].itertuples()]}


def deviation_table(fit, storage, normal_hdd):
    """Table B: HDD/day away from normal -> Bcf/d, Bcf per 7-day week, and
    (Midwest) storage change per week. normal_hdd = typical Nov-Mar HDD/day."""
    rows = []
    for dev in DEVIATIONS:
        row = {"deviation_hdd_per_day": dev, "bcfd_change": round(fit["slope_bcfd_per_hdd"] * dev, 2),
               "bcf_per_week": round(7 * fit["slope_bcfd_per_hdd"] * dev, 1),
               "demand_at_normal_bcfd": round(fit["intercept_bcfd"] + fit["slope_bcfd_per_hdd"] * (normal_hdd + dev), 2)}
        if storage:
            row["storage_change_per_week_bcf"] = round(7 * storage["bcf_per_weekly_hdd"] * dev, 1)
        rows.append(row)
    return rows


def main():
    with db.connect() as conn:
        daily = observed_hdd(conn)
        climo = db.load_climatology(conn)
    print(f"observed daily HDD: {len(daily)} days, {daily.index.min().date()} to {daily.index.max().date()}")
    out = {"generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
           "method": f"OLS on heating months (mean HDD >= {HEATING_MIN_HDD}) of the last {FIT_WINTERS} winters; "
                     "EIA monthly residential + commercial consumption per state vs observed airport HDD",
           "regions": {}}
    storage_fit = None
    for key, states in STATES.items():
        bcfd = monthly_bcfd(states)
        fit = fit_region(key, daily, bcfd)
        # typical Nov-Mar HDD/day from the station normals, with this index's weights
        cities = hdd.national_cities() if key == "us" else hdd.region_cities(next(r for r in REGIONS if r["key"] == key))
        normals = np.array([climo[st] for st, _ in cities])                       # [city, doy]
        winter_doys = list(range(305, 366)) + list(range(1, 91))                  # Nov 1 .. Mar 31
        normal_hdd = float(np.mean([hdd.weighted_hdd(normals[:, d - 1], [w for _, w in cities]) for d in winter_doys]))
        if key in STORAGE_SERIES:
            storage_fit = fit_storage(daily, weekly_storage(STORAGE_SERIES[key]))
            fit["storage"] = storage_fit
        fit["normal_winter_hdd_per_day"] = round(normal_hdd, 1)
        fit["deviation_table"] = deviation_table(fit, fit.get("storage"), normal_hdd)
        out["regions"][key] = fit
        print(f"\n{key}: Bcf/d = {fit['intercept_bcfd']:.2f} + {fit['slope_bcfd_per_hdd']:.3f} x HDD/day  "
              f"(r2 {fit['r2']}, resid sd {fit['resid_sd_bcfd']} Bcf/d, {fit['fit_months']} months, winters {fit['fit_winters']}); "
              f"all {fit['all_winters']['from']}-{fit['all_winters']['to']}: slope {fit['all_winters']['slope']:.3f}; "
              f"summer floor {fit['floor_bcfd']} Bcf/d; normal winter HDD {normal_hdd:.1f}/day")
        if "storage" in fit:
            s = fit["storage"]
            print(f"  storage: weekly change = {s['intercept_bcf']:.1f} + {s['bcf_per_weekly_hdd']:.3f} x weekly HDD "
                  f"(r2 {s['r2']}, sd {s['resid_sd_bcf']} Bcf, {s['fit_weeks']} weeks) -> per HDD/day over a week: "
                  f"{7 * s['bcf_per_weekly_hdd']:.2f} Bcf; storage/consumption ratio {7 * s['bcf_per_weekly_hdd'] / (7 * fit['slope_bcfd_per_hdd']):.2f}")
    DOCS.mkdir(exist_ok=True)
    (DOCS / "demand.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("\nwrote docs/demand.json")


if __name__ == "__main__":
    main()
