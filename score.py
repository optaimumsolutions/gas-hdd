#!/usr/bin/env python3
"""
score.py - score every past forecast in history/ against observed weather.

For every run CSV in history/ and every forecast date in it that now has a
row in history/observed.csv, one score row per index (us_hdd, <region>_hdd,
<region>_bcf): lead day, forecast median / p10 / p90, the observed value, the
error (median - observed) and whether the observed value fell inside p10-p90.

Outputs, all rebuilt from scratch each time so they always agree with history/:
  history/scores.csv      one row per run x date x index
  docs/accuracy.json      aggregates by lead day and by month, the last 14 days,
                          the CPC comparison and the Thursday storage checks
  docs/accuracy_lead.png  MAE and band coverage by lead day
  docs/accuracy_bands.png last 30 days: observed vs the day-1 / 7 / 14 bands

A statistic is only reported from MIN_N = 5 scores; below that it is null and
the page says "insufficient data". No network. Usage: python score.py
"""

import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from regions import REGIONS

ROOT = Path(__file__).resolve().parent
HISTORY = ROOT / "history"
DOCS = ROOT / "docs"
MIN_N = 5                 # fewer scores than this and a statistic is "insufficient data"
TARGET_COVERAGE = 0.80    # a p10-p90 band should contain the truth 80% of the time

# index name -> (history file prefix, observed.csv column)
INDEXES = {"us_hdd": ("hdd", "us_hdd")}
for _r in REGIONS:
    INDEXES[f"{_r['key']}_hdd"] = (f"{_r['key']}_hdd", f"{_r['key']}_hdd")
    INDEXES[f"{_r['key']}_bcf"] = (f"{_r['key']}_bcf", f"{_r['key']}_bcfd")
HDD_INDEXES = [k for k in INDEXES if k.endswith("_hdd")]   # same units, share one chart axis
SCORE_COLUMNS = ["run", "date", "lead_day", "index", "median", "p10", "p90",
                 "observed", "error", "abs_error", "inside_band"]


def load_observed(history=HISTORY):
    path = history / "observed.csv"
    return pd.read_csv(path, index_col="date") if path.exists() else pd.DataFrame()


def forecast_runs(prefix, history=HISTORY):
    """run datetime -> DataFrame[date rows, member columns] for every
    history/<prefix>_YYYYMMDD_HHz.csv (observed/scores/validation never match)."""
    runs = {}
    for path in sorted(history.glob(f"{prefix}_*.csv")):
        try:
            run = datetime.strptime(path.stem, f"{prefix}_%Y%m%d_%Hz")
        except ValueError:
            continue
        runs[run] = pd.read_csv(path, index_col="date")
    return runs


def build_scores(history=HISTORY):
    """One row per run x date x index for every forecast date that has an
    observation. lead_day 1 is the first forecast day of the run."""
    observed = load_observed(history)
    rows = []
    for index, (prefix, obs_col) in INDEXES.items():
        if obs_col not in observed.columns:
            continue
        obs = observed[obs_col].dropna()
        for run, frame in forecast_runs(prefix, history).items():
            members = frame.to_numpy(dtype=float)                       # [day, member]
            p10, median, p90 = np.percentile(members, [10, 50, 90], axis=1)
            for i, date in enumerate(frame.index):
                if date not in obs.index:
                    continue
                o, m = float(obs[date]), float(median[i])
                rows.append({"run": f"{run:%Y-%m-%d %H}Z", "date": date, "lead_day": i + 1, "index": index,
                             "median": round(m, 3), "p10": round(float(p10[i]), 3), "p90": round(float(p90[i]), 3),
                             "observed": round(o, 3), "error": round(m - o, 3), "abs_error": round(abs(m - o), 3),
                             "inside_band": bool(p10[i] <= o <= p90[i])})
    scores = pd.DataFrame(rows, columns=SCORE_COLUMNS)
    return scores.sort_values(["index", "date", "run"]).reset_index(drop=True)


# ---------------------------------------------------------------- aggregates
def stats(group):
    """n, MAE, bias (mean of median - observed), RMSE and band coverage for a
    set of score rows. The statistics are None below MIN_N rows."""
    n = int(len(group))
    if n < MIN_N:
        return {"n": n, "mae": None, "bias": None, "rmse": None, "coverage": None}
    err = group["error"].to_numpy(dtype=float)
    return {"n": n, "mae": round(float(np.abs(err).mean()), 2), "bias": round(float(err.mean()), 2),
            "rmse": round(float(np.sqrt((err ** 2).mean())), 2),
            "coverage": round(float(group["inside_band"].mean()), 3)}


def _num(x):
    """A JSON-safe float, or None for a blank/NaN cell."""
    return None if x is None or pd.isna(x) else round(float(x), 2)


def last_days(scores, observed, today, n_days=14, leads=(1, 5, 10)):
    """For each observed date in the last n_days: the observed values and what
    the latest run with lead day 1 / 5 / 10 said (its median)."""
    if observed.empty:
        return []
    start = (today - timedelta(days=n_days)).isoformat()
    rows = []
    for date in sorted(d for d in observed.index if d >= start):
        row = {"date": date, "observed": {index: _num(observed.loc[date, col])
                                          for index, (_, col) in INDEXES.items() if col in observed.columns}}
        for k in leads:
            picks = scores[(scores["date"] == date) & (scores["lead_day"] == k)].sort_values("run")
            row[f"day{k}"] = {index: (_num(g.iloc[-1]["median"]) if len(g) else None)
                              for index, g in ((i, picks[picks["index"] == i]) for i in INDEXES)}
        rows.append(row)
    return rows


def cpc_comparison(observed):
    """Our observed Midwest HDD against CPC's East/West North Central gas-weighted
    HDD: a check on the city weights, not on the forecast."""
    if observed.empty or {"midwest_hdd", "cpc_midwest_hdd"} - set(observed.columns):
        return {"n": 0, "mae": None, "bias": None}
    both = observed[["midwest_hdd", "cpc_midwest_hdd"]].dropna()
    err = (both["midwest_hdd"] - both["cpc_midwest_hdd"]).to_numpy(dtype=float)
    n = int(len(both))
    if n < MIN_N:
        return {"n": n, "mae": None, "bias": None}
    return {"n": n, "mae": round(float(np.abs(err).mean()), 2), "bias": round(float(err.mean()), 2)}


def storage_checks(history=HISTORY):
    """Rows logged by validate_log.py (the Thursday EIA check) and their hit rate."""
    path = history / "validation.csv"
    if not path.exists():
        return {"n": 0, "hit_rate": None, "rows": []}
    v = pd.read_csv(path)
    rows = [{k: (bool(x) if k == "inside_range" else (_num(x) if isinstance(x, (int, float, np.number)) else x))
             for k, x in r.items()} for r in v.to_dict("records")]
    n = int(len(v))
    hit = round(float(v["inside_range"].astype(bool).mean()), 2) if n >= MIN_N else None
    return {"n": n, "hit_rate": hit, "rows": rows}


def aggregate(scores, observed, history=HISTORY, now=None):
    """Everything docs/accuracy.json holds."""
    now = now or datetime.now(timezone.utc)
    by_lead, by_month, overall = {}, {}, {}
    for index in INDEXES:
        s = scores[scores["index"] == index]
        by_lead[index] = [dict(lead_day=k, **stats(s[s["lead_day"] == k])) for k in range(1, 16)]
        by_month[index] = [dict(month=m, **stats(g)) for m, g in s.groupby(s["date"].str[:7])]
        overall[index] = stats(s)
    return {
        "generated_utc": now.strftime("%Y-%m-%d %H:%M UTC"),
        "min_n": MIN_N,
        "target_coverage": TARGET_COVERAGE,
        "runs_scored": int(scores["run"].nunique()) if len(scores) else 0,
        "dates_observed": int(len(observed)),
        "observed_through": max(observed.index) if len(observed) else None,
        "indexes": list(INDEXES),
        "overall": overall,
        "by_lead": by_lead,
        "by_month": by_month,
        "last_14_days": last_days(scores, observed, now.date()),
        "cpc": cpc_comparison(observed),
        "storage": storage_checks(history),
        "images": {"lead": "accuracy_lead.png", "bands": "accuracy_bands.png"},
        "sources": {
            "observations": "Iowa Environmental Mesonet ASOS hourly temperatures, max/min per UTC day "
                            "at each city's airport; daily mean = (Tmax + Tmin) / 2",
            "cpc": "NOAA CPC daily gas-weighted HDD, East + West North Central, consumption-weighted average",
            "storage": "history/validation.csv, logged by validate_log.py from the Thursday EIA report",
        },
    }


# -------------------------------------------------------------------- charts
def png(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def draw_lead_chart(by_lead, path):
    """Left: MAE by lead day for the HDD indexes. Right: band coverage by lead
    day for every index, with the 80% target line."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), dpi=150)
    drew = False
    for index in INDEXES:
        pts = [(r["lead_day"], r["mae"], r["coverage"]) for r in by_lead[index] if r["mae"] is not None]
        if not pts:
            continue
        drew = True
        leads, mae, cov = zip(*pts)
        if index in HDD_INDEXES:
            ax1.plot(leads, mae, marker="o", label=index)
        ax2.plot(leads, cov, marker="o", label=index)
    ax2.axhline(TARGET_COVERAGE, ls="--", color="gray", label=f"target {TARGET_COVERAGE:.0%}")
    ax1.set_title("Mean absolute error of the median forecast")
    ax1.set_ylabel("MAE (HDD)")
    ax2.set_title("Share of days inside the p10-p90 band")
    ax2.set_ylim(0, 1)
    for ax in (ax1, ax2):
        ax.set_xlabel("lead day")
        ax.set_xlim(0.5, 15.5)
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="best", fontsize=8)
    if not drew:
        ax1.text(0.5, 0.5, f"insufficient data\n(fewer than {MIN_N} scores per lead day)",
                 ha="center", va="center", transform=ax1.transAxes, color="gray")
    png(fig, path)


def draw_bands_chart(scores, observed, today, path, n_days=30):
    """One panel per HDD index: observed values as black dots, and the p10-p90
    band of the latest run at lead day 1 (darkest), 7 and 14 (lightest)."""
    fig, axes = plt.subplots(len(HDD_INDEXES), 1, figsize=(10, 3.4 * len(HDD_INDEXES)), dpi=150, squeeze=False)
    start = (today - timedelta(days=n_days)).isoformat()
    for ax, index in zip(axes[:, 0], HDD_INDEXES):
        col = INDEXES[index][1]
        if len(observed) and col in observed.columns:
            o = observed[col].dropna()
            o = o[o.index >= start]
            if len(o):
                ax.plot(pd.to_datetime(o.index), o.values, "o", color="black", label="observed", zorder=5)
        s = scores[(scores["index"] == index) & (scores["date"] >= start)]
        for k, alpha in [(14, 0.12), (7, 0.25), (1, 0.45)]:
            g = s[s["lead_day"] == k].sort_values("run").groupby("date").last()   # latest run per date
            if len(g):
                ax.fill_between(pd.to_datetime(g.index), g["p10"], g["p90"], color="tab:blue",
                                alpha=alpha, lw=0, label=f"day-{k} p10-p90")
        ax.set_title(f"{index}: observed vs forecast bands, last {n_days} days")
        ax.set_ylabel("HDD")
        ax.set_ylim(bottom=0)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.grid(alpha=0.3)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="upper left", fontsize=8)
        else:
            ax.text(0.5, 0.5, "no scored days yet", ha="center", va="center", transform=ax.transAxes, color="gray")
    png(fig, path)


def main():
    scores = build_scores()
    observed = load_observed()
    now = datetime.now(timezone.utc)
    summary = aggregate(scores, observed, now=now)
    DOCS.mkdir(exist_ok=True)
    scores.to_csv(HISTORY / "scores.csv", index=False)
    draw_lead_chart(summary["by_lead"], DOCS / "accuracy_lead.png")
    draw_bands_chart(scores, observed, now.date(), DOCS / "accuracy_bands.png")
    (DOCS / "accuracy.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(f"scored {len(scores)} forecast-days from {summary['runs_scored']} runs "
          f"against {summary['dates_observed']} observed dates (through {summary['observed_through']})")
    for index in INDEXES:
        o = summary["overall"][index]
        print(f"  {index:14s} n={o['n']:3d}  MAE={o['mae']}  bias={o['bias']}  coverage={o['coverage']}")
    print("wrote history/scores.csv, docs/accuracy.json, docs/accuracy_lead.png, docs/accuracy_bands.png")


if __name__ == "__main__":
    main()
