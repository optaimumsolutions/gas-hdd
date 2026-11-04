#!/usr/bin/env python3
"""
calibrate.py - fit the EMOS correction per airport and lead day (PRD-calibration.md, M5).

For every (airport, lead day) the archived runs give samples of the ensemble
mean m and spread s of the daily mean temperature, and the thermometer's
value y. EMOS (nonhomogeneous Gaussian regression) fits

    y ~ N(a + b*m,  c^2 + d^2*s^2)

by minimising the mean CRPS over the training window. The fit is applied to
the 31 members as a shift and scale (member -> mu + (member - m) * sigma / s),
so members keep their rank structure across airports and the index and Bcf
maths downstream stay member by member (PRD D16).

Gate (per airport and lead): the fit must lower the CRPS on the holdout
window; otherwise it is stored as the identity (0, 1, 0, 1) and listed.
Gate (overall, reported, not enforced): holdout p10-p90 coverage of the US and
Midwest HDD indexes within 0.75-0.85 at every lead and corrected MAE <= raw.

Writes gefs.emos_fit (Neon), docs/calibration.json and docs/calibration.png.

Usage:
  python calibrate.py                                   # train 2024-25, hold out 2025-26
  python calibrate.py --train-to 2026-11-01             # M7: train on everything, holdout still reported
"""

import argparse
import json
from datetime import date, datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

import db
import hdd
import score

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
IDENTITY = (0.0, 1.0, 0.0, 1.0)
MIN_SAMPLES = 30          # fewer training samples than this -> identity
COVERAGE_BAND = (0.75, 0.85)


# ------------------------------------------------------------------ the fit
def crps_gaussian(mu, sigma, y):
    z = (y - mu) / sigma
    return sigma * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1 / np.sqrt(np.pi))


def emos_sigma(c, d, s):
    return np.sqrt(c * c + d * d * s * s) + 1e-6


def fit_emos(m, s, y):
    """(a, b, c, d) minimising mean CRPS. Starts from the identity with c set
    to the residual spread."""
    def loss(p):
        a, b, c, d = p
        return crps_gaussian(a + b * m, emos_sigma(c, d, s), y).mean()

    x0 = [0.0, 1.0, float(np.std(y - m)), 1.0]
    res = minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 5000})
    a, b, c, d = res.x
    return float(a), float(b), float(abs(c)), float(abs(d))


def adjust(members, fit):
    """Shift and scale members[day, member] (daily mean F) with one fit."""
    a, b, c, d = fit
    m = members.mean(axis=1, keepdims=True)
    s = members.std(axis=1, keepdims=True)
    mu = a + b * m
    sigma = emos_sigma(c, d, s)
    return mu + (members - m) * np.where(s > 0, sigma / np.maximum(s, 1e-9), 0.0)


def adjust_run(stations, fits):
    """Apply the per-airport, per-lead fits to one run's station data.
    Returns station data whose Tmax and Tmin both equal the corrected daily
    mean, so hdd.index_members sees exactly the corrected means."""
    out = {}
    for st, (tmax, tmin) in stations.items():
        mean = (tmax + tmin) / 2                                   # [day, member]
        adj = mean.copy()
        for day in range(mean.shape[0]):
            fit = fits.get((st, day + 1))
            if fit is not None:
                adj[day] = adjust(mean[day:day + 1], fit)[0]
        out[st] = (adj, adj)
    return out


# -------------------------------------------------------------- the samples
def sample_table(forecasts, observed_means):
    """One row per run x date x airport: lead_day, m, s, y."""
    rows = []
    for run, fc in forecasts.items():
        for i, day in enumerate(fc["dates"]):
            obs = observed_means.get(day)
            if not obs:
                continue
            for st, (tmax, tmin) in fc["stations"].items():
                if st not in obs:
                    continue
                mean = (tmax[i] + tmin[i]) / 2
                rows.append((run, day, i + 1, st, float(mean.mean()), float(mean.std()), float(obs[st])))
    return pd.DataFrame(rows, columns=["run", "date", "lead_day", "station", "m", "s", "y"])


def fit_all(table, train_from, train_to, holdout):
    """Fit every (airport, lead); gate each on the holdout CRPS. Returns
    ({(station, lead): fit}, rows for gefs.emos_fit, list of identity keys)."""
    train = table[(table["date"] >= train_from) & (table["date"] <= train_to)]
    fits, rows, identity = {}, [], []
    for (st, lead), g in train.groupby(["station", "lead_day"]):
        h = holdout[(holdout["station"] == st) & (holdout["lead_day"] == lead)]
        fit, n = IDENTITY, len(g)
        if n >= MIN_SAMPLES:
            cand = fit_emos(g["m"].to_numpy(), g["s"].to_numpy(), g["y"].to_numpy())
            if len(h) < MIN_SAMPLES:
                fit = cand                                           # nothing to gate on
            else:
                hm, hs, hy = h["m"].to_numpy(), h["s"].to_numpy(), h["y"].to_numpy()
                raw = crps_gaussian(hm, np.maximum(hs, 1e-3), hy).mean()
                new = crps_gaussian(cand[0] + cand[1] * hm, emos_sigma(cand[2], cand[3], hs), hy).mean()
                fit = cand if new < raw else IDENTITY
        if fit == IDENTITY:
            identity.append(f"{st} d{lead}")
        fits[(st, lead)] = fit
        rows.append((st, lead, *fit, int(n), fit == IDENTITY, train_from, train_to))
    return fits, rows, identity


# ----------------------------------------------------------- index holdout
def index_holdout(forecasts, observed, fits, holdout_from, holdout_to):
    """Per index and lead: raw vs corrected MAE and p10-p90 coverage over the
    holdout runs. observed = score.observed_indexes()."""
    out = {index: {k: [] for k in range(1, 16)} for index in score.HDD_INDEXES}
    for run, fc in forecasts.items():
        if not (holdout_from <= run.date() <= holdout_to):
            continue
        corrected = adjust_run(fc["stations"], fits)
        for index in score.HDD_INDEXES:
            raw = score.index_value(fc["stations"], index)
            cal = score.index_value(corrected, index)
            for i, day in enumerate(fc["dates"]):
                key = day.isoformat()
                if key not in observed.index:
                    continue
                o = float(observed.loc[key, index])
                r10, r50, r90 = np.percentile(raw[:, i], [10, 50, 90])
                c10, c50, c90 = np.percentile(cal[:, i], [10, 50, 90])
                out[index][i + 1].append((abs(r50 - o), r10 <= o <= r90, abs(c50 - o), c10 <= o <= c90))
    result = {}
    for index, leads in out.items():
        result[index] = []
        for k, rows in leads.items():
            if not rows:
                result[index].append({"lead_day": k, "n": 0})
                continue
            arr = np.array(rows, dtype=float)
            result[index].append({"lead_day": k, "n": len(rows),
                                  "raw_mae": round(float(arr[:, 0].mean()), 2), "raw_coverage": round(float(arr[:, 1].mean()), 3),
                                  "cal_mae": round(float(arr[:, 2].mean()), 2), "cal_coverage": round(float(arr[:, 3].mean()), 3)})
    return result


def gate(holdout):
    """The overall gate from the PRD: coverage in band and MAE not worse, at every lead."""
    reasons = []
    for index, rows in holdout.items():
        for r in rows:
            if r["n"] < score.MIN_N:
                continue
            if not (COVERAGE_BAND[0] <= r["cal_coverage"] <= COVERAGE_BAND[1]):
                reasons.append(f"{index} day {r['lead_day']}: coverage {r['cal_coverage']:.2f}")
            if r["cal_mae"] > r["raw_mae"]:
                reasons.append(f"{index} day {r['lead_day']}: MAE {r['cal_mae']} > raw {r['raw_mae']}")
    return {"passed": not reasons, "failures": reasons}


def draw(holdout, path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), dpi=150)
    for index, rows in holdout.items():
        rows = [r for r in rows if r["n"] >= score.MIN_N]
        if not rows:
            continue
        leads = [r["lead_day"] for r in rows]
        line, = ax1.plot(leads, [r["raw_mae"] for r in rows], ls="--", alpha=0.6, label=f"{index} raw")
        ax1.plot(leads, [r["cal_mae"] for r in rows], marker="o", color=line.get_color(), label=f"{index} corrected")
        line, = ax2.plot(leads, [r["raw_coverage"] for r in rows], ls="--", alpha=0.6, label=f"{index} raw")
        ax2.plot(leads, [r["cal_coverage"] for r in rows], marker="o", color=line.get_color(), label=f"{index} corrected")
    ax2.axhspan(*COVERAGE_BAND, color="green", alpha=0.08, label="target 75-85%")
    ax1.set_title("Holdout MAE of the median, raw vs corrected")
    ax1.set_ylabel("MAE (HDD)")
    ax2.set_title("Holdout p10-p90 coverage, raw vs corrected")
    ax2.set_ylim(0, 1)
    for ax in (ax1, ax2):
        ax.set_xlabel("lead day")
        ax.set_xlim(0.5, 15.5)
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=7)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--train-from", type=date.fromisoformat, default=date(2024, 11, 1))
    parser.add_argument("--train-to", type=date.fromisoformat, default=date(2025, 3, 31))
    parser.add_argument("--holdout-from", type=date.fromisoformat, default=date(2025, 11, 1))
    parser.add_argument("--holdout-to", type=date.fromisoformat, default=date(2026, 3, 31))
    args = parser.parse_args()

    with db.connect() as conn:
        forecasts = db.load_forecasts(conn)
        observed_means = db.load_observed_means(conn)
        table = sample_table(forecasts, observed_means)
        holdout = table[(table["date"] >= args.holdout_from) & (table["date"] <= args.holdout_to)]
        print(f"{len(table)} samples from {len(forecasts)} runs; training {args.train_from}..{args.train_to}, "
              f"holdout {args.holdout_from}..{args.holdout_to}")
        fits, rows, identity = fit_all(table, args.train_from, args.train_to, holdout)
        print(f"fitted {len(rows)} (airport, lead) pairs; {len(identity)} kept as identity")
        observed = score.observed_indexes(observed_means, {})
        result = index_holdout(forecasts, observed, fits, args.holdout_from, args.holdout_to)
        verdict = gate(result)
        db.write_emos(conn, rows)

    now = datetime.now(timezone.utc)
    summary = {"fitted_at": now.strftime("%Y-%m-%d %H:%M UTC"), "method": "EMOS / NGR per airport per lead, CRPS-minimised, "
               "applied to members as shift-and-scale (ECC)",
               "train_from": args.train_from.isoformat(), "train_to": args.train_to.isoformat(),
               "holdout_from": args.holdout_from.isoformat(), "holdout_to": args.holdout_to.isoformat(),
               "fits": len(rows), "identity": identity, "holdout": result, "gate": verdict,
               "image": "calibration.png"}
    DOCS.mkdir(exist_ok=True)
    (DOCS / "calibration.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    draw(result, DOCS / "calibration.png")
    for index, leads in result.items():
        print(f"\n{index}: lead  n   raw MAE  cal MAE  raw cov  cal cov")
        for r in leads:
            if r["n"]:
                print(f"  {r['lead_day']:>4} {r['n']:>4} {r['raw_mae']:>8} {r['cal_mae']:>8} {r['raw_coverage']:>8} {r['cal_coverage']:>8}")
    print(f"\ngate {'PASSED' if verdict['passed'] else 'FAILED'}" + ("" if verdict["passed"] else ": " + "; ".join(verdict["failures"][:8])))
    print("wrote gefs.emos_fit, docs/calibration.json, docs/calibration.png")


if __name__ == "__main__":
    main()
