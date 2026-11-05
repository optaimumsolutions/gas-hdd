"""
regions.py - regional definitions for the gas-hdd dashboard. Pure data plus the
HDD -> Bcf demand model. No I/O.

Regions are EIA's five weekly-storage regions (not Census), because that is the
only regional cut of the gas market with a free weekly observation: the Thursday
storage report. Phase 1 defines the Midwest only.

Each region:
  cities  (city, state, lat, lon, weight, station); weights sum to 1.0 within
          the region. Coordinates are approximate airport (ASOS) locations and
          station is that airport's ICAO code, which observe.py scores against.
  demand  Bcf/d = max(floor, intercept + slope * W + slope * kink_frac * max(0, W - kink))
  storage weekly EIA change (Bcf) = storage_intercept + storage_per_hdd * (week's summed W)

Slope, intercept, floor and the storage fit were fitted on 2026-11-05 by
fit_demand.py (EIA monthly residential + commercial consumption per state and
EIA weekly Midwest working gas, against this index's own observed airport HDD;
last 5 winters 2021-22..2025-26 for the numbers, 10 winters checked). Re-run
fit_demand.py and paste the printed numbers here to refit. The cold kink is
still EIA's 2014 STEO Supplement note.
"""

from datetime import timedelta

import numpy as np

# ------------------------------------------------------------------- Midwest
# EIA storage Midwest = IL, IN, IA, KY, MI, MN, MO, TN, WI.
# (Ohio is East; Kansas is South Central; NE and the Dakotas are Mountain.)
#
# City weights = each state's share of the 9-state residential + commercial gas
# consumption, 2024 + 2025 combined (EIA Natural Gas Consumption by End Use,
# NG_CONS_SUM_A_EPG0_VRS_MMCF_A and _VCS_, released 2026-10-30):
#
#   state   res+comm Bcf   share
#   IL         1,190.1     0.26375
#   MI           944.4     0.20930   -> Detroit 60% / Grand Rapids 40%
#   WI           490.8     0.10877
#   MN           465.4     0.10315
#   IN           455.4     0.10093
#   MO           317.9     0.07046   -> St. Louis 60% / Kansas City 40%
#   TN           262.0     0.05805   (estimate: 2024 123.6 + 2025 138.3)
#   IA           219.8     0.04872
#   KY           166.4     0.03687
#   total      4,512.2     1.00000
#
# The MI and MO 60/40 splits are by metro population, not utility customer
# counts: a judgment call, flagged for refinement.
#
# Slope 0.354 Bcf/d per HDD/day, intercept 0.96, floor 1.82: fit_demand.py
# 2026-11-05, r2 0.993, residual sd 0.33 Bcf/d over 37 heating months of the
# last 5 winters (all 10 winters: slope 0.349). EIA's 2014 census-division
# regrouping gave 0.40 / 0 / 1.9, the placeholders this replaced.
MIDWEST = {
    "name": "Midwest",
    "key": "midwest",
    "states": ["IL", "IN", "IA", "KY", "MI", "MN", "MO", "TN", "WI"],
    "eia_storage_region": "Midwest",
    "cities": [
        ("Chicago",      "IL", 41.98, -87.90, 0.26375, "ORD"),
        ("Detroit",      "MI", 42.21, -83.35, 0.12558, "DTW"),   # MI 60% (DTE territory)
        ("Grand Rapids", "MI", 42.88, -85.52, 0.08372, "GRR"),   # MI 40% (Consumers territory)
        ("Milwaukee",    "WI", 42.95, -87.90, 0.10877, "MKE"),
        ("Minneapolis",  "MN", 44.88, -93.22, 0.10315, "MSP"),
        ("Indianapolis", "IN", 39.72, -86.29, 0.10093, "IND"),
        ("St. Louis",    "MO", 38.75, -90.37, 0.04228, "STL"),   # MO 60%
        ("Kansas City",  "MO", 39.30, -94.71, 0.02818, "MCI"),   # MO 40%
        ("Nashville",    "TN", 36.12, -86.68, 0.05805, "BNA"),
        ("Des Moines",   "IA", 41.53, -93.66, 0.04872, "DSM"),
        ("Louisville",   "KY", 38.17, -85.73, 0.03687, "SDF"),
    ],
    "slope_bcfd_per_hdd": 0.354,    # fit_demand.py 2026-11-05 (EIA 2014 placeholder was 0.40)
    "intercept_bcfd": 0.96,         # fit_demand.py 2026-11-05 (was 0)
    "floor_bcfd": 1.82,             # median Bcf/d of months with < 1 HDD/day (was 1.9)
    "cold_kink_hdd": 45,            # daily mean below 20 F: demand steepens (EIA extreme-cold note)
    "cold_kink_extra_frac": 0.15,   # +15% slope above the kink; still EIA's note, not fitted
    # EIA weekly Midwest change (Bcf) = 72.7 - 0.471 x the storage week's summed HDD;
    # residual sd 16.9 Bcf, r2 0.84, 120 winter weeks of the last 5 winters (fit_demand.py 2026-11-05).
    # Replaces the national storage/consumption ratio 1.12-1.26 (the fitted ratio is 1.33).
    "storage_intercept_bcf": 72.7,
    "storage_per_weekly_hdd": -0.471,
    "storage_sd_bcf": 16.9,
    "basins": ["Appalachia / Pittsburgh"],   # names from hdd.BASINS
    "refit": True,                  # band (calibrate.py) and demand (fit_demand.py) both fitted; kink is not
}

REGIONS = [MIDWEST]   # Phase 2 appends EAST, SOUTH_CENTRAL, MOUNTAIN, PACIFIC

# Res+comm slopes (Bcf/d per HDD/d) per EIA storage region, from the same EIA
# Table 2 census-division slopes regrouped into storage regions. They sum to the
# EIA national res+comm slope of 1.59. Not used yet: Phase 2 copies each one
# into its region dict.
PHASE2_SLOPES = {"east": 0.62, "midwest": 0.40, "south_central": 0.23, "mountain": 0.14, "pacific": 0.20}

for _r in REGIONS:
    assert abs(sum(c[4] for c in _r["cities"]) - 1) < 1e-4, f"{_r['key']} weights do not sum to 1"


def demand_bcfd(region, W):
    """Weighted HDD -> res/comm gas demand in Bcf/d. Works on any array shape.
    Linear above the floor, with the slope 15% steeper beyond the cold kink.
    Nonlinear, so percentiles must be taken AFTER converting each member."""
    slope = region["slope_bcfd_per_hdd"]
    extra = slope * region["cold_kink_extra_frac"] * np.maximum(0.0, W - region["cold_kink_hdd"])
    return np.maximum(region["floor_bcfd"], region["intercept_bcfd"] + slope * W + extra)


def scenario_10pct(region, W_median_path):
    """Bcf change over the window if every day were 10% colder / warmer (in HDD)
    than the median path. Linear in the slope, as EIA quotes it."""
    delta = 0.10 * region["slope_bcfd_per_hdd"] * float(np.sum(W_median_path))
    return {"colder_bcf": round(delta, 1), "warmer_bcf": round(-delta, 1)}


Z90 = 1.2816  # standard-normal quantile for p90 (p10 = -Z90)


def implied_withdrawal(region, W, days):
    """Implied EIA weekly storage change per storage week (Saturday-Friday)
    that overlaps the forecast, from the fitted weekly relation
    change = storage_intercept_bcf + storage_per_weekly_hdd * weekly HDD.
    For each week: each member's HDD summed over the week's days inside the
    forecast, the median of those sums, through the fit; low/high = the median
    -/+ 1.2816 x the fit's residual sd (the p10-p90 of past weeks).
    EIA's sign: negative = withdrawal, positive = injection.
    W[member, day] = weighted HDD; days = list of dates. A week only partly
    inside the forecast is scaled to 7 days."""
    weeks = {}
    for i, d in enumerate(days):
        week_ending = d + timedelta(days=(4 - d.weekday()) % 7)   # the Friday
        weeks.setdefault(week_ending, []).append(i)
    rows = []
    for week_ending, idx in sorted(weeks.items()):
        hdd_week = float(np.median(W[:, idx].sum(axis=1))) * 7 / len(idx)
        change = region["storage_intercept_bcf"] + region["storage_per_weekly_hdd"] * hdd_week
        rows.append({"week_ending": week_ending.isoformat(), "days_in_forecast": len(idx),
                     "hdd_week": round(hdd_week, 1), "median_bcf": round(change, 1),
                     "low_bcf": round(change - Z90 * region["storage_sd_bcf"], 1),
                     "high_bcf": round(change + Z90 * region["storage_sd_bcf"], 1)})
    return rows
