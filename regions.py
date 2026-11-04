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

Model source for every slope/intercept/ratio below: EIA, "STEO Supplement:
Weather Sensitivity in Natural Gas Markets" (Oct 2014).
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
# Slope: EIA STEO Supplement Table 2, res+comm Bcf/d per HDD/d by census
# division: East North Central 0.33, West North Central 0.12, East South
# Central 0.09. Scaled to the storage Midwest by each division's share inside
# it: ENC minus Ohio (0.33 x 0.80) + WNC's IA/MN/MO (0.12 x 0.70) + KY/TN
# (0.09 x 0.60) = 0.40. A 2010-14 fit, not fitted to this index: placeholder
# until the Phase 3 re-fit.
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
    "slope_bcfd_per_hdd": 0.40,     # EIA Table 2: ENC 0.33x0.80 + WNC 0.12x0.70 + ESC 0.09x0.60
    "intercept_bcfd": 0.0,          # EIA national res (-2.05) + comm (+2.28) intercepts net to ~0
    "floor_bcfd": 1.9,              # summer res/comm baseload, scaled from the 13-state 2.3
    "cold_kink_hdd": 45,            # daily mean below 20 F: demand steepens (EIA extreme-cold note)
    "cold_kink_extra_frac": 0.15,   # +15% slope above the kink; placeholder until re-fit
    "storage_ratio": (1.12, 1.26),  # EIA storage slope 1.78-2.00 vs res+comm 1.59 Bcf/d per HDD
    "basins": ["Appalachia / Pittsburgh"],   # names from hdd.BASINS
    # The band width is no longer a placeholder: calibrate.py fits it per
    # airport and per lead from the forecast archive (PRD-calibration.md).
    # The slope and kink above are still EIA's 2010-14 numbers.
    "refit": True,
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


def implied_withdrawal(region, D, days):
    """Implied storage withdrawal per EIA storage week (Saturday-Friday) that
    overlaps the forecast. For each week: each member's demand summed over the
    week's days inside the forecast, the median of those member sums, times the
    storage_ratio low and high. D[member, day] in Bcf/d; days = list of dates."""
    weeks = {}
    for i, d in enumerate(days):
        week_ending = d + timedelta(days=(4 - d.weekday()) % 7)   # the Friday
        weeks.setdefault(week_ending, []).append(i)
    lo, hi = region["storage_ratio"]
    rows = []
    for week_ending, idx in sorted(weeks.items()):
        median = float(np.median(D[:, idx].sum(axis=1)))
        rows.append({"week_ending": week_ending.isoformat(), "days_in_forecast": len(idx),
                     "low_bcf": round(lo * median, 1), "high_bcf": round(hi * median, 1)})
    return rows
