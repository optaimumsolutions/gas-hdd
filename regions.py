"""
regions.py - regional definitions for the gas-hdd dashboard. Pure data plus the
HDD -> Bcf demand model. No I/O.

Regions are EIA's five weekly-storage regions (not Census), because that is the
only regional cut of the gas market with a free weekly observation: the Thursday
storage report. Phase 1 defines the Midwest only.

Each region:
  cities  (city, state, lat, lon, weight); weights sum to 1.0 within the region.
          Coordinates are approximate airport (ASOS) locations.
  demand  Bcf/d = max(floor, intercept + slope * weighted_HDD)
"""

import numpy as np

Z90 = 1.2816  # standard-normal quantile for p90 (p10 = -Z90)

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
# Slope: EIA STEO Supplement "Weather Sensitivity in Natural Gas Markets"
# (Oct 2014), res+comm Bcf/d per HDD/d by census division: East North Central
# 0.33, West North Central 0.12, East South Central 0.09. Scaled to the storage
# Midwest: ENC minus Ohio (0.33 x 0.80 = 0.26) + WNC's IA/MN/MO (0.12 x 0.65 =
# 0.08) + KY/TN (0.09 x 0.70 = 0.06) = 0.38-0.40. Plausible range 0.33-0.45.
# A 2010-14 fit, not fitted to this index: placeholder until the Phase 3 re-fit.
MIDWEST = {
    "name": "Midwest",
    "key": "midwest",
    "states": ["IL", "IN", "IA", "KY", "MI", "MN", "MO", "TN", "WI"],
    "eia_storage_region": "Midwest",
    "cities": [
        ("Chicago",      "IL", 41.98, -87.90, 0.26375),
        ("Detroit",      "MI", 42.21, -83.35, 0.12558),   # MI 60% (DTE territory)
        ("Grand Rapids", "MI", 42.88, -85.52, 0.08372),   # MI 40% (Consumers territory)
        ("Milwaukee",    "WI", 42.95, -87.90, 0.10877),
        ("Minneapolis",  "MN", 44.88, -93.22, 0.10315),
        ("Indianapolis", "IN", 39.72, -86.29, 0.10093),
        ("St. Louis",    "MO", 38.75, -90.37, 0.04228),   # MO 60%
        ("Kansas City",  "MO", 39.30, -94.71, 0.02818),   # MO 40%
        ("Nashville",    "TN", 36.12, -86.68, 0.05805),
        ("Des Moines",   "IA", 41.53, -93.66, 0.04872),
        ("Louisville",   "KY", 38.17, -85.73, 0.03687),
    ],
    "slope_bcfd_per_hdd": 0.38,     # see above; RE-FIT in Phase 3
    "intercept_bcfd": 0.0,          # hinge model: the floor carries the baseload
    "floor_bcfd": 1.9,              # summer res/comm baseload, scaled from the 13-state 2.3
    "slope_sigma_frac": 0.10,       # 10% until re-fit, then 5%
    # GEFS 2 m temperature is under-dispersive, worse at longer leads (NOAA
    # GEFSv12 docs). Deviation-from-median multipliers for days 1-15; a guess
    # until Phase 3 fits spread-skill = 1 per lead day.
    "spread_inflation": [1.0] * 5 + [1.15] * 2 + [1.3] * 8,
    "basins": ["Appalachia / Pittsburgh"],   # names from hdd.BASINS
    "refit": False,                 # True once validate.py has replaced the placeholders
}

REGIONS = [MIDWEST]   # Phase 2 appends EAST, SOUTH_CENTRAL, MOUNTAIN, PACIFIC

for _r in REGIONS:
    assert abs(sum(c[4] for c in _r["cities"]) - 1) < 1e-4, f"{_r['key']} weights do not sum to 1"
    assert len(_r["spread_inflation"]) == 15, f"{_r['key']} needs 15 inflation factors"


def demand_bcfd(region, W):
    """Weighted HDD -> res/comm gas demand in Bcf/d. Works on any array shape."""
    return np.maximum(region["floor_bcfd"],
                      region["intercept_bcfd"] + region["slope_bcfd_per_hdd"] * W)


def calibrated_band(region, D):
    """Calibrated p10/p90 from per-member daily demand D[member, day] (Bcf/d).

    1. Stretch each member's deviation from the daily ensemble median by that
       day's spread_inflation factor (floor still applies).
    2. sigma_spread = std across members of the stretched values.
    3. sigma_slope = slope_sigma_frac x median, for the slope's own uncertainty.
    4. band = median -/+ 1.2816 x sqrt(sigma_spread^2 + sigma_slope^2).

    Returns daily lo/hi arrays (Bcf/d) and the 15-day total lo/hi (Bcf). The
    medians are the raw ensemble medians; only the width is calibrated."""
    median_d = np.median(D, axis=0)
    inflated = np.maximum(region["floor_bcfd"],
                          median_d + (D - median_d) * np.array(region["spread_inflation"]))
    frac = region["slope_sigma_frac"]

    sigma_d = np.hypot(inflated.std(axis=0), frac * median_d)
    median_s = np.median(D.sum(axis=1))
    sigma_s = np.hypot(inflated.sum(axis=1).std(), frac * median_s)
    return {
        "daily_lo": np.maximum(region["floor_bcfd"], median_d - Z90 * sigma_d),
        "daily_hi": median_d + Z90 * sigma_d,
        "total_lo": median_s - Z90 * sigma_s,
        "total_hi": median_s + Z90 * sigma_s,
    }
