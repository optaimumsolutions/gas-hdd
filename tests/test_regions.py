"""Tests for regions.py. Run with: python -m pytest"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from regions import MIDWEST, PHASE2_SLOPES, demand_bcfd  # noqa: E402


def test_member_sums_not_converted_percentiles():
    """With the cold kink, converting the p90 HDD path to Bcf is NOT the p90 of
    the members' Bcf. Member p90 by 15-day HDD total is a cold-snap member whose
    600 HDD include one day above the kink; a steady member with the same 600
    HDD would burn less gas. So Bcf percentiles must come from member sums."""
    W = np.full((31, 15), 20.0)                 # 27 mild members, 300 HDD each
    W[27] = 37.5                                # the cold-snap member: 14 x 37.5 ...
    W[27, 7] = 75.0                             # ... + one day at 75 HDD (kink is 45) = 600
    W[28:] = 50.0                               # 3 cold members, 750 HDD each

    p90_of_demand = np.percentile(demand_bcfd(MIDWEST, W).sum(axis=1), 90)
    p90_hdd_total = np.percentile(W.sum(axis=1), 90)                  # 600 = the cold-snap member
    demand_of_p90 = 15 * demand_bcfd(MIDWEST, p90_hdd_total / 15)     # as a steady 40 HDD/day

    assert p90_hdd_total == 600
    assert demand_of_p90 != p90_of_demand
    kink_bcf = MIDWEST["slope_bcfd_per_hdd"] * MIDWEST["cold_kink_extra_frac"] * (75 - 45)
    assert np.isclose(p90_of_demand - demand_of_p90, kink_bcf)


def test_phase2_slopes_sum_to_eia_national():
    assert np.isclose(sum(PHASE2_SLOPES.values()), 1.59)
    assert PHASE2_SLOPES["midwest"] == MIDWEST["slope_bcfd_per_hdd"]


def test_weights_sum_to_one():
    assert np.isclose(sum(c[4] for c in MIDWEST["cities"]), 1.0, atol=1e-4)


def test_floor_holds_in_summer():
    assert demand_bcfd(MIDWEST, np.array(0.0)) == MIDWEST["floor_bcfd"]
