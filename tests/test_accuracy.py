"""Tests for observe.py and score.py. Run with: python -m pytest"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hdd  # noqa: E402
import observe  # noqa: E402
import score  # noqa: E402


def test_scoring_synthetic_run_gives_expected_mae_and_coverage(tmp_path):
    """31 members evenly spread 10..20 every day: median 15, p10 11, p90 19.
    Observed [15, 10, 20, 16, 14] -> errors [0, 5, -5, -1, 1]: MAE 2.4,
    bias 0, RMSE sqrt(52/5), and 3 of 5 days inside the band."""
    dates = [f"2026-11-{d:02d}" for d in range(1, 6)]
    members = np.tile(np.linspace(10, 20, 31), (5, 1))          # [day, member]
    frame = pd.DataFrame(members, index=pd.Index(dates, name="date"),
                         columns=[f"m{m:02d}" for m in range(31)])
    frame.to_csv(tmp_path / "hdd_20261101_00z.csv")
    pd.DataFrame({"us_hdd": [15.0, 10.0, 20.0, 16.0, 14.0]},
                 index=pd.Index(dates, name="date")).to_csv(tmp_path / "observed.csv")

    scores = score.build_scores(tmp_path)
    us = scores[scores["index"] == "us_hdd"].sort_values("lead_day")
    assert list(us["lead_day"]) == [1, 2, 3, 4, 5]
    assert np.allclose(us["median"], 15) and np.allclose(us["p10"], 11) and np.allclose(us["p90"], 19)
    assert np.allclose(us["error"], [0, 5, -5, -1, 1])
    assert list(us["inside_band"]) == [True, False, False, True, True]

    summary = score.aggregate(scores, score.load_observed(tmp_path), history=tmp_path)
    overall = summary["overall"]["us_hdd"]
    assert overall["n"] == 5
    assert np.isclose(overall["mae"], 2.4) and overall["bias"] == 0
    assert np.isclose(overall["rmse"], round(np.sqrt(52 / 5), 2))
    assert np.isclose(overall["coverage"], 0.6)
    # one score per lead day is below MIN_N, so every per-lead statistic is "insufficient data"
    assert all(r["n"] == 1 and r["mae"] is None for r in summary["by_lead"]["us_hdd"][:5])
    assert summary["storage"]["n"] == 0 and summary["cpc"]["n"] == 0


def test_observed_hdd_is_the_forecast_code_path(monkeypatch):
    """Both the forecast (hdd.city_hdd) and the observation (observe.observed_row)
    must go through hdd.weighted_hdd; a spy on it has to see both calls."""
    calls = []
    real = hdd.weighted_hdd

    def spy(daily_mean_f, weights):
        calls.append(np.shape(daily_mean_f))
        return real(daily_mean_f, weights)

    monkeypatch.setattr(hdd, "weighted_hdd", spy)

    # forecast path: a 2x2 grid at a constant 50 F (283.15 K), one UTC day
    lat, lon = np.array([41.5, 42.0]), np.array([-88.0, -87.5])
    TX = np.full((len(hdd.MEMBERS), 5, 2, 2), 283.15, dtype=np.float32)
    forecast = hdd.city_hdd(TX, TX, lat, lon, [[0, 1, 2, 3]], [(41.9, -87.6, 1.0)])
    assert np.allclose(forecast, 15.0)                          # 65 - 50

    # observed path: every station at 50 F
    row = observe.observed_row({s: 50.0 for s in observe.STATIONS})
    assert row["us_hdd"] == 15.0 and row["midwest_hdd"] == 15.0
    assert len(calls) >= 2                                      # one from each path


def test_refetched_observation_is_logged_not_overwritten(capsys):
    cols = ["us_hdd", "midwest_hdd", "cpc_midwest_hdd", "source", "fetched_utc"]
    existing = pd.DataFrame([[10.0, 12.0, np.nan, "x", "t0"]],
                            index=pd.Index(["2026-11-01"], name="date"), columns=cols)
    new = pd.DataFrame([[11.0, 12.2, 13.0, "x", "t1"], [9.0, 9.5, np.nan, "x", "t1"]],
                       index=pd.Index(["2026-11-01", "2026-11-02"], name="date"), columns=cols)

    merged, added = observe.merge(existing, new)

    assert added == ["2026-11-02"]
    assert merged.loc["2026-11-01", "us_hdd"] == 10.0            # first value kept
    assert merged.loc["2026-11-01", "midwest_hdd"] == 12.0
    assert merged.loc["2026-11-01", "cpc_midwest_hdd"] == 13.0   # only a blank reference cell is filled
    assert merged.loc["2026-11-02", "us_hdd"] == 9.0
    out = capsys.readouterr().out
    assert "revision ignored: 2026-11-01 us_hdd stored 10.00, refetched 11.00" in out
    assert "midwest_hdd stored" not in out                        # 0.2 HDD is inside the tolerance
