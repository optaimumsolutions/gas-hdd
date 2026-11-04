# PRD: Calibrated HDD forecast (Neon archive + EMOS post-processing)

Status: decided 2026-11-04, not started. Owner: founder. Builder: Claude Code sessions.

## 0. Status ledger

| # | Milestone | State | Verified by |
|---|-----------|-------|-------------|
| M0 | `DATABASE_URL` secret in the gas-hdd repo | DONE 2026-11-04 15:23 UTC | `gh secret list -R optaimumsolutions/gas-hdd` shows it; `select 1` from this machine OK |
| M1 | `gefs` schema in Neon + `hdd.py` writes per-airport members to it | DONE 2026-11-04 | 2026-11-03 00Z run: 330 rows in `gefs.forecast_station`; Midwest member CSV rebuilt from Neon to 0.0005 |
| M2 | Backfill workflow: Nov 2024–Mar 2025 and Nov 2025–today, 00Z daily, Actions matrix | next | `SELECT count(DISTINCT run) FROM gefs.forecast_station` ≈ 300 |
| M3 | `observe.py` writes to Neon; truth backfilled for the same dates; 10-year climatology per airport | open | `gefs.observed` covers every backfilled forecast date |
| M4 | `score.py` reads from Neon; adds skill vs climatology; Bcf rows at the floor excluded | open | Accuracy tab shows per-lead n ≥ 200 |
| M5 | `calibrate.py`: EMOS per airport per lead, trained 2024-25, validated 2025-26 | open | holdout coverage 0.75–0.85 at every lead, MAE below raw |
| M6 | `hdd.py` applies the fit; chart shows corrected band + raw median grey line | open | page shows "calibrated" badge; `refit: true` in latest.json |
| M7 | Nightly refit on both winters; Accuracy tab reports out-of-sample from go-live | open | fit timestamp on page advances daily |
| M8 | ECMWF ENS second line (gated: after M6 and a few days of Pivotal checks) | open | second score column on Accuracy tab |

## 1. Goal

The dashboard's 15-day gas-weighted HDD forecast should have a median that is not
systematically off and a p10–p90 band that contains the observed value about 80% of
the time at every lead day. Both facts must be measured on the Accuracy tab, not
assumed. That is what "a spread inside a high-confidence interval, competitive with
the models traders already use" means in numbers.

### Success criteria (per index: `us_hdd`, `midwest_hdd`; per lead day 1–15)

| Metric | Target | Where |
|--------|--------|-------|
| Coverage of p10–p90 | 0.75 to 0.85 (target 0.80) | `docs/accuracy.json` `by_lead[].coverage` |
| Bias of the median | within ±0.5 HDD/day at days 1–7 | `by_lead[].bias` |
| Skill vs climatology | MAE below climatology MAE through at least day 10 | new `by_lead[].skill` |
| Holdout proof | the above hold on winter 2025-26 with a fit trained only on 2024-25 | `docs/accuracy.json` `holdout` |

## 2. Decisions (locked 2026-11-04)

| # | Decision | Chosen | Rejected |
|---|----------|--------|----------|
| D1 | Target band | p10–p90 at 80% coverage (as coded, `TARGET_COVERAGE`) | p5–p95; both bands |
| D2 | Archive of record | gas-hdd's own Neon project (`ep-shy-bread-bea9j2nu`, us-east-2, database `neondb`, Postgres 18, empty on 2026-11-04), schema `gefs` | ng-storage's Neon project; dual-write CSVs; read-only mirror |
| D3 | Backfill depth | two winters: Nov 2024–Mar 2025 and Nov 2025–now, one 00Z run per day (~300 runs) | one winter; all winters since 2020 |
| D4 | Archive grain | per airport, per member, per day: forecast Tmax and Tmin | index-only (today's CSVs); full grid |
| D5 | Fit method | EMOS / nonhomogeneous Gaussian regression per airport per lead, applied as a shift-and-scale of the 31 members | bias + stretch; index-level only |
| D6 | Fit window | train 2024-25, validate 2025-26, then nightly refit on everything archived | rolling 60 days; fit once and freeze |
| D7 | Baseline | 10-year daily normal per airport computed from IEM (same stations, same UTC-day rule) | NOAA 1991–2020 normals; persistence |
| D8 | ECMWF | after calibration ships; keep the Pivotal gate | add to backfill now; drop |
| D9 | Chart | corrected band primary, raw GEFS median as a thin grey line | corrected only; two bands |
| D10 | Backfill host | GitHub Actions matrix (20 jobs); nothing stored on the founder's machine | local machine; Render scheduler |
| D11 | Secret | `DATABASE_URL` repo secret (set 2026-11-04 from the connection the founder pasted) | separate Neon role (can be added later) |
| D12 | Schema ownership | gas-hdd owns schema `gefs` and its DDL (`sql/gefs.sql`, idempotent); the database holds nothing else | ng-storage numbered migration |

Decided by the builder, flagged here so they can be overruled:

- D13 Member storage: one row per (run, airport, date) with `tmax float4[31]`, `tmin float4[31]`.
  330 rows and about 100 KB per run; two winters ≈ 30 MB; a year of four runs a day ≈ 150 MB.
- D14 Bcf scoring: rows where the observed demand sits at the floor (1.9 Bcf/d) are excluded
  from coverage and MAE. A zero-width band at the floor scores 100% while saying nothing.
- D15 `history/*.csv` stop being written once M1 is verified. Existing files are loaded into
  Neon once (M1) and deleted from the repo in M4. `docs/latest.json` and the PNGs stay: the
  page is static and reads no database.
- D16 EMOS applied through the members (ensemble copula coupling): members at each airport are
  shifted and scaled so their mean and variance equal the fitted ones, ranks preserved. The
  HDD and Bcf conversions stay member by member, so cross-airport correlation and the cold
  kink keep working.

## 3. Why the data is needed, and where it lives

The errors being corrected are 1–3°F, different at Boston than at Houston, growing with
lead day. Seeing that through weather noise needs hundreds of forecast days per lead.
Two winters of daily runs is the minimum.

What moves where, per run:

| Stage | Size | Where | Kept? |
|-------|------|-------|-------|
| GRIB download (31 members × 65 leads, 2 m temperature) | ~300 MB | Actions runner `data/` cache | no, discarded with the runner |
| Lower-48 temperature cube in memory | ~1 GB RAM | runner | no |
| Per-airport Tmax/Tmin per member per day | ~100 KB | Neon `gefs.forecast_station` | yes, forever |
| `docs/latest.json` + PNGs | ~2 MB | git, served by Pages | yes, overwritten each run |

Nothing lands on the founder's machine. The 90 GB backfill figure is runner download
traffic from NOAA's free bucket, never stored.

## 4. Architecture

```
NOAA S3 (GEFS GRIB) ──hdd.py──► Neon gefs.forecast_station ◄──backfill.yml (matrix)
                         │
IEM ASOS hourly ──observe.py──► Neon gefs.observed, gefs.climatology
                         │
                         ▼
        score.py ──► gefs.scores + docs/accuracy.json + PNGs
        calibrate.py ──► gefs.emos_fit (per airport, per lead)
                         │
        hdd.py (next run) reads gefs.emos_fit ──► corrected members ──► docs/latest.json, PNGs
```

Scripts stay flat Python files, no framework. One new module `db.py` (connection from
`DATABASE_URL`, `psycopg[binary]`, schema-qualified names, no `SET search_path` through
the pooler). One new script `calibrate.py`. One new workflow `backfill.yml`. `run.yml`
and `observe.yml` gain the secret and lose the `history/*` commit pattern.

### Neon schema `gefs` (DDL in `sql/gefs.sql`, `CREATE ... IF NOT EXISTS`)

```sql
forecast_station(run timestamptz, station text, date date,
                 tmax real[31], tmin real[31], model text default 'gefs',
                 ingested_at timestamptz, primary key (model, run, station, date))
observed(station text, date date, tmax real, tmin real, hours int, source text,
         fetched_at timestamptz, primary key (station, date))      -- never updated
climatology(station text, doy int, tmean_normal real, tmean_sd real, years int,
            primary key (station, doy))
emos_fit(model text, station text, lead_day int, a real, b real, c real, d real,
         n int, train_from date, train_to date, fitted_at timestamptz,
         primary key (model, station, lead_day))
scores(model text, run timestamptz, date date, lead_day int, index text,
       median real, p10 real, p90 real, observed real, raw_median real,
       primary key (model, run, date, index))                       -- rebuilt each run
```

Indexes derived from `forecast_station` in Python (weights from `hdd.CITIES` and
`regions.REGIONS`), so a weight change rescored against the whole archive is one rerun of
`score.py`. ng-storage lives in a different Neon project, so a Midwest Bcf vs weekly
storage comparison reads both databases from Python; there is no cross-database JOIN.

## 5. Milestones in build order

**M0 — Secret. DONE 2026-11-04.** `DATABASE_URL` is set on the repo (pooler host,
`sslmode=require&channel_binding=require`). Local runs read the same value from a
git-ignored `.env` in the repo root.

**M1 — Neon write path (half a day).** `db.py`, `sql/gefs.sql`, `hdd.py` gains
`station_extremes()` returning per-airport Tmax/Tmin per member per day and writes them
to `forecast_station` in step 7 (upsert, idempotent). `--no-docs` flag skips charts and
`latest.json` for backfill jobs. One-time loader for the six existing `history/*.csv`
runs is NOT possible at airport grain (the CSVs hold only indexes): those six runs are
re-ingested from S3 by the backfill instead. Test: run `hdd.py --run 2026-11-04T00`
locally against Neon, then `SELECT count(*)` = 22 × 15 = 330.

**M2 — Backfill (half a day to build, ~1.5 h to run).** `backfill.yml` with
`workflow_dispatch` inputs `from`, `to`; a first job lists the 00Z runs in range and
splits them into 20 shards; the matrix runs `hdd.py --run <r> --no-docs` per run with
`continue-on-error` so one throttled run never kills a shard. A final job prints which
runs are missing from Neon so the founder can re-dispatch just those. Observed truth
must not be fetched here (IEM rate limits): that is M3.

**M3 — Truth and climatology (half a day).** `observe.py --from --to` writes
`gefs.observed` in chunks of 30 days with a pause (IEM asks for that). Climatology:
one pull per station of 10 winters of hourly data → daily (Tmax+Tmin)/2 → mean and sd
per day-of-year, smoothed with a 15-day window, into `gefs.climatology`. Test: every
date in `forecast_station` has a row in `observed` for all 22 stations.

**M4 — Scoring from Neon (half a day).** `score.py` reads `forecast_station` +
`observed`, computes indexes with the current weights, writes `gefs.scores` and
`docs/accuracy.json`. Adds `climatology_mae` and `skill = 1 − MAE/MAE_clim` per lead.
Applies D14. Deletes `history/*.csv` from the repo. Test: `tests/test_accuracy.py`
rewritten against a temporary schema; Accuracy tab shows n ≥ 200 per lead.

**M5 — EMOS fit (one day).** `calibrate.py`: for each (station, lead_day), fit
`obs ~ N(a + b·ens_mean, c + d·ens_var)` on daily-mean temperature by minimising CRPS
(scipy, closed-form Gaussian CRPS). Train on 2024-11-01…2025-03-31, score on
2025-11-01…today, write `gefs.emos_fit` plus a `holdout` block in `accuracy.json`
(raw vs corrected MAE and coverage per lead). Gate to M6: holdout coverage inside
0.75–0.85 at every lead and corrected MAE ≤ raw MAE at every lead. If a station fails the
gate its fit is stored with `b=1, d=1, a=0, c=0` (identity) and listed on the page.

**M6 — Apply in the forecast (half a day).** `hdd.py` loads `emos_fit`, transforms the
members per station per lead (D16), computes everything downstream from the corrected
members, and also keeps the raw index median for the grey line (D9). `latest.json` gains
`calibration: {fitted_at, train_from, train_to, stations_identity: [...]}`;
`regions.py` `spread_inflation` and `slope_sigma_frac` placeholders are retired and
`refit` becomes true. README changelog entry, as was done for the Tmax/Tmin change.

**M7 — Nightly refit (1 hour).** `observe.yml` runs `calibrate.py` after `score.py`,
training on every archived winter day up to yesterday. The Accuracy tab's headline
statistics from go-live onward are out of sample by construction (each day's forecast
was corrected with a fit that had not seen that day).

**M8 — ECMWF ENS (one day, gated).** herbie `model="ecmwf", product="enfo"`, 0.25°,
51 members, 15 days, `model='ecmwf'` rows in the same tables, own EMOS fit, second line
on the chart, second column on the Accuracy tab. Backfill limited by ECMWF open-data
retention (checked at build time). Starts only after M6 is live and the founder has
compared the GEFS line with Pivotal for a few days.

## 6. Risks and mitigations

- **Neon unreachable during a run.** Forecast step still writes `latest.json` and PNGs
  and the job commits them; the Neon write is retried 5× then logged as a warning, and
  the next run re-ingests the missed run (upsert). Calibration falls back to the last
  fit in `latest.json` is NOT done: without Neon the run is published uncorrected and
  says so in the yellow note.
- **IEM rate limits on a 10-winter pull.** Climatology pull is one station per request
  with a 2 s pause, run once, cached in Neon; never part of the nightly job.
- **Pooler session state.** Only schema-qualified names, no `SET search_path`
  (ng-storage gotcha, 2026-10-26).
- **EMOS overfits at long leads.** 300 samples per (station, lead) is adequate for four
  parameters; the holdout gate (M5) is the guard. Day 13–15 fits typically collapse
  toward climatology, which is the correct answer there.
- **Backfill cost.** Public repo: Actions minutes are free. NOAA bucket: free egress.
  Neon storage: ~30 MB for the backfill.

## 7. Out of scope

Hourly-resolution demand, power burn, other regions (Phase 2 of regions.py, unchanged),
any paid data, any server beyond the existing Actions + Pages + Neon.

## 8. Open items for the founder

1. Confirm D13–D16 stand, or say which to change, before M1 starts.
