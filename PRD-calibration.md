# PRD: Calibrated HDD forecast (Neon archive + EMOS post-processing)

Status: decided 2026-11-04, not started. Owner: founder. Builder: Claude Code sessions.

## 0. Status ledger

| # | Milestone | State | Verified by |
|---|-----------|-------|-------------|
| M0 | `DATABASE_URL` secret in the gas-hdd repo | DONE 2026-11-04 15:23 UTC | `gh secret list -R optaimumsolutions/gas-hdd` shows it; `select 1` from this machine OK |
| M1 | `gefs` schema in Neon + `hdd.py` writes per-airport members to it | DONE 2026-11-04 | 2026-11-03 00Z run: 330 rows in `gefs.forecast_station`; Midwest member CSV rebuilt from Neon to 0.0005 |
| M2 | Backfill workflow: Nov 2024–Mar 2025 and Nov 2025–today, 00Z daily, Actions matrix | DONE 2026-11-04 21:00 UTC: 337 runs, 0 missing, 0 failed shards | `SELECT count(DISTINCT run) FROM gefs.forecast_station` ≈ 300 |
| M3 | `observe.py` writes to Neon; truth backfilled for the same dates; 10-year climatology per airport | DONE 2026-11-04 (3177091) | `gefs.observed` covers every backfilled forecast date |
| M4 | `score.py` reads from Neon; adds skill vs climatology; Bcf rows at the floor excluded | DONE 2026-11-04 | first Neon scoring: 236 runs, 10,060 forecast-days, n ≥ 221 per lead |
| M5 | `calibrate.py`: EMOS per airport per lead + index spread stage, trained 2024-25, validated 2025-26 | BUILT 2026-11-04; holdout: 20 of 30 lead checks inside 75–85%, 7 over by ≤0.02, US days 1–3 at 71–75% (see M5 notes); founder decision on shipping | holdout coverage 0.75–0.85 at every lead, MAE below raw |
| M6 | `hdd.py` applies the fit; chart shows corrected band + raw median grey line | DONE 2026-11-04 (founder: "ship") | local run 2026-11-03 00Z: 330 fits + 30 multipliers applied, raw grey line, `calibration` block in latest.json, `refit: true` |
| M7 | Nightly refit on both winters; Accuracy tab reports out-of-sample from go-live | DONE 2026-11-04 (`calibrate.py --nightly` in observe.yml 09:00 UTC; identity verdicts inherited from the holdout run) | fit timestamp on page advances daily |
| M8 | ECMWF ENS second line (gate lifted 2026-11-05, grill G3) | IN PROGRESS: ecmwf_ens.py + ecmwf.yml + backfill --model + per-model score/calibrate + second line on charts shipped; first GitHub run, two-winter backfill, ECMWF fit and scorecard pending |
| M9 | Band scorecard: grade the published calibrated band (coverage, error/spread, MAE) with raw alongside; `gefs.forecast_index` written per run; "as published" block accrues from go-live | DONE 2026-11-05 |
| M10 | Totals scorecard (15-day total + EIA-week sums: coverage, error/spread, MAE), a fitted spread multiplier for the 15-day total (key (index, 0) in spread_fit), and no-skill shading on the charts | DONE 2026-11-05: total holdout US 0.54→0.73, Midwest 0.60→0.77; weeks graded unstretched (0.77/0.79 week 1) |
| M11 | Demand fitted to data: `fit_demand.py` (EIA monthly res+comm per state + weekly Midwest storage vs 10 winters of observed airport HDD); regions.py refit (Midwest 0.96 + 0.354·HDD, floor 1.82; storage change = 72.7 − 0.471·weekly HDD); reference tables A/B on the US and Midwest tabs | DONE 2026-11-05 | second score column on Accuracy tab |

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

**M2 — Backfill (built 2026-11-04).** `backfill.yml` with `workflow_dispatch` inputs
`from`, `to`, `shards`; the `plan` job asks Neon which 00Z runs in the range are missing
and splits them into shards; the matrix runs `backfill.py run <runs>` (each run =
`hdd.py --run <r> --no-docs`, GRIB cache deleted after each, THREADS=4, failures
collected, `fail-fast: false`); the `report` job lists what is still missing, so
re-dispatching the same range retries only the leftovers. Dispatch windows (summer
skipped: no HDD signal): 2024-11-01..2025-03-31, 2025-11-01..2026-03-31,
2026-10-01..today. Observed truth is not fetched here (IEM rate limits): that is M3.

**M3 — Truth and climatology (built 2026-11-04).** `observe.py` now writes every
complete station-day to `gefs.observed` (insert-only, judged per station: one short
station no longer discards the other 21 for that day) and still appends `observed.csv`
until M4. `--from/--to` backfill a range in 31-day chunks with a 2 s pause; both
workflows pass the secret. `--climatology`: one IEM pull per airport, 10 years of hourly
data → daily (Tmax+Tmin)/2 → mean and sd per day of year (Feb 29 folded into day 60),
15-day circular smoothing, upserted into `gefs.climatology`. Truth windows loaded:
2024-11-01..2025-04-14, 2025-11-01..2026-04-14, 2026-10-01..2026-11-03. Test: every
date in `forecast_station` has a row in `observed` for all 22 stations.

**M4 — Scoring from Neon (built 2026-11-04).** `hdd.index_members(data, cities)` is now
the one definition of every index: the live forecast, `score.py` and `validate_log.py`
all call it on airport-grain data (live arrays or `db.load_forecasts`). Side effect: the
US index now samples the grid at the airports (O'Hare, not downtown Chicago) like the
Midwest index always did; a small step change, logged in the README. `score.py` reads
`forecast_station` + `observed` + `climatology`, writes `gefs.scores` and
`docs/accuracy.json` with `clim_mae` and `skill = 1 − MAE/clim_mae` per lead and month;
D14 applied (`excluded_at_floor` reported). Run-to-run deltas in `hdd.py` come from Neon.
`history/` keeps only `observed.csv` (home of the CPC reference) and `validation.csv`.
First result, two winters: US HDD skill 0.77 at day 1 falling to 0.04 at day 15, so the
forecast carries no information beyond climatology after day 13; raw p10–p90 coverage
0.26 at day 1 rising to ~0.70, never reaching 0.80.

**M5 — EMOS fit (one day).** `calibrate.py`: for each (station, lead_day), fit
`obs ~ N(a + b·ens_mean, c + d·ens_var)` on daily-mean temperature by minimising CRPS
(scipy, closed-form Gaussian CRPS). Train on 2024-11-01…2025-03-31, score on
2025-11-01…today, write `gefs.emos_fit` plus a `holdout` block in `accuracy.json`
(raw vs corrected MAE and coverage per lead). Gate to M6: holdout coverage inside
0.75–0.85 at every lead and corrected MAE ≤ raw MAE at every lead. If a station fails the
gate its fit is stored with `b=1, d=1, a=0, c=0` (identity) and listed on the page.

M5 first result (2026-11-04 19:40 UTC, training winter only 74 of 151 runs loaded, mostly
Nov–Dec 2024): station-level holdout coverage rose from 0.30 to 0.70 at day 1 and MAE fell
10%; index-level US coverage 0.27 → 0.64 at day 1, 0.65 → 0.83 at day 6, unchanged ≥ day 9
because the per-station CRPS gate kept 149 of 330 pairs at identity (long leads). Two
likely causes, both testable once the backfill lands: the fit saw no January–February
errors (larger than November's), and index coverage is capped by the raw ensemble's
cross-airport correlation, which ECC preserves. Full-winter refit (21:00 UTC): airport fits cover 80% in-sample and 75–81% on the holdout at
days 1–4, but beyond day 5 the CRPS-fitted Gaussian covers only 71–79% even in-sample
(heavy tails), and the index band is narrower still (70–75% in-sample at days 2–8) because
airport errors are more correlated than the members are. **Second stage added:** one spread
multiplier k per HDD index per lead, fitted by bisection to 80% training coverage after the
airport stage, each lead pooled with its two neighbours (k = 1.07–1.41), stored in
`gefs.spread_fit`, inherited by the Bcf conversion (stretch first, convert after). MAE test
tolerance 2%. Holdout result (n = 150 per lead):

| Lead | US raw cov | US cal cov | US raw MAE | US cal MAE | Midwest raw cov | Midwest cal cov |
|---|---|---|---|---|---|---|
| 1 | 0.27 | 0.71 | 1.03 | 1.04 | 0.37 | 0.78 |
| 3 | 0.45 | 0.73 | 1.32 | 1.31 | 0.48 | 0.73 |
| 5 | 0.55 | 0.81 | 1.66 | 1.68 | 0.57 | 0.86 |
| 8 | 0.71 | 0.79 | 2.63 | 2.62 | 0.70 | 0.82 |
| 10 | 0.71 | 0.87 | 3.14 | 3.14 | 0.69 | 0.82 |
| 15 | 0.68 | 0.75 | 4.56 | 4.50 | 0.65 | 0.79 |

20 of 30 (index, lead) checks inside 0.75–0.85; 7 over by at most 0.02 (within the ±0.033
sampling error at n = 150); 3 under: US days 1 and 3 and Midwest day 3 at 0.71–0.73. Day 1
is the one lead where the two winters differ (holdout airport errors 19% wider in sigma
units than the training fit), which the nightly refit on both winters (M7) will absorb.

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

**M8 — ECMWF ENS (built 2026-11-05).** `ecmwf_ens.py`: herbie `model="ecmwf", product="enfo",
resolution="0p25"`, only the `2t` messages (51 members) of the 6-hourly steps 0–360 h,
read at the 22 airports step by step and deleted (a run never sits on disk); daily
Tmax/Tmin = max/min of the 6-hourly values valid 00–24Z (the open-data ensemble has no
6-h max/min fields; the per-airport EMOS fit absorbs the definition difference). 00Z runs
give 15 days, 12Z runs 14. Rows carry `model='ecmwf'`; `score.py --model ecmwf` and
`calibrate.py --model ecmwf` write `accuracy_ecmwf.json` / `calibration_ecmwf.json` and
their own `emos_fit` / `spread_fit` rows. `hdd.py` draws the newest archived ECMWF run,
corrected with its own fit when one exists, as a dashed red median on every chart and
quotes its shared-day total next to GEFS's. Triggers: `ecmwf.yml` 08:45 / 20:45 UTC plus
the Render dispatch job. Backfill: `backfill.yml` with `model=ecmwf` over the same windows
(bucket holds `enfo` 0.25° from 2024-01). Gotcha: the `ecmwf-forecasts` bucket rate-limits
a single IP hard (503 SlowDown); THREADS stays at 3 and the real test runs on GitHub.

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

## 9. Grill of 2026-11-05 (accuracy and what to improve)

| # | Decision | Chosen |
|---|----------|--------|
| G1 | Judged on | the 15-day total and the EIA-week sums, scored like the daily values (M10) |
| G2 | No-skill days | keep the 15-day total; shade days with skill < 20% on the chart (M10) |
| G3 | ECMWF | gate lifted; same backfill (bucket has enfo 0.25° from Jan 2024), own calibration and scorecard, second line; blend only once scored (M8) |
| G4 | Demand fit basis | EIA monthly res+comm per state and weekly Midwest storage vs this index's observed HDD, 10 winters stored (observe.py --climatology now writes the daily rows), last 5 winters set the numbers (M11) |
| G5 | Fit storage | `fit_demand.py` prints; `regions.py` updated by hand with date and source; refit is a deliberate act |
| G6 | Tables | A: HDD level → Bcf/d with past range; B: deviation from normal → ΔBcf/d, ΔBcf/week, Δstorage/week; on the US and every region tab |
| G7 | Storage card | fitted weekly relation replaces the 1.12–1.26 ratio; EIA's sign convention (negative = withdrawal) |
