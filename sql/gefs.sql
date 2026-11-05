-- Schema for the gas-hdd forecast archive (PRD-calibration.md, D2/D4/D12/D13).
-- Idempotent: db.ensure_schema() runs this file on every connection.
-- Names are always schema-qualified in code; nothing relies on search_path.

CREATE SCHEMA IF NOT EXISTS gefs;

-- One row per run x airport x forecast day. tmax/tmin hold every member's
-- daily max / min 2 m temperature in F (GEFS: 31, index 0 = control; ECMWF: 51).
-- Daily mean = (tmax + tmin) / 2, the same definition observe.py uses.
CREATE TABLE IF NOT EXISTS gefs.forecast_station (
    model       text        NOT NULL DEFAULT 'gefs',
    run         timestamptz NOT NULL,
    station     text        NOT NULL,          -- ICAO airport code, e.g. ORD
    date        date        NOT NULL,          -- UTC calendar day the values cover
    lead_day    smallint    NOT NULL,          -- 1 = first full forecast day of the run
    tmax        real[]      NOT NULL,          -- 31 members, F
    tmin        real[]      NOT NULL,          -- 31 members, F
    ingested_at timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (model, run, station, date)
);
-- GEFS has 31 members, ECMWF 51: the old fixed-size check is dropped.
ALTER TABLE gefs.forecast_station DROP CONSTRAINT IF EXISTS forecast_station_check;

-- Where each airport sits on the model grid, so a coordinate change is visible.
CREATE TABLE IF NOT EXISTS gefs.station (
    station text PRIMARY KEY,
    lat     real NOT NULL,
    lon     real NOT NULL
);

-- Observed truth per airport and UTC day from IEM hourly ASOS: max and min of
-- the hourly values, how many hourly reports the day had. Rows are never
-- updated (observe.py inserts with ON CONFLICT DO NOTHING).
CREATE TABLE IF NOT EXISTS gefs.observed (
    station    text        NOT NULL,
    date       date        NOT NULL,
    tmax       real        NOT NULL,          -- F
    tmin       real        NOT NULL,          -- F
    hours      smallint    NOT NULL,
    source     text        NOT NULL,
    fetched_at timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (station, date)
);

-- Climatology per airport and day of year (1..366, Feb 29 folded into 60):
-- mean and sd of the daily mean temperature (F) over `years` years, smoothed
-- with a 15-day circular window. The skill baseline for score.py.
CREATE TABLE IF NOT EXISTS gefs.climatology (
    station      text     NOT NULL,
    doy          smallint NOT NULL,
    tmean_normal real     NOT NULL,
    tmean_sd     real     NOT NULL,
    years        smallint NOT NULL,
    built_at     timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (station, doy)
);

-- One row per run x forecast date x index, rebuilt from scratch by score.py
-- (observed, clim = the climatology baseline's value for that date).
CREATE TABLE IF NOT EXISTS gefs.scores (
    model     text    NOT NULL,
    run       timestamptz NOT NULL,
    date      date    NOT NULL,
    lead_day  smallint NOT NULL,
    index     text    NOT NULL,
    median    real    NOT NULL,
    p10       real    NOT NULL,
    p90       real    NOT NULL,
    observed  real    NOT NULL,
    clim      real,
    PRIMARY KEY (model, run, date, index)
);

-- EMOS / nonhomogeneous Gaussian regression per airport and lead day
-- (PRD D5): observed daily mean ~ N(a + b*m, c^2 + d^2*s^2) with m, s the
-- ensemble mean and standard deviation of the daily mean. Identity = (0,1,0,1).
CREATE TABLE IF NOT EXISTS gefs.emos_fit (
    model      text     NOT NULL,
    station    text     NOT NULL,
    lead_day   smallint NOT NULL,
    a          real     NOT NULL,
    b          real     NOT NULL,
    c          real     NOT NULL,
    d          real     NOT NULL,
    n          int      NOT NULL,           -- training samples
    identity   boolean  NOT NULL DEFAULT false,   -- fit failed the holdout gate; stored as identity
    train_from date     NOT NULL,
    train_to   date     NOT NULL,
    fitted_at  timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (model, station, lead_day)
);

-- Second-stage spread multiplier per index and lead day (calibrate.py):
-- index members are stretched about their median by k before percentiles.
CREATE TABLE IF NOT EXISTS gefs.spread_fit (
    model      text     NOT NULL,
    index      text     NOT NULL,
    lead_day   smallint NOT NULL,
    k          real     NOT NULL,
    n          int      NOT NULL,
    train_from date     NOT NULL,
    train_to   date     NOT NULL,
    fitted_at  timestamptz NOT NULL DEFAULT now(),
    PRIMARY KEY (model, index, lead_day)
);

-- What each run PUBLISHED per index and day (calibrated and raw), written by
-- hdd.py: the numbers score.py grades, so the scorecard is for the band that
-- was on the page, out of sample from go-live (2026-11-04 22:00 UTC).
CREATE TABLE IF NOT EXISTS gefs.forecast_index (
    model      text        NOT NULL DEFAULT 'gefs',
    run        timestamptz NOT NULL,
    index      text        NOT NULL,          -- us_hdd, midwest_hdd, midwest_bcf
    date       date        NOT NULL,
    lead_day   smallint    NOT NULL,
    median     real NOT NULL, p10 real NOT NULL, p90 real NOT NULL, sd real NOT NULL,
    raw_median real NOT NULL, raw_p10 real NOT NULL, raw_p90 real NOT NULL, raw_sd real NOT NULL,
    PRIMARY KEY (model, run, index, date)
);

ALTER TABLE gefs.scores ADD COLUMN IF NOT EXISTS sd real;
ALTER TABLE gefs.scores ADD COLUMN IF NOT EXISTS raw_median real;
ALTER TABLE gefs.scores ADD COLUMN IF NOT EXISTS raw_p10 real;
ALTER TABLE gefs.scores ADD COLUMN IF NOT EXISTS raw_p90 real;
ALTER TABLE gefs.scores ADD COLUMN IF NOT EXISTS published boolean NOT NULL DEFAULT false;
