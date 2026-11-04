-- Schema for the gas-hdd forecast archive (PRD-calibration.md, D2/D4/D12/D13).
-- Idempotent: db.ensure_schema() runs this file on every connection.
-- Names are always schema-qualified in code; nothing relies on search_path.

CREATE SCHEMA IF NOT EXISTS gefs;

-- One row per run x airport x forecast day. tmax/tmin hold the 31 members'
-- daily max / min 2 m temperature in F (index 0 = control, 1..30 = perturbed).
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
    PRIMARY KEY (model, run, station, date),
    CHECK (cardinality(tmax) = 31 AND cardinality(tmin) = 31)
);

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
