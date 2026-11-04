"""
db.py - the one place that talks to the Neon archive (PRD-calibration.md).

DATABASE_URL comes from the environment, or from a git-ignored .env file in the
repo root for local runs. The schema in sql/gefs.sql is applied on every
connect, so a fresh database needs no setup step. Every table name is written
schema-qualified; search_path is never set (Neon's pooler would leak it).
"""

import os
from datetime import timezone
from pathlib import Path

import numpy as np
import psycopg

ROOT = Path(__file__).resolve().parent
SQL = ROOT / "sql" / "gefs.sql"


def database_url():
    """DATABASE_URL from the environment, else from .env; None if neither has it."""
    url = os.environ.get("DATABASE_URL")
    if url:
        return url
    env = ROOT / ".env"
    if env.exists():
        for line in env.read_text(encoding="utf-8").splitlines():
            if line.startswith("DATABASE_URL="):
                return line.split("=", 1)[1].strip().strip("'\"")
    return None


def connect():
    """Open a connection and make sure the gefs schema exists. Raises if
    DATABASE_URL is not set anywhere."""
    url = database_url()
    if not url:
        raise RuntimeError("DATABASE_URL is not set (environment or .env)")
    conn = psycopg.connect(url, connect_timeout=30)
    with conn.cursor() as cur:
        cur.execute(SQL.read_text(encoding="utf-8"))
    conn.commit()
    return conn


def write_forecast(conn, run, days, stations, tmax, tmin, model="gefs"):
    """Upsert one run. stations = {icao: (lat, lon)} in the order of the third
    axis of tmax/tmin [member, day, station] (F). Returns the row count."""
    codes = list(stations)
    rows = []
    for di, day in enumerate(days):
        for si, code in enumerate(codes):
            rows.append((model, run, code, day, di + 1,
                         [float(x) for x in tmax[:, di, si]], [float(x) for x in tmin[:, di, si]]))
    with conn.cursor() as cur:
        cur.executemany(
            "INSERT INTO gefs.station (station, lat, lon) VALUES (%s, %s, %s) "
            "ON CONFLICT (station) DO UPDATE SET lat = EXCLUDED.lat, lon = EXCLUDED.lon",
            [(c, float(stations[c][0]), float(stations[c][1])) for c in codes])
        cur.executemany(
            "INSERT INTO gefs.forecast_station (model, run, station, date, lead_day, tmax, tmin) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s) "
            "ON CONFLICT (model, run, station, date) DO UPDATE SET "
            "lead_day = EXCLUDED.lead_day, tmax = EXCLUDED.tmax, tmin = EXCLUDED.tmin, ingested_at = now()",
            rows)
    conn.commit()
    return len(rows)


def count_forecast(conn, run, model="gefs"):
    with conn.cursor() as cur:
        cur.execute("SELECT count(*) FROM gefs.forecast_station WHERE model = %s AND run = %s", (model, run))
        return cur.fetchone()[0]


def write_observed(conn, rows):
    """Insert observed days; rows = (station, date, tmax, tmin, hours, source).
    Existing rows are left untouched. Returns how many were new."""
    with conn.cursor() as cur:
        cur.executemany(
            "INSERT INTO gefs.observed (station, date, tmax, tmin, hours, source) "
            "VALUES (%s, %s, %s, %s, %s, %s) ON CONFLICT (station, date) DO NOTHING",
            rows, returning=False)
        new = cur.rowcount
    conn.commit()
    return new


def write_climatology(conn, rows):
    """Upsert climatology rows (station, doy, tmean_normal, tmean_sd, years)."""
    with conn.cursor() as cur:
        cur.executemany(
            "INSERT INTO gefs.climatology (station, doy, tmean_normal, tmean_sd, years) "
            "VALUES (%s, %s, %s, %s, %s) ON CONFLICT (station, doy) DO UPDATE SET "
            "tmean_normal = EXCLUDED.tmean_normal, tmean_sd = EXCLUDED.tmean_sd, "
            "years = EXCLUDED.years, built_at = now()",
            rows)
    conn.commit()
    return len(rows)


# ------------------------------------------------------------------ readers
def _naive_utc(ts):
    return ts.astimezone(timezone.utc).replace(tzinfo=None)


def list_runs(conn, model="gefs"):
    """Archived runs as naive UTC datetimes, oldest first."""
    with conn.cursor() as cur:
        cur.execute("SELECT DISTINCT run FROM gefs.forecast_station WHERE model = %s ORDER BY run", (model,))
        return [_naive_utc(r[0]) for r in cur.fetchall()]


def _assemble(rows):
    """rows (run, station, date, tmax, tmin) sorted by run, date -> {run:
    {"dates": [date...], "stations": {code: (tmax[day, member], tmin[day, member])}}}."""
    out = {}
    for run, station, day, tmax, tmin in rows:
        r = out.setdefault(_naive_utc(run), {"dates": [], "stations": {}})
        if not r["dates"] or r["dates"][-1] != day:
            r["dates"].append(day)
        s = r["stations"].setdefault(station, ([], []))
        s[0].append(tmax)
        s[1].append(tmin)
    for r in out.values():
        r["stations"] = {code: (np.array(tx, dtype=float), np.array(tn, dtype=float))
                         for code, (tx, tn) in r["stations"].items()}
    return out


def load_forecasts(conn, model="gefs", run=None):
    """Every archived run (or one run) at airport grain; see _assemble."""
    sql = "SELECT run, station, date, tmax, tmin FROM gefs.forecast_station WHERE model = %s"
    params = [model]
    if run is not None:
        sql += " AND run = %s"
        params.append(run.replace(tzinfo=timezone.utc))
    with conn.cursor() as cur:
        cur.execute(sql + " ORDER BY run, date, station", params)
        return _assemble(cur.fetchall())


def load_observed_means(conn):
    """{date: {station: daily mean F}} for every observed station-day."""
    out = {}
    with conn.cursor() as cur:
        cur.execute("SELECT station, date, tmax, tmin FROM gefs.observed")
        for station, day, tmax, tmin in cur.fetchall():
            out.setdefault(day, {})[station] = (tmax + tmin) / 2
    return out


def load_climatology(conn):
    """{station: array[365] of normal daily mean F, index 0 = day of year 1}."""
    out = {}
    with conn.cursor() as cur:
        cur.execute("SELECT station, doy, tmean_normal FROM gefs.climatology ORDER BY station, doy")
        for station, doy, normal in cur.fetchall():
            out.setdefault(station, np.full(365, np.nan))[doy - 1] = normal
    return out


def write_scores(conn, rows, model="gefs"):
    """Replace every score row for the model; rows = (run, date, lead_day,
    index, median, p10, p90, observed, clim)."""
    with conn.cursor() as cur:
        cur.execute("DELETE FROM gefs.scores WHERE model = %s", (model,))
        cur.executemany(
            "INSERT INTO gefs.scores (model, run, date, lead_day, index, median, p10, p90, observed, clim) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
            [(model, r[0].replace(tzinfo=timezone.utc), *r[1:]) for r in rows])
    conn.commit()
    return len(rows)


def write_emos(conn, rows, model="gefs"):
    """Upsert fits; rows = (station, lead_day, a, b, c, d, n, identity, train_from, train_to)."""
    with conn.cursor() as cur:
        cur.executemany(
            "INSERT INTO gefs.emos_fit (model, station, lead_day, a, b, c, d, n, identity, train_from, train_to) "
            "VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) "
            "ON CONFLICT (model, station, lead_day) DO UPDATE SET a = EXCLUDED.a, b = EXCLUDED.b, c = EXCLUDED.c, "
            "d = EXCLUDED.d, n = EXCLUDED.n, identity = EXCLUDED.identity, train_from = EXCLUDED.train_from, "
            "train_to = EXCLUDED.train_to, fitted_at = now()",
            [(model, *r) for r in rows])
    conn.commit()
    return len(rows)


def load_emos(conn, model="gefs"):
    """{(station, lead_day): (a, b, c, d)} plus the fit's metadata, or (None, None) if no fit."""
    with conn.cursor() as cur:
        cur.execute("SELECT station, lead_day, a, b, c, d, train_from, train_to, fitted_at, identity "
                    "FROM gefs.emos_fit WHERE model = %s", (model,))
        rows = cur.fetchall()
    if not rows:
        return None, None
    fits = {(st, ld): (a, b, c, d) for st, ld, a, b, c, d, *_ in rows}
    meta = {"train_from": min(r[6] for r in rows).isoformat(), "train_to": max(r[7] for r in rows).isoformat(),
            "fitted_at": max(r[8] for r in rows).strftime("%Y-%m-%d %H:%M UTC"),
            "identity": sorted({f"{st} d{ld}" for st, ld, *_, ident in rows if ident})}
    return fits, meta
