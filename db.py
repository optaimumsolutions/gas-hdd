"""
db.py - the one place that talks to the Neon archive (PRD-calibration.md).

DATABASE_URL comes from the environment, or from a git-ignored .env file in the
repo root for local runs. The schema in sql/gefs.sql is applied on every
connect, so a fresh database needs no setup step. Every table name is written
schema-qualified; search_path is never set (Neon's pooler would leak it).
"""

import os
from pathlib import Path

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
