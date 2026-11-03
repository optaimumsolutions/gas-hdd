# gas-hdd — US natural-gas weather dashboard (GEFS)

**Live page:** https://optaimumsolutions.github.io/gas-hdd/

One Python script, one GitHub Actions workflow, one static web page. Four times a day it
downloads NOAA's GEFS ensemble forecast, turns it into US gas-weighted heating degree days
(HDD), draws a chart and three maps, and commits the results back to this repo. GitHub Pages
serves the `docs/` folder as the dashboard. No servers, no accounts, no paid data.

## What the dashboard shows and how to read it

The chart is the next 15 days of US gas-weighted heating demand. Each of the 31 GEFS
ensemble members is a faint grey line: one possible version of the weather. The bold blue
line is the median of those 31, and the shaded band runs from the 10th to the 90th
percentile. A narrow band means the members agree and the forecast is trustworthy; a wide
band means they do not. The number that matters most is the 15-day total, shown above the
chart with how much it moved since the previous run and since the run 24 hours earlier.
Green arrow up means the forecast got colder (more gas burned); red arrow down means warmer.
The change is measured over the days both runs forecast, so a day falling off the front of
the 15-day window does not count as a "change".

The three maps show *where* the cold is and *how sure* the model is. Colour is the
ensemble-mean 2 m temperature averaged over days 1–5, 6–10 and 11–15. The bold black line is
32°F. Hatching marks places where fewer than 70% of the members agree on which side of 32°F
the period average falls: solid colour is a confident forecast, hatched is a coin flip. Black
dots are the 16 HDD cities. Triangles are the gas-producing basins: a triangle turns red when
30% or more of the members take that basin's coldest moment in the period below its
freeze-off threshold (20°F for Permian, Haynesville, Anadarko and Eagle Ford; 0°F for
Appalachia, DJ and Bakken). The label next to each triangle is that percentage.

## Regional tabs (Midwest, in Bcf)

The **US** tab is the original population-weighted index (legacy). Each regional tab, starting
with **Midwest**, works differently:

- **Region = EIA weekly-storage region**, not a Census region, so it can be checked against the
  Thursday storage report. Midwest = IL, IN, IA, KY, MI, MN, MO, TN, WI (Ohio is East, Kansas
  is South Central).
- **Cities are weighted by residential + commercial gas consumption** (EIA, 2024 + 2025), not
  population. The weights, their sources and the state table are in `regions.py`.
- **HDD is converted to gas demand**: `Bcf/d = max(floor, intercept + slope × weighted HDD)`.
  Midwest: slope 0.38 Bcf/d per HDD, floor 1.9 Bcf/d (summer baseload).
- **Member by member**: each of the 31 members is converted to Bcf/d per day and summed over
  15 days *before* percentiles are taken. That is why the 15-day median is not exactly the
  sum of the daily medians.
- **Two ranges.** *Raw* p10–p90 is the spread of the 31 members. *Calibrated* p10–p90 widens
  that for GEFS's known under-dispersion (deviations from the median stretched ×1.0 days 1–5,
  ×1.15 days 6–7, ×1.3 days 8–15) and adds 10% slope uncertainty, then takes median ±
  1.2816 σ. The calibrated range is the honest one.
- Regional deltas are in **median 15-day Bcf**, over the days both runs share, same rule as
  the US tab.

**Sanity check on the slope.** A cold January in the Midwest runs about 40 weighted HDD/day
→ 0.38 × 40 ≈ 15 Bcf/d. EIA monthly data for IL + MI alone in Jan 2025 averaged 7.3 Bcf/d,
and IL + MI is about 47% of the 9-state region → ≈ 15.5 Bcf/d. The slope is in the right
place.

**Adding a region** = one more dict in `regions.py` appended to `REGIONS`. The page builds
its tabs from `latest.json`, so `index.html` does not change.

### Verifying the regional numbers

1. **HDD vs NOAA CPC.** CPC publishes free daily *gas-weighted* HDD by census division at
   `ftp.cpc.ncep.noaa.gov/htdocs/degree_days/weighted/daily_data/<year>/UtilityGas.Heating.txt`
   (observed, so check a day after it has happened). Day-1 Midwest HDD here should be within
   about 10% of the consumption-weighted average of CPC's East North Central and West North
   Central values for the same date. Different cities and weights: shape and level should
   match, not the exact number. Do this for three days running before trusting the tab.
2. **Bcf vs EIA storage.** Each Thursday, sum last week's forecast median Midwest demand over
   the EIA storage week (Saturday–Friday) and compare with the reported Midwest implied
   withdrawal. Expect the withdrawal to run 10–20% *above* res/comm demand (storage also feeds
   power burn and industrial). Keep a log of both numbers in `history/validation.csv`.
3. **Chart vs Pivotal Weather.** Same check as the US tab, but expect the Midwest HDD to run
   *higher* than the national number, since the warm cities are gone.

### What is still a placeholder

The page says so until `refit` is set to `True` in `regions.py`.

- Slope 0.38 is scaled from a 2010–14 EIA fit (STEO supplement, Oct 2014), not fitted to this
  index. Plausible range 0.33–0.45.
- Spread-inflation factors (1.0 / 1.15 / 1.3) are a guess at GEFS under-dispersion by lead.
- MI 60/40 (Detroit / Grand Rapids) and MO 60/40 (St. Louis / Kansas City) splits are by
  metro population, not utility customer counts.
- Daily mean is the average of four 6-hourly values, not (Tmax + Tmin)/2. Switch when the
  Phase 3 re-fit uses observed data, so forecast and observed HDD share one definition.

Phase 3 (`validate.py`, not built yet) re-fits the slope against weekly EIA Midwest
withdrawals and the inflation factors against this repo's own forecast archive.

## One-time setup (all in the browser)

1. **Settings → Pages**: under "Build and deployment" choose *Deploy from a branch*, branch
   `main`, folder `/docs`, then Save. The page will be at
   https://optaimumsolutions.github.io/gas-hdd/ within a minute or two of the first commit.
2. **Actions → "GEFS HDD dashboard" → Run workflow → Run workflow**. The job takes about
   10–20 minutes depending on how busy NOAA's S3 bucket is. Green check = the page has
   been updated. Red X = open the job and read the last lines of the "Build the dashboard"
   step; the error names the step that failed. A one-off red run is not a disaster: the
   next scheduled run six hours later starts fresh.

After that it runs itself at 05:30, 11:30, 17:30 and 23:30 UTC (about 5.5 hours after each
GEFS cycle starts, which is when the full 31-member run is normally on S3). If the newest
run is not complete yet, the script falls back to the previous one and says so in a yellow
note at the top of the page.

## How to verify it is correct

Do this a few times over the first week before trusting any number.

1. **HDD chart vs Pivotal Weather.** Go to pivotalweather.com → Models → GEFS and open the
   population-weighted HDD chart for the same run time shown at the top of this dashboard.
   Compare the 15-day total and the shape of the curve. They will not match exactly (Pivotal
   uses a different, larger set of cities and weights) but they should be within about 10%
   and show cold snaps on the same days.
2. **Maps vs Tropical Tidbits.** Go to tropicaltidbits.com → Models → GEFS → 2 m temperature
   (ensemble mean), same run, and step through the days. The cold and warm areas, and the
   position of the freeze line, should look the same as the maps here.
3. **The log.** Every Actions run prints the summary table and the file counts. If the job
   says it processed 2015 files (31 members × 65 lead times) and the table has 15 days, the
   inputs were complete.

If this dashboard and those websites disagree by a lot, this script is wrong, not the
websites. Open an issue or fix `hdd.py`.

## How the numbers are made

Everything lives in `hdd.py`, with comments at each step.

- **Source**: GEFS 0.5° "atmos.5" files on the NOAA open-data S3 bucket, every member
  (control + 30 perturbed), lead times 0–384 h every 6 h. Only the `TMP:2 m above ground`
  field is downloaded from each file, using the `herbie` library.
- **Daily mean temperature** per city per member = average of the four 6-hourly values in
  that UTC calendar day, at the nearest 0.5° grid cell, converted Kelvin → °F.
- **HDD** = max(0, 65 − daily mean °F).
- **US HDD** per day per member = Σ over cities of (weight × city HDD), weights normalised to
  sum to 1. The cities and weights are the table at the top of `hdd.py`.
- **Median / p10 / p90** are taken across the 31 members for each day. The 15-day total is
  the median of each member's own 15-day sum (close to, but not exactly, the sum of the daily
  medians).
- **Map agreement** = fraction of members whose period-average temperature at a grid cell is
  below 32°F; "uncertain" (hatched) where that fraction is between 0.3 and 0.7.
- **Freeze-off risk** = fraction of members whose period *minimum* 6-hourly temperature at
  the basin's grid cell is below the basin's threshold; flagged red at ≥ 0.3.
- **Days**: the 15-day window is the first 15 UTC calendar days that have all four 6-hourly
  values. For the 00Z run that starts today; for other cycles it starts tomorrow.

## Files

```
hdd.py                      the whole pipeline
regions.py                  region definitions, city weights, HDD -> Bcf model (pure data)
requirements.txt            python packages
.github/workflows/run.yml   the schedule
docs/index.html             the page (GitHub Pages serves this folder)
docs/latest.json            run time, table numbers, deltas (the page reads this)
docs/hdd_forecast.png       the chart, overwritten each run
docs/map_d1-5.png, map_d6-10.png, map_d11-15.png
history/hdd_YYYYMMDD_HHz.csv   per-member daily HDDs, one file per run, never deleted
docs/midwest_hdd.png, midwest_bcf.png, midwest_map_d1-5.png, ...   the Midwest tab
history/midwest_hdd_YYYYMMDD_HHz.csv, midwest_bcf_...csv   per-member daily HDD and Bcf/d
```

## Known limitations

- **GEFS only.** No ECMWF ensemble yet, so you cannot see when the two big models disagree.
- **0.5° resolution** (about 35 miles). City temperatures come from the nearest grid cell,
  which is fine for a national demand number but not for any one city.
- **Population weights are rough** (US tab only) and cover 16 cities. Real gas-weighted HDD uses state
  gas-consumption weights; expect a steady offset versus commercial products.
- **UTC days, not local days.** A "day" runs 00Z–00Z (7 pm–7 pm Eastern), so a cold front
  that arrives in the evening can land on a different day than in US-local products.
- **Freeze-off uses 6-hourly snapshots**, so the true overnight minimum can be a few degrees
  lower than what the script sees.
- GitHub pauses scheduled workflows in a repo with no commits by a person for 60 days; if the
  page goes stale, open the Actions tab and re-enable the workflow.

## Running it on your own computer (optional)

```
pip install -r requirements.txt
python hdd.py                      # newest complete run
python hdd.py --run 2026-11-03T00  # a specific run
```

Downloads about 300 MB per run into `data/` (git-ignored) and takes 5–10 minutes.
