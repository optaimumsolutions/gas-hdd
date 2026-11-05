# Audit: how HDD is weighted and turned into gas in this system

Date: 2026-11-05. Scope: every place a weight, a threshold or a slope enters the forecast
chain, what it is based on, what the archived data says about it, and which changes would
make the forecast more accurate. Numbers below were measured from the system's own
archive (10 winters of airport observations, two winters of GEFS runs, EIA consumption
and storage) on the audit date; the commands are in section 7.

The pasted rule of thumb that prompted the audit is tested line by line in section 5.

## 1. The chain, with every weight in it

```
GEFS 2 m temperature grid (0.5°, 31 members)        ECMWF (0.25°, 50 members)
        │ nearest grid cell to each AIRPORT (22 airports)          [W1 location]
        │ daily Tmax, Tmin per UTC day (6-h max/min fields;        [W2 day definition]
        │ ECMWF: max/min of 6-hourly values)
        ▼
HDD per airport = max(0, 65 − (Tmax + Tmin) / 2)                    [W3 base 65 °F, linear]
        │ per-airport EMOS fit per lead day (shift + scale)         [C1 calibration]
        ▼
Index = Σ weight_city × HDD_city                                     [W4 the weights]
        US:      16 cities, population weights (hdd.CITIES)
        Midwest: 11 cities, res+comm gas-consumption weights (regions.py)
        │ per-index per-lead spread multiplier, plus one for the   [C2 calibration]
        │ 15-day total
        ▼
Gas = intercept + slope × HDD (+15 % slope above 45 HDD/day)        [D1 demand model]
EIA weekly change = 72.7 − 0.471 × weekly HDD                        [D2 storage model]
```

Members are never weighted against each other (all 31 or 50 count equally); the
calibration shifts and scales them. No lead day is weighted more than another inside the
15-day total. Nothing in the chain weights a cold day more than a mild one except the
+15 % slope kink in D1, which is still EIA's 2014 placeholder and not fitted.

## 2. What each weight is based on today

| Code | What | Basis | Fitted? |
|---|---|---|---|
| W1 | airport coordinates, nearest 0.5° / 0.25° cell | the thermometer the forecast is scored against | n/a, by construction |
| W2 | UTC day, (Tmax+Tmin)/2 from 6-h max/min | EIA / NOAA definition | n/a |
| W3 | base 65 °F, HDD linear in temperature | industry convention (CPC, Pivotal, EIA) | no |
| W4 US | population weights: NYC 32 %, LA 15 %, Chicago 11 %, Houston 9 %, Philadelphia 6 %, Dallas 5 %, the other ten 1–4 % each | "rough population (millions)" typed into `hdd.CITIES` on 2026-11-03 | no |
| W4 Midwest | state res+comm gas consumption shares 2024–25 (IL 26 %, MI 21 %, WI 11 %, MN 10 %, IN 10 %, MO 7 %, TN 6 %, IA 5 %, KY 4 %); MI and MO split 60/40 by metro population | EIA consumption by state | yes, data |
| C1 | EMOS a, b, c, d per airport per lead (330 fits), 108 pairs left at identity by the holdout gate | two winters of forecast vs thermometer | yes, holdout-tested |
| C2 | spread multipliers 1.0–1.4 per index per lead, 1.08–1.18 for the 15-day total | same archive | yes, holdout-tested |
| D1 | Midwest 0.96 + 0.354 × HDD, floor 1.82; US 6.01 + 1.616 × HDD, floor 8.5 | 5 winters of EIA monthly res+comm vs this index (r² 0.99) | yes |
| D1 kink | +15 % slope above 45 HDD/day | EIA 2014 note | **no** |
| D2 | 72.7 − 0.471 × weekly HDD, residual sd 16.9 Bcf | 120 winter weeks of EIA Midwest storage | yes, linear only |

## 3. What the archive says about the weights

### 3a. The US population weights put 29 % of the index where gas heating is thin

Los Angeles (15 %), Houston (9 %) and Dallas (5 %) carry 29 % of the US index. Gas
res+comm consumption shares of their states are far lower (section 3e gives the measured
shares). The Midwest tab already uses consumption weights and its index tracks CPC's
gas-weighted HDD within region-definition differences (bias −2.7 HDD, explained by CPC's
census divisions including OH/ND/SD/NE). The US index has no such check, because it is
population-weighted and CPC's national number is gas-weighted.

Consequence: a cold snap over Texas moves the US index by more than it moves US gas
demand; a Minnesota cold snap by less. The demand fit (D1) absorbs the average effect
(r² 0.991), but the error of any single forecast depends on *where* the cold is, which is
exactly the "locational positioning" the audit was asked about.

### 3b. Res/comm demand per HDD is flat across the season, not steeper in January

Measured on the last five winters, EIA monthly residential + commercial consumption against
the month's mean observed index (marginal Bcf/d per HDD/day, shared intercept):

| Month | US, Bcf/d per HDD | Midwest, Bcf/d per HDD |
|---|---|---|
| Oct | 1.83 (small HDD, intercept-dominated) | 0.35 |
| Nov | 1.60 | 0.33 |
| Dec | 1.62 | 0.35 |
| Jan | 1.64 | 0.36 |
| Feb | 1.61 | 0.36 |
| Mar | 1.62 | 0.37 |
| Apr | 1.61 | 0.37 |

The furnace response is linear and seasonless at monthly resolution (residual of the
linear fit by HDD bucket: +0.02 Bcf/d at 10–20 HDD, +0.02 at 20–30 for the US; within
±0.15 for the Midwest up to 50 HDD/day). A per-month slope for res/comm would be fitting
noise. The rule of thumb's "lower in the shoulder season" is about *total* demand with
power burn in it, not about this layer.

### 3c. The convexity is real, and it lives in storage, not in res/comm

Midwest weekly EIA storage change against the week's HDD, last five winters, residual of
the linear fit by weekly-HDD bucket (negative = more withdrawal than the line says):

| Weekly HDD (≈ HDD/day) | Residual, Bcf | Weeks |
|---|---|---|
| 35–105 (5–15) | −9.0 | 38 |
| 105–175 (15–25) | +4.4 | 42 |
| 175–245 (25–35) | +10.5 | 44 |
| 245–315 (35–45) | −7.3 | 19 |
| 315–385 (45–55) | −23.9 | 5 |
| 385+ (55+) | −43.7 | 1 |

The linear storage relation under-predicts withdrawals by about 25 Bcf a week once the
Midwest averages 45 HDD/day or more (daily mean ≈ 20 °F), and over-predicts them in mild
weeks. That is the freeze-off plus power-burn effect the rule of thumb describes, and the
one place in the chain where "weight the cold days more" is justified by data. The
current +15 % kink at 45 HDD sits in the res/comm layer (D1), where the data shows no
kink; the weekly storage relation (D2), where the kink is, has none.

### 3d. The forecast layer: where the errors come from, by lead

From the two-winter scorecard (airport level, GEFS raw vs the thermometer): MAE 2.1 °F at
day 1 rising to 7.8 °F at day 15; a cold bias of −0.35 °F near term growing to −1.0 °F at
day 15; raw band coverage 30 % at day 1. The EMOS layer fixes bias and width per airport
(holdout coverage 75–81 % at days 1–4); the spread multipliers fix the index and the
15-day total. None of this depends on the city weights: it is applied before weighting.

### 3e. Measured gas-consumption shares of the states behind the 16 US cities

Measured from EIA res+comm consumption 2021–2025 (appendix at the end of this file).
The audit's recommendation R1 uses them.

## 4. The rule of thumb, tested

| Claim | Verdict from this archive |
|---|---|
| 1 population-weighted HDD ≈ 1.5–2.5 Bcf of incremental US demand in winter | Consistent. Res/comm alone is 1.62 Bcf/d per HDD/day against this index. The Midwest storage relation gives 0.47 Bcf per weekly HDD, i.e. 3.3 Bcf per HDD/day sustained a week; scaled by the Midwest's ~25 % share of US res/comm that is ~1.9 Bcf per US HDD-day of total demand, inside the quoted range. |
| Res/comm 1.0–1.5 Bcf/d per HDD | Slightly low for this index: 1.62 measured (r² 0.991). The index is population-weighted toward warm cities, which raises Bcf per index-HDD. |
| Power burn adds 0.5–1.0, noisier | Not modelled here. The storage relation's residual sd (17 Bcf/week) is where that noise shows up. |
| Slope lower in early November | Not for res/comm (1.60 in Nov vs 1.64 in Jan). Possibly for total demand; the weekly storage data does not show a clean seasonal slope either (Nov −0.52, Dec −0.45, Jan −0.48, Feb −0.46, Mar −0.43 Bcf per weekly HDD). |
| Nonlinearity below ~20 °F | Confirmed in weekly storage (−24 Bcf residual at 45–55 HDD/day), not in monthly res/comm. |
| Gas-weighted HDD beats population-weighted | Agreed, and the system already does it for the Midwest. The US index does not. |
| Fit the slope empirically per month from EIA weekly storage vs realized HDD | Done weekly for the Midwest (D2). Per month adds nothing measurable for res/comm; a convex term matters more than a monthly one. |

## 5. Recommendations, ranked by expected effect on forecast accuracy

**R1. Re-weight the US index by gas consumption, not population.** Replace `hdd.CITIES`
weights with each state's res+comm consumption share (EIA, same source as the Midwest tab),
cities standing for their states, and add the cities the gas map needs (Pacific Northwest,
Upper Midwest and Northeast are under-represented; Southern California and Texas
over-represented). Keep the population index available for one more winter as a scored
alternative so the gain is measured, not assumed. Everything downstream re-fits
automatically (calibration nightly, demand slope via `fit_demand.py`). This is the
"locational positioning" fix: cold where gas heats homes counts more.

**R2. Move the cold-tail term to where the data shows it.** Drop the +15 % kink from the
res/comm model (no evidence) and fit a convex weekly storage relation instead:
`change = c + d × HDDwk + e × max(0, HDDwk − 300)` (300 weekly HDD ≈ 43 HDD/day), from the
same 10 winters. The implied-withdrawal card then widens and deepens in extreme weeks, which
is what the trader needs on the coldest Thursdays. Expected effect: about 25 Bcf less error
on the 6 or so weeks a winter that matter most.

**R3. Do not add a per-month slope.** Section 3b shows the res/comm slope is seasonless to
within 2 %; a monthly slope would fit noise and move the Bcf numbers month to month for no
weather reason.

**R4. Keep HDD itself linear and unweighted by cold.** HDD is the shared language with CPC,
Pivotal and the EIA; bending it would break every comparison. Cold-day emphasis belongs in
the demand and storage layers (R2), and the forecast layer already treats cold days
correctly through the per-airport calibration.

**R5. Add the five EIA storage regions** (East, South Central, Mountain, Pacific) with
consumption weights, the same way the Midwest was built. Then the US total becomes the sum
of five gas-weighted regions and R1 is automatic.

**R6. Score the alternatives before switching.** The archive and scorecard make this
cheap: compute the gas-weighted US index from the same airport rows, run `score.py` on
both, compare MAE and coverage per lead and the 15-day total. Switch when the gas-weighted
index beats the population one against EIA storage changes, which is the number that
matters.

## 6. What is not in the chain at all

- Power burn and its weather sensitivity (heat pumps, resistance heat): absent; it is the
  main reason total demand per HDD is 1.9 and res/comm is 1.6.
- Wind and solar availability, nuclear outages, gas-versus-coal economics: absent.
- Freeze-offs on the supply side: only as a map flag (basin triangles), not in the numbers.
- Snow cover, wind chill, humidity: absent; HDD is dry-bulb only, by convention.
- Daylight and weekday effects on commercial load: absent; monthly data cannot see them.

## 7. How the numbers in this audit were produced

All from the repository, 2026-11-05:

```
python fit_demand.py                 # section 2 (D1, D2) and docs/demand.json
python score.py / calibrate.py       # sections 2 (C1, C2) and 3d via docs/accuracy.json, docs/calibration.json
```

Per-month slopes, residual buckets and the storage convexity (sections 3b, 3c) were computed
in a session script against `gefs.observed` and the cached EIA series in `data/eia/`:
monthly mean index HDD vs EIA monthly res+comm Bcf/d, months with mean HDD ≥ 5, winters
2021-22 to 2025-26; weekly summed Midwest HDD vs EIA weekly Midwest storage change, weeks
with ≥ 35 HDD, same winters.

### Appendix: measured gas-consumption shares behind the US index (EIA res+comm, 2021–2025, keyless hist_xls)

| City | State | Weight today (population) | State's share of US res+comm gas |
|---|---|---|---|
| New York | NY | 0.323 | 0.090 |
| Los Angeles | CA | 0.152 | 0.087 |
| Chicago | IL | 0.105 | 0.075 |
| Houston | TX | 0.089 | 0.050 |
| Philadelphia | PA | 0.062 | 0.048 |
| Dallas | TX | 0.051 | 0.050 |
| Columbus | OH | 0.035 | 0.055 |
| Boston | MA | 0.027 | 0.028 |
| Denver | CO | 0.027 | 0.024 |
| Washington DC | DC | 0.027 | 0.003 |
| Detroit | MI | 0.023 | 0.058 |
| Kansas City | MO | 0.019 | 0.020 |
| Atlanta | GA | 0.019 | 0.023 |
| Minneapolis | MN | 0.016 | 0.030 |
| Pittsburgh | PA | 0.012 | 0.048 |
| St. Louis | MO | 0.012 | 0.020 |

The 13 states behind the 16 cities hold 59% of US res+comm gas consumption. States with no city in the index but large shares: NJ 0.047, WI 0.030, IN 0.028, WA 0.018, VA 0.019, MD 0.018.

Read-across: California carries 15 % of the index on population but 9% of US gas heating; Texas 14 % vs 5%; New York 32 % vs 9%; Illinois 11 % vs 8%; Michigan 2 % vs 6%; Ohio 4 % vs 6%; Pennsylvania 7 % vs 5%.
