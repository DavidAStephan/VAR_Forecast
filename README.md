# soe-bvar-suite

A reproducible R pipeline that forecasts the key Australian macroeconomic
aggregates 1–12 quarters ahead with a **suite of block-exogenous Bayesian
VARs** combined into pooled **density** forecasts, evaluated in a
pseudo-real-time recursive out-of-sample exercise.

Built for a small-open-economy workflow: every VAR partitions the system into
a foreign block (world activity, commodity prices, foreign policy rate) and a
domestic block, and block exogeneity — domestic variables never feed back into
the foreign block — is imposed in estimation and verified numerically.

This is the single reference document for the project. The per-model
specifications and live forecast-performance tables are in the auto-generated
[`reports/model_scorecard.md`](reports/model_scorecard.md) (regenerated on every
run); the rendered narrative report with figures is
[`reports/report.html`](reports/report.html).

## Contents

- [Quick start](#quick-start)
- [What it produces](#what-it-produces)
- [The model suite](#the-model-suite)
- [Correctness properties](#correctness-properties)
- [Repository layout](#repository-layout)
- [Data sources and audit](#data-sources-and-audit)
- [Modelling decisions](#modelling-decisions) (D1–D18)

## Quick start

```sh
git clone <repo>
cd soe-bvar-suite
echo 'FRED_API_KEY=your_key_here' > .Renviron   # free key: https://fred.stlouisfed.org/docs/api/api_key.html
Rscript run_all.R
```

That restores pinned dependencies (`renv`), then runs `targets::tar_make()`.

**No `targets`?** Use the drop-in plain runner instead — same steps, same
outputs, no `targets` dependency (you lose only DAG-level caching of the cheap
downstream steps; the expensive OOS step stays disk-cached):

```sh
Rscript run_all_plain.R
```

This is the entry point for locked-down environments (e.g. the RBA work setup).
Its hard requirements are `scoringRules`, `stochvol`, `coda` plus the usual
data/compute stack (`readrba`, `readabs`, `fredr`, `furrr`, `future`, `ggplot2`,
`yaml`, `digest`, `glue`); it fails fast listing any that are missing.
`logger` is **optional** — logging runs through a facade (`R/aaa_logging.R`) that
uses `logger` when installed and falls back to base-R console messages (same
layout) when it is not. `rdbnomics` is only needed for the `alt_foreign`
variable set, and the `quarto` R package only for the optional HTML report.

**The default is real data** (`data.source: real` in `config/config.yml`),
pulling live observations through the latest complete quarter:

- **Australian** data, key-free: RBA via `readrba` (cash rate, trimmed-mean
  CPI, real TWI, 10y yield, commodity prices) and ABS via `readabs` (national
  accounts, labour force, WPI).
- **Foreign block** from **FRED** (free API key required, kept in a gitignored
  `.Renviron`): US real GDP (`GDPC1`) for world activity and the effective fed
  funds rate (`FEDFUNDS`).

Pulls are cached under `data/raw/` (gitignored) with a manifest recording the
provider, series ID, and last observation; the cache re-downloads automatically
when a series ID changes. The current balanced panel runs **1997Q4–2026Q1**
(114 quarters). A full run is ~10 min plus first-time downloads.

### Offline / CI fallback

For a fully offline run with no network or keys, set:

```yaml
data:
  source: synthetic
```

This simulates a block-exogenous SOE data-generating process (used by the test
suite and the block-exogeneity ground-truth diagnostic). Synthetic outputs are
clearly labelled and validate the machinery, not Australian dynamics; **no
synthetic data ever enters a real run** — a missing real series stops the
pipeline with an explicit error rather than substituting.

A complete real-data run (foreign block GDPC1/FEDFUNDS, panel to 2026Q1, all §9
diagnostics green) is archived as `reports/report_real_data.html`.

## What it produces

| Where | What |
|---|---|
| `reports/model_scorecard.md` (+ `.pdf`) | **the scorecard**: spec table for every model; forecast performance on the **level** error (CRPS / log score / RMSE by variable and horizon — cumulative level for GDP/inflation, the rate level for unemployment/cash rate) **and on quarterly growth** (each target's single-quarter outcome); and a **per-model profile** (spec, role, strengths/failure modes, and where each model actually ranks on both views) for every member |
| `reports/report.html` | rendered Quarto report (methodology + results narrative, fan charts, PIT calibration) |
| `output/tables/` | scores by horizon (log score, CRPS, RMSE), DM tests, diagnostics, transform spec, combination weights |
| `output/figures/` | fan charts, PIT calibration histograms, weight-evolution plots, CRPS-by-horizon plots |
| `output/forecasts/` | combined + per-member point/interval forecast table |

The complete record of every modelling choice and its rationale is the
[Modelling decisions](#modelling-decisions) section below.

## The model suite

Seven block-exogenous VARs spanning prior family (Minnesota vs Villani
steady-state), size (8-var SOE core vs 13-var medium), tightness
(marginal-likelihood-selected λ vs fixed tight/loose), lag length and
volatility (constant vs stochastic), plus four univariate anchors (random
walk, AR(4), UCSV, unconditional mean). Combination schemes: equal weights,
recursive log-score weights with forgetting, the Hall–Mitchell/Geweke–Amisano
optimal pool (all shrunk toward equal, per variable × horizon bucket), and
BMA reported as a diagnostic. Full per-model detail — specification, role,
strengths and failure modes, and where each model ranks in the evaluation — is
in [`reports/model_scorecard.md`](reports/model_scorecard.md) §5; the design
rationale for the roster is decision [D8](#d8-the-suite-roster).

## Correctness properties

Enforced and tested, not aspirational:

- **Block exogeneity binds** — max |posterior-mean domestic coefficient in a
  foreign equation| is logged per fit and gated; engines with structural
  zeros report exactly 0.
- **No look-ahead** — forecasts at origin *t* are produced from `td[1:t, ]`
  alone (including data-driven hyperparameter selection); an automated test
  corrupts all post-origin data and requires bit-identical forecasts.
  Combination weights at *t* use only scores realized by *t*.
- **Reproducibility** — one master seed, deterministically derived per task;
  same seed ⇒ identical draws, regardless of parallel scheduling.

## Repository layout

```
config/config.yml   every knob: variables, transforms, suite, MCMC, evaluation, combination
R/                  data_sources, transforms, priors, engines, forecast, benchmarks,
                    evaluate, combine, covid, report, scorecard, utils;
                    aaa_logging (logging facade -> logger optional)
_targets.R          pipeline DAG (targets); parallel over origins via future/furrr
run_all.R           clone-to-result wrapper (renv restore + targets::tar_make)
run_all_plain.R     same pipeline without targets (for environments lacking it)
tests/testthat/     unit tests (priors, block exogeneity, IW sampler, recursion,
                    scoring, pooling, optimiser, no-look-ahead, reproducibility)
cache/              per-(member, origin) OOS results keyed by config hash
```

Run the tests with `Rscript -e 'testthat::test_dir("tests/testthat")'`.

## Data sources and audit

Method: every configured series was pulled from its provider and its **metadata**
(series description, table, series type, unit, frequency, coverage) checked
against the intended economic concept; every dlog-transformed series was
additionally tested for residual quarterly seasonality (F-test of quarter
dummies on quarterly log changes); constructed series were round-trip tested.
Sections: verified series → corrections made → proxies and constructed data →
synthetic data inventory → simplifying assumptions.

**Production configuration (this revision).** `data.source: real` is the
default. The foreign block is sourced from **FRED** (key required, in a
gitignored `.Renviron`): `f_act` = US **real GDP** (GDPC1), `f_rate` = US
**effective fed funds** (FEDFUNDS). Australian data comes key-free from RBA
(`readrba`) and ABS (`readabs`); commodity prices from RBA. The balanced panel
runs **1997Q4–2026Q1** (114 quarters, all 13 series real). The synthetic
generator is retained as an explicit opt-in (`data.source: synthetic`) for
tests/CI and the block-exogeneity ground-truth diagnostic only.

### 1. Verified correct (provider metadata, as of this audit)

| Variable | Series | Verified metadata |
|---|---|---|
| `gdp_growth` | ABS A2304402X (5206.0 Tab 1) | "Gross domestic product: **Chain volume measures**", **Seasonally Adjusted**, $m, quarterly, 1959Q3– |
| `cons_growth` | ABS A2304081W (5206.0 Tab 2) | "Households; Final consumption expenditure", **chain volume**, **SA**, $m, quarterly |
| `tot` | ABS A2304200A (5206.0 Tab 1) | "Terms of trade: Index", **SA**, quarterly (see §5 note on residual Q1 pattern) |
| `unemp_rate` | ABS A84423050A (6202.0 Tab 1) | "Unemployment rate; Persons", **SA**, %, monthly 1978– |
| `emp_growth` | ABS A84423043C (6202.0 Tab 1) | "Employed total; Persons", **SA**, '000, monthly |
| `wpi_growth` | ABS **A2713849C** (6345.0 Tab 1) | "Total hourly rates of pay excluding bonuses; Australia; Private and Public; All industries", Quarterly Index, **SA** — *corrected in this audit, see §2* |
| `cpi_inflation` | RBA GCPIOCPMTMQP (G1) | "Consumer price index; **Trimmed mean**; Quarterly change (per cent)" — the RBA/ABS trimmed mean is computed from seasonally adjusted CPI components by construction |
| `cash_rate` | RBA FIRMMCRTD (F1) | "**Cash Rate Target**", daily |
| `rtwi` | RBA FRERTWI (F15) | "**Real** trade-weighted index … adjusted for relative consumer price levels", **quarterly native** |
| `bond10y` | RBA FCMYGBAG10D (F2) | "Australian Government 10 year bond" yield, interpolated, daily |
| `f_comm` | RBA GRCPAISDR (I2) | "Index of commodity prices; All items; **SDR**", monthly, 1982– (SDR strips the endogenous AUD — correct for the price-taker foreign block) |
| `f_rate` | **FRED FEDFUNDS** | "Effective Federal Funds Rate", %, monthly, NSA (no seasonal concept for a rate); used as level, delta=1 |
| `f_act` | **FRED GDPC1** | "Real Gross Domestic Product", **chain-volume**, **Seasonally Adjusted**, Bil. Chained 2017 $, **quarterly-native** (1947–); used as dlog (qtr % growth), delta=0 — the same real-GDP concept as the domestic `gdp_growth`, so the foreign/domestic activity comovement is consistent |

Empirical seasonality F-tests on quarterly log changes (real data): all clean
(p > 0.05) except `tot` (p = 0.003; see §5). `wpi_growth` was failing
(p ≈ 1e-10) before the §2 correction. `f_act` = GDPC1 is published already
seasonally adjusted and quarterly-native, so the seasonality test is trivially
clean and within-quarter aggregation is a no-op for it.

### 2. Corrections made

1. **Foreign block moved to FRED (this revision).** `f_act` changed from
   DBnomics IMF/IFS US **industrial production** to **FRED GDPC1 (US real
   GDP)** — a better world-activity proxy for a quarterly macro VAR: it is
   quarterly-native (no frequency conversion) and the *same concept* (real
   GDP) as the domestic activity variable, rather than the narrower
   industrial-production index. `f_rate` changed from DBnomics
   `FED/H15/RIFSPFF_N.M` to **FRED FEDFUNDS** (same effective-rate concept,
   primary not fallback, fresher). This removes the previous "panel ends
   2024Q4" staleness (the IFS series lagged ~a year): the panel now reaches
   2026Q1. DBnomics remains only for the parked `alt_foreign` trade-weighted
   variant.
2. **WPI was the wrong series type.** The config used `A2603609J`, the
   **Original (non-seasonally-adjusted)** quarterly WPI index — confirmed by
   ABS metadata (`series_type: Original`) and a glaring seasonal in its log
   changes (F-test p ≈ 1e-10). Replaced with `A2713849C` (identical concept —
   total hourly rates of pay ex bonuses, Australia, private & public, all
   industries — **Seasonally Adjusted**).
3. **Commodity and terms-of-trade IDs.** Commodity index `GRCPAISAD` (A$/bulk
   -spot, starts 2009, embeds the endogenous AUD) → **`GRCPAISDR`**
   (SDR-denominated world price, 1982–). ToT `A2303731T` (did not resolve) →
   verified SA index `A2304200A`.

### 3. Proxies and constructed series (real-data mode)

These are **real published data**, but not literally the named concept:

1. **`f_act` "World activity" is US real GDP (a US-only proxy).** GDPC1 is the
   right *concept* (real GDP, the same as the domestic activity variable) but
   is US-only, a deliberate, documented proxy (see D1 below): the domain
   brief flags that the US is an imperfect proxy for Australia's Asia-weighted
   trading partners. The intended trade-weighted alternative (`f_act_tw`, OECD
   G20 GDP via DBnomics) is **stale on DBnomics (ends 2023Q3)** and is parked
   in the unused `alt_foreign` set, not the default suite — no fresh, machine
   -readable trade-weighted partner-GDP series was found (the RBA does not
   publish its trading-partner GDP aggregate as a statistical-table series).
   This US-only-world proxy is the principal remaining modelling simplification
   and is acceptable for production but flagged.
2. **`cpi_inflation` index level is constructed.** The RBA publishes the
   trimmed mean as a quarterly % change; the data layer cumulates it to an
   index (`100·cumprod(1+q/100)`) purely so the uniform dlog transform
   applies. Round-trip error ≈ 8e-14 (exact to float precision). The model
   forecasts log changes, ~0.2bp below the published arithmetic % change at
   typical inflation rates — negligible and consistent across all dlog series.
3. **Frontier and publication lags.** With the foreign block on FRED (GDPC1
   published ~1 month after the quarter; FEDFUNDS ~1 month) and ABS/RBA at
   their normal cadence, the balanced panel reaches **2026Q1** — the binding
   (stalest) series are the quarterly-native ones (ABS national accounts, the
   RBA trimmed mean, FRERTWI, GDPC1), all of which have a 2026Q1 print. The
   `end: 2026-03-01` config pin and the coverage-aware quarterly aggregation
   (§5.7) together ensure no partial (sub-quarter) frontier observation enters
   the panel. The old key-free "ends 2024Q4 / a year behind" limitation no
   longer applies.

### 4. Synthetic data inventory

1. **`data.source: synthetic` (an explicit opt-in, no longer the default)**
   simulates the *entire* panel from a block-exogenous DGP. It is retained so
   the pipeline can run fully offline (tests/CI) and so the block-exogeneity
   diagnostic has a known ground truth — it is still exercised by the test
   suite (test-pipeline.R, test-covid.R). It is clearly labelled in every
   output (the report states the data source and warns that synthetic results
   validate machinery, not Australian dynamics). **No synthetic values are
   ever mixed into the real-data path** — `get_raw_data()` branches strictly
   on `data.source`, and if a real series cannot be obtained the pipeline
   stops with an explicit error rather than substituting.
2. **Steady-state prior anchors** (`ss_mean`: inflation 2.5 % p.a. target
   midpoint, potential growth 2.8 %, NAIRU 4.5 %, neutral cash rate 3.5 %)
   are judgment numbers entering through priors, not data; sd's are
   configured and documented (see D5 below).
3. **External forecasts hook** (`read_external_forecasts`) returns NULL
   unless the user supplies a real file — no placeholder data.

### 5. Simplifying assumptions (accepted and documented)

1. **Quarterly aggregation = within-quarter mean** for daily/monthly series
   (cash rate, bond yield, fed funds, commodity index, unemployment,
   employment). Quarterly-average rates are the standard VAR convention;
   end-of-quarter is the main alternative and would slightly change rate
   dynamics. Native-quarterly series (US GDP/GDPC1, AU GDP/consumption/ToT,
   trimmed-mean CPI, WPI, real TWI) are unaffected.
2. **Cash rate = the target, not the realized interbank rate.** Identical in
   practice over the modelled sample (1997Q4–), where the target is announced
   and the interbank rate tracks it within a basis point.
3. **Final-vintage data, not real-time vintages.** The pseudo-real-time
   evaluation uses today's published history at every origin. Revisions to
   GDP/employment are material in Australia; rankings sensitive to revisions
   (especially GDP at h=1–2) should be read with care. Stated loudly in the
   report; building a vintage database is the natural extension.
4. **Balanced-panel trim.** The common sample is the intersection across all
   active variables: it starts **1997Q4** (the first dlog of WPI, whose level
   starts 1997Q3 — the binding series; GDPC1 from 1947 and FEDFUNDS from 1954
   do not bind) and ends at the stalest series (2026Q1, see §3.3). Alternative:
   drop WPI from the medium set and gain ~7 years of history — rejected for now
   because wages are central to the inflation block, but it is a one-line
   config change.
5. **`tot` residual Q1 pattern (p = 0.003).** The series is the ABS SA index;
   the pattern (Q1 mean +1.5pp vs ~0, sd 3.4, lag-1 autocorr 0.21) reflects
   bulk-commodity contract repricing historically clustered in Q1 — a price
   phenomenon, not a removable seasonal in the SA sense, and the simple
   F-test overstates significance under autocorrelation. Kept as is.
6. **The default `start: 1990` is aspirational for the real panel** — the
   actual start is determined by the balanced trim (1997Q3). Config start
   only truncates from the left.
7. **No seasonal adjustment is performed in-pipeline.** All series are sourced
   already-SA (or are concepts that need no SA: rates, SDR commodity prices,
   real TWI). If a future series is only available Original, adjust at source
   or add an SA step — do not dlog an Original series (the WPI bug in §2 is
   the cautionary tale).
8. **Coverage-aware quarterly aggregation (this revision).** `to_quarterly()`
   drops (sets NA) any quarter with materially fewer source observations than
   a full quarter (< 80% of the fullest quarter's count), so a partial
   frontier quarter becomes NA and is trimmed rather than entering as a
   1–2-month average. Combined with the loud interior-NA error in
   `transform_data()` and the pct_change interior-NA guard, the real path
   never silently substitutes or shortens: a gap either trims cleanly
   (trailing) or errors (interior). The cache is keyed by provider+series_id
   so a changed series ID re-downloads rather than serving stale data.

## Modelling decisions

Running log, in the spirit of "decide, document, proceed". Each entry: the
choice, the rationale, the rejected alternative, and the config key that
controls it.

### D1. Variable roster and blocks

**Choice.** Foreign block: world activity (`f_act`, **US real GDP** — FRED
GDPC1; see D18), commodity prices (`f_comm`, RBA Index of Commodity Prices,
SDR), foreign policy rate (`f_rate`, US effective fed funds — FRED FEDFUNDS),
plus an alternative trade-weighted activity measure (`f_act_tw`, G20 GDP)
reserved for an alt-foreign variant. Domestic small core: GDP growth,
trimmed-mean CPI inflation, unemployment rate, cash rate, real TWI. Medium set
adds terms of trade, WPI growth, employment growth, consumption growth,
10-year AGS yield (13 variables total).

**Why.** The small core covers the four policy targets plus the exchange rate
(the SOE adjustment channel). The medium set follows Bańbura–Giannone–Reichlin
(2010): medium systems often forecast better with enough shrinkage. The US
proxy is the standard caveat for Australia (Asia matters more now); the
trade-weighted alternative is wired into the spec (`sets: [alt_foreign]`) as a
documented suite-diversity hook.

**Rejected.** A large (20+) system — recursive re-estimation cost grows
quadratically and the medium set already exercises the shrinkage machinery.
**Config:** `variables`, each member's `set`.

### D2. Transforms and Minnesota deltas

**Choice.** Growth-type series enter as quarterly 100·Δlog with `delta = 0`;
interest/unemployment rates as levels with `delta = 1`; real TWI and ToT as
100·log levels with `delta = 1`. Trimmed-mean CPI (published as a quarterly %
change) is cumulated to an index in the data layer so all series arrive as raw
levels and one transform layer serves both data paths.

**Why.** `delta = 0` for mean-reverting growth (white-noise prior centre) and
`delta = 1` for persistent levels is the textbook Minnesota assignment
(Litterman 1986); modelling rates in levels keeps the cash rate and the policy
question in interpretable units, with the I(1)-discipline handled by priors
(D5) rather than differencing away level information.

**Rejected.** Differencing everything (loses cointegration-ish level
relations); modelling log CPI level (long-memory inflation drifts wreck a
fixed-mean VAR over a 35-year window). **Config:** `variables.*.transform`,
`variables.*.delta`.

### D3. Block exogeneity as a first-class property

**Choice.** Every VAR member imposes the foreign/domestic partition
(Cushman–Zha 1997; Zha 1999): foreign equations exclude domestic information.
Three mechanisms, by engine: (i) `gibbs`/`ss`: near-zero prior mean **and**
variance (`mcmc.block_exog_prior_sd: 1e-4`) on domestic-lag coefficients in
foreign equations — possible only because these engines use an independent
(non-Kronecker) prior; (ii) `conj_br`: exact, by estimating the foreign block
as its own VAR (block-recursive, the Bloor–Matheson RBNZ approach) with the
domestic block conditioned on contemporaneous foreign values; (iii) `sv`:
exact, by dropping domestic columns from foreign equations in the
triangularised system. The max |posterior-mean domestic coefficient in a
foreign equation| is logged per fit, checked against synthetic ground truth,
and gated in `assert_diagnostics()`.

**Why.** Without it the model implies Australia moves world demand and
commodity prices, which contaminates medium-horizon dynamics; RBA RDP 2013-06
shows the restriction materially raises the role of foreign shocks.

**Scope — dynamics, not the steady-state means.** Block exogeneity constrains
the **lag dynamics** (the autoregressive coefficient matrix: foreign equations
carry no domestic lags). It does *not* make the steady-state *means*
block-exogenous: in the `ss` (Villani 2009) member the unconditional means ψ are
drawn from their joint GLS full conditional using the full, non-block-diagonal
error covariance Σ, so the foreign and domestic means are statistically coupled
through the contemporaneous error correlations. This is the correct posterior —
a foreign-only ψ update would be an *incorrect* full conditional — so the
foreign steady states the domestic forecast inherits are informed by the joint
residual covariance, by design (see the note in `fit_ss`, R/engines.R). The
block-exogeneity diagnostic in `assert_diagnostics()` therefore (correctly)
checks the lag coefficients, not ψ.

**Rejected.** A pure Kronecker-conjugate suite (cannot represent asymmetric
prior variances — this is exactly why the spec requires the Gibbs engine).
**Config:** `variables.*.block`, `mcmc.block_exog_prior_sd`.

### D4. Tightness: GLP-style marginal-likelihood selection

**Choice.** Overall tightness λ is selected by maximising the closed-form
conjugate NIW marginal likelihood over a grid (`glp.lambda_grid`), re-run *at
every forecast origin inside the recursive loop*, and shared by the members
whose `prior.lambda` is `auto`. Two members keep fixed λ (0.05 tight, 0.4
loose) as suite diversity.

**Why.** Giannone–Lenza–Primiceri (2015): treat shrinkage hyperparameters as
parameters, not folklore. The conjugate marginal likelihood is cheap and
closed-form; a full hierarchical MCMC over λ inside every origin × member
would multiply compute for little forecast gain. Re-selecting per origin keeps
the no-look-ahead property exact.

**Rejected.** Hand-set λ everywhere (indefensible); full hierarchical NUTS/MH
over hyperparameters (cost). **Config:** `glp.*`, member `prior.lambda`.

### D5. Long-run priors: SOC + DIO, and a steady-state member

**Choice.** Sum-of-coefficients (Doan–Litterman–Sims) and
dummy-initial-observation (Sims–Zha) artificial observations on the
`delta = 1` level variables in selected members
(`prior.soc/dio`); and a Villani (2009) steady-state member with informative
priors on the unconditional means: trimmed-mean inflation 2.5 % p.a. (target
midpoint, sd 0.6pp), GDP growth 2.8 % p.a., NAIRU 4.5 %, neutral cash rate
3.5 %, foreign anchors likewise; loose sample-mean anchors for the TWI/ToT.

**Why.** Iterated forecasts at h = 8–12 converge to the model's unconditional
mean, which plain Minnesota pins down poorly; the steady-state
reparametrisation lets policy-relevant long-run information (and the foreign
block's steady states, which the domestic forecast inherits) enter directly.

**Rejected.** SOC/DIO on growth-rate variables (their sample means are already
well-estimated; the dummies would just bias short-run dynamics).
**Config:** member `prior.soc/dio/soc_mu/dio_delta`, `variables.*.ss_mean/ss_sd`.

### D6. Stochastic volatility via exact triangular equation-by-equation

**Choice.** The SV engine factorises the reduced form as a recursive system
(foreign block ordered first), estimating each conditional equation
independently: Bayesian regression with SV errors, β drawn by weighted least
squares given the log-variance path, the path updated with
`stochvol::svsample_fast_cpp`. This is the Carriero–Clark–Marcellino (2019)
equation-by-equation idea; with independent per-equation priors the posterior
factorises *exactly*, so per-equation estimation is not an approximation.

**Why.** Keeps the recursive OOS loop feasible (the medium SV member estimates
in seconds) while delivering the density-calibration benefits of SV
(Clark 2011) and protecting against the 2020 outliers.

**Rejected.** Joint VAR-SV with a full covariance Gibbs (cost explodes in the
recursive loop); Lenza–Primiceri outlier dummies (SV with the COVID quarter in
the synthetic DGP already absorbs it; noted as an extension).
**Config:** members with `engine: sv`.

### D7. Iterated multi-step forecasts

**Choice.** All members iterate the one-step model forward, simulating shocks
and parameter draws jointly (h = 1…12); no direct/local-projection members.

**Why.** Marcellino–Stock–Watson (2006): iterated wins at long horizons when
the one-step model is adequate, and only iterated forecasts give internally
coherent joint densities for year-ended transforms and fan charts.

**Rejected.** A direct-forecast robustness member — 12 extra estimations per
member-origin for a known-weaker density forecaster; revisit if point RMSE at
h ≥ 8 looks broken. **Config:** none (structural).

### D8. The suite roster

Seven VARs + four univariate anchors (see `suite` in config):
`small_minn` (Gibbs, GLP-λ, SOC+DIO), `small_ss` (steady-state),
`small_sv` (SV), `small_loose_p5` (λ=0.4, p=5), `medium_minn` (13-var SV,
p=2), `medium_conj` (13-var block-recursive conjugate, SOC+DIO),
`small_tight` (λ=0.05 conjugate), plus `rw`, `ar4`, `ucsv`, `ucmean`.
Diversity axes per the brief: prior family, size, tightness, lag length,
volatility model. The members genuinely fail differently: the tight small
model is robust at long horizons, the medium SV member wins short-horizon
density, the steady-state member anchors h ≥ 8, the benchmarks discipline
everything.

**Rejected.** An explicit US-proxy vs trade-weighted pair of otherwise
identical members in the default roster (the alt foreign series is a DBnomics
pull that fails without network; the variant is wired in config via the
`alt_foreign` set but not in the default suite, so the offline run stays
green). **Config:** `suite`, `benchmarks`.

### D9. UCSV applied per target variable, outlier-robust specification

**Choice.** A UC-SV model (local level + FFBS) is run for *each* target
variable, not just inflation, with three deliberate deviations from the
canonical twin-SV Stock–Watson form, all forced by COVID-sized outliers in the
evaluation window:
(i) **t-distributed transitory errors** (`nu ~ Exp(0.1)` via stochvol) — the
Gaussian version goes sticky (trend ESS ≈ 2) at origins where the outlier is
the sample endpoint;
(ii) **constant trend-shock variance with a conjugate IG update** instead of a
second SV process — the trend-vol dimension is where the endpoint
trend-vs-noise bimodality piles up, and a direct Gibbs draw removes it;
(iii) **tight vol-of-vol prior** (`sigma2 ~ Gamma(2, 50)`, mean 0.04) — the
Stock–Watson tradition *fixes* the vol smoothing parameter near 0.2; a loose
prior lets an endpoint outlier inflate vol-of-vol and the 12-step-ahead
variance forecast explodes through exp(h/2) compounding (observed: 99.5%
unemployment quantiles in the hundreds). Updated 2026-06-12 after the first
full synthetic run failed the convergence gate at origins 117–120.

**Why per target.** It is the canonical inflation benchmark, and as a generic
trend-plus-noise-with-SV model it is a sensible univariate density anchor for
the other targets too. **Config:** `benchmarks`, `mcmc.bench_*`.

**Convergence gating.** The trend/noise *split* of a near-white series is
weakly identified in any UC model (the classic variance pile-up): chains for
the decomposition parameters (trend variance, SV level) and for the trend
endpoint mix slowly even in very long runs, while the predictive *sum* —
the only thing the forecasting suite consumes — is well identified. UCSV
convergence is therefore gated on a Monte-Carlo-precision criterion for the
forecast-relevant quantity: the MCSE of the trend endpoint must be below 15 %
of the one-step predictive sd, with an adaptive triple-thinning retry before
failure. Raw ESS values are still computed and reported (trend-endpoint ESS
appears as `ess_min` in the diagnostics table; it can be low at origins where
an outlier sits at the boundary, which is expected and documented rather than
hidden).

### D10. Combination: linear pools, per variable × horizon bucket, shrunk

**Choice.** Density-level linear pools with weights per target variable and
per horizon bucket (1–4, 5–8, 9–12): equal; recursive log-score weights
(Jore–Mitchell–Vahey 2010) with forgetting 0.95; optimal pool
(Hall–Mitchell / Geweke–Amisano) maximising the discounted historical pool log
score, softmax-parametrised BFGS; BMA via the shortest-horizon-in-bucket
predictive likelihood (h = 1 / 5 / 9 for near / medium / far),
**reported as a diagnostic only**. Non-BMA weights are shrunk toward equal
with κ = 0.3; weights at origin t use only scores with s + h ≤ t; equal
weights until 8 training origins exist.

**Why.** Bucketing by horizon recognises that the best member changes with
horizon (the design premise of the suite); shrinkage + forgetting are the
standard guards against weight overfitting (the combination puzzle — equal
weights are the benchmark to beat, so they are also a scheme). Geweke–Amisano:
the optimal pool does not degenerate to one model; BMA does, which is why it
is reported but not recommended. Pool log scores use log-sum-exp; pool CRPS
uses draws resampled from the mixture; pool PIT is the weighted member PIT
(exact mixture CDF).

**Rejected as default.** Del Negro–Hasegawa–Schorfheide dynamic pools — the
forgetting-factor log-score scheme already provides time variation; a full
dynamic pool is the natural next extension. **Config:** `combination.*`.

### D11. Evaluation design

**Choice.** Expanding window (rolling available via `evaluation.window`),
origins = the last `max_origins` (36) quarters ending at T−1, final-vintage
data with the caveat stated loudly in the report. Scoring: RMSE, log score
(kernel density from draws), CRPS, PIT, by horizon and variable, quarterly and
year-ended. DM tests with Harvey correction and Newey–West (h−1) variance vs
the AR(4) and RW anchors.

**Why.** Expanding windows match how a central bank actually re-estimates;
36 origins balances test power against compute. Real-time vintages for
Australia are not reliably available key-free; pretending otherwise would be
worse than documenting the caveat. DM rather than Giacomini–White: with an
expanding window and these sample sizes the conditional GW adds machinery
without changing conclusions; noted as an extension.
**Config:** `evaluation.*`.

### D12. Synthetic DGP

**Choice.** A stationary block-exogenous VAR(2) in the transformed units with
a structured (not random) coefficient pattern — foreign→domestic activity
loadings, a Taylor-rule-ish cash-rate equation, monetary transmission with the
right signs — Gaussian shocks whose contemporaneous loadings run strictly
foreign→domestic, a COVID-style outlier quarter with partial rebound, and
integration back to raw levels so the transform layer runs identically to the
real path.

**Why.** Random off-diagonal coefficients compound badly under high
persistence (rate series wandered to ±10 %); a structured DGP keeps every
series in realistic ranges, gives the block-exogeneity diagnostic a known
ground truth, and exercises the same code paths as real data.
**Config:** `synthetic.*`, `data.source`.

### D13. Seeds and parallelism

**Choice.** One master seed; every stochastic task (member × origin, final
forecasts, mixture resampling) derives its own seed by hashing the master seed
with a task key, so results are identical regardless of worker scheduling.
Parallelism is over origins inside each member target (`furrr`), because the
recursive loop, not any single fit, is the bottleneck.

**Rejected.** future's built-in L'Ecuyer streams alone (reproducible only for
a fixed worker count). **Config:** `master_seed`, `parallel.workers`.

### D14. Caching

**Choice.** Raw downloads cached in `data/raw/` with a manifest and freshness
window; per-(member, origin) OOS results cached in `cache/oos_<confighash>/`
keyed by a hash of the estimation-relevant config sections, so combination
experiments re-use estimation; `targets` provides DAG-level caching on top.

**Why.** The recursive loop is the expensive stage; weight-scheme tweaks should
cost seconds, not an hour. **Config:** `data.cache_max_age_days`; hash covers
`master_seed, data, synthetic, variables, horizons, mcmc, glp, suite,
benchmarks, evaluation`.

### D15. Conditional forecasts

**Choice.** Hard conditioning by substitution (overwrite the conditioned
variable's path each step as it feeds the recursion), exposed as an optional
argument, documented as an approximation.

**Rejected.** Waggoner–Zha conditional simulation — correct but heavy, and
conditional forecasting is an optional capability in this spec, not a scored
deliverable.

### D16. Real-data details

**Choice.** Australian series via `readrba`/`readabs` (key-free); the US
foreign block via FRED (`fredr`, key required — see D18); commodity prices via
RBA. Mixed-frequency series are quarterly-averaged with a coverage guard (a
partial frontier quarter is set NA, not averaged from <3 months). Failures stop
with a message pointing at the synthetic path (the always-runnable guarantee
lives there, not in fragile retries or silent series substitution).

**Caveat recorded.** Trimmed-mean CPI is published as a quarterly rate and is
cumulated to an index purely so the uniform dlog transform reproduces it; the
round trip is exact up to float error, and an interior gap errors loudly rather
than NA-poisoning the cumulative product. The same `pre: pct_change` mechanism
handles any source published as a % change. The long-history RBA commodity
index is `GRCPAISDR` (from 1982); the bulk-spot variant `GRCPAISAD` only starts
in 2009 and silently truncated the panel — found and fixed earlier.

### D18. Production migration to real-first data (added 2026-06-13)

**Choice.** `data.source: real` is the **default** (was synthetic); the
foreign block is sourced from **FRED** — `f_act` = **US real GDP (GDPC1)**,
`f_rate` = **US effective fed funds (FEDFUNDS)** — and Australian data extends
to **2026Q1**. The FRED key lives in a gitignored `.Renviron`
(`FRED_API_KEY`). `fred` is a first-class provider in `fetch_series`; the
legacy DBnomics→FRED silent-fallback machinery was removed (it could swap an
unconfigured series ID); the on-disk cache is keyed by provider+series_id so a
changed ID re-downloads; quarterly aggregation is coverage-aware (DATA_AUDIT
§5.8). The synthetic generator stays as an explicit opt-in for tests/CI and
the block-exogeneity diagnostic.

**Why f_act = GDPC1 rather than US IP.** For a *quarterly* macro VAR whose
domestic block already contains Australian real GDP growth, US real GDP is the
conceptually consistent world-activity measure (same concept, quarterly-native,
no frequency conversion) — superior to industrial production (monthly, a
shrinking ~15–20% slice of US output). GDPC1 is chain-volume, seasonally
adjusted, published ~1 month after the quarter, so the panel reaches the
current frontier rather than lagging ~a year as the old key-free IMF/IFS US-IP
series did. The foreign activity steady-state anchor is correspondingly US real
GDP growth (`ss_mean` 0.58/qtr ≈ 2.3% p.a., the 1997Q4–2026Q1 realized mean),
distinct from the old IP-index growth anchor.

**Why FRED for the US block.** With a (free) key, FRED is the authoritative,
fresh, reliable US source; it removes the dependency on DBnomics (which was the
flaky link — the FED/G17 series failed and the IMF/IFS proxy was stale).
DBnomics remains only for the parked `alt_foreign` trade-weighted variant.

**Rejected.** Keeping synthetic as the default (the institutional
offline-runnable requirement is satisfied by retaining it as an opt-in, not as
the default for a production forecasting suite); INDPRO as a silent fallback
for GDPC1 (a different concept — a failure should error, not substitute);
a non-US trade-weighted world proxy in the default suite (no fresh,
machine-readable partner-GDP aggregate is available — the US-only proxy remains
the documented simplification). **Config:** `data.source`, `data.end`,
`data.fred_api_key_env`, `variables.f_act/f_rate.source`.

### D17. COVID-period estimation (added 2026-06-13 after the code/literature audit)

**Problem.** Three quarters of 2020 data (GDP −7 %, unemployment +3 pp, the
2020Q3 rebound) are 5–30 standard-deviation events. Untreated, they (i) wreck
constant-volatility BVAR coefficient estimates — Lenza & Primiceri (2022 JAE)
measure a 965-log-point marginal-likelihood gap and show the GLP tightness
posterior shifting from λ≈0.13 to λ≈0.5; (ii) contaminate the Minnesota scale
calibration itself, since σ-scalings from AR residual RMSDs are not
outlier-robust (Hartwig 2024 SNDE); and (iii) in plain persistent-SV models,
get misread as a *persistent* volatility shift that inflates predictive bands
for years (CCMM 2024 REStat; Hartwig).

**Choice (three coordinated mechanisms, config `covid:`):**
1. **Constant-volatility engines (`gibbs`, `ss`, `conj_br`) + AR/RW/mean
   benchmarks: Lenza–Primiceri volatility scaling.** y_t = c + B(L)y_{t−1} +
   s_t ε_t with a free scale per configured COVID quarter (default 2020Q1–Q3)
   and geometric decay of the excess sd afterwards. Implemented exactly as LP:
   rows of (Y, X) divided by s_t (GLS), which preserves every engine's
   machinery including the Kronecker conjugate; (s, ρ) estimated at *every
   forecast origin* by maximizing the closed-form conjugate marginal
   likelihood including the Jacobian −M·Σlog s_t, coordinate-descent on grids
   — the GLP hyperparameter treatment, fully recursive (no look-ahead).
   Predictive shocks at COVID-era origins carry the decaying s-path
   (LP's density-forecast payoff: dropping data "vastly underestimates
   uncertainty"). The Minnesota σ-calibration and the GLP λ-selection use the
   same weighted data (closing Hartwig's contamination channel).
2. **SV engine: t-distributed errors** (stochvol `nu ~ Exp(0.1)`, scale
   mixture in the WLS step, t-draws in prediction) instead of row weighting —
   the CCMM SV-t specification, which in their comparison performs within
   0–2 % of their preferred SVO-t and avoids misattributing the spike to the
   persistent volatility component. Config `covid.sv_t_errors`.
3. **`treatment: dummy` option**: scales → 10³ at the COVID quarters,
   equivalent to dummying-out/dropping those likelihood rows (Schorfheide–Song
   2024; Cascaldi-Garcia's Pandemic Priors at φ→0), with future shocks at
   s = 1. Provided as the comparison/robustness option; its documented cost is
   predictive bands stuck at pre-COVID widths (CCMM).

**What is deliberately NOT removed.** COVID observations remain in the
forecast jump-off (lagged regressors / conditioning information) — the
literature is explicit that only the likelihood contribution should be
treated (CCMM's AR(1) jump-off discussion; Schorfheide–Song filter through
the pandemic for the origin state).

**Evaluation interaction.** Mean log scores over windows containing 2020
realizations are dominated by a few quarters; the evaluation now also writes
score tables excluding the COVID realization dates, and CRPS (less
outlier-sensitive) is reported alongside log scores throughout.

**Rejected.** Full CCMM SVO (outlier mixture states) — SV-t achieves nearly
the same per their own results at a fraction of the implementation; Ng's
(2021) epidemiological covariates — wrong data dependencies for an
always-runnable pipeline; sample truncation at 2019Q4 (RBA RDP practice) —
forfeits 5 years of data and the post-COVID inflation episode, the most
informative period in the evaluation window.

**Validation.** Unit tests: ML detects the synthetic outlier (scales ≥ 8);
LP-weighted coefficient posteriors are materially closer to clean-sample
posteriors than untreated ones on synthetic ground truth; the dummy limit is
insensitive to the exact tiny weight; no-look-ahead holds with the treatment
active. **Config:** `covid.*`.

**Empirical result (real Australian data, 2026-06-13).** Recursive ML
estimates of the scales stabilise at (1, 8, 4) for 2020Q1–Q3 — 2020Q1 was
normal in Australia, the residual sd in 2020Q2 was ~8× normal, 2020Q3 ~4×.
Comparing `lp_scaling` vs `none` over the same 36-origin evaluation
(`reports/covid_treatment_comparison.csv`): at COVID-era forecast origins
(2020Q1–2022Q4) CRPS improves for *every* treated member — ratios 0.89
(loose member) to 0.97, GDP-growth CRPS −7 % — and the pooled combinations
gain 0.16–0.25 nats of mean log predictive density; at all other origins the
treatment is a no-op (ratios 0.98–1.00), confirming the adjustment only
binds where it should (Álvarez–Odendahl's pre-specified-dates principle).
Two constant-volatility members lose a little log score at COVID origins
while gaining CRPS — wider LP tails cost log density when the realization is
central; this is the known CRPS/log-score divergence on outlier windows and
is why both are reported.

