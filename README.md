# soe-bvar-suite

A reproducible R pipeline that forecasts the key Australian macroeconomic
aggregates 1–12 quarters ahead with a **suite of block-exogenous Bayesian
VARs** combined into pooled **density** forecasts, evaluated in a
pseudo-real-time recursive out-of-sample exercise.

Built for a small-open-economy workflow: every VAR partitions the system into
a foreign block (world activity, commodity prices, foreign policy rate) and a
domestic block, and block exogeneity — domestic variables never feed back into
the foreign block — is imposed in estimation and verified numerically.

## Quick start

```sh
git clone <repo>
cd soe-bvar-suite
echo 'FRED_API_KEY=your_key_here' > .Renviron   # free key: https://fred.stlouisfed.org/docs/api/api_key.html
Rscript run_all.R
```

That restores pinned dependencies (`renv`), then runs `targets::tar_make()`.
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
| `output/tables/` | scores by horizon (log score, CRPS, RMSE), DM tests, diagnostics, transform spec, combination weights |
| `output/figures/` | fan charts, PIT calibration histograms, weight-evolution plots, CRPS-by-horizon plots |
| `output/forecasts/` | combined + per-member point/interval forecast table |
| `reports/model_scorecard.md` | **the scorecard**: spec table for every model; forecast performance (CRPS / log score / RMSE, plus year-ended and cumulative-level) by variable and horizon; and a **per-model profile** (spec, role, strengths/failure modes, and where each model actually ranks) for every member |
| `reports/report.html` | rendered Quarto report (methodology + results narrative, fan charts, PIT calibration) |
| `DECISIONS.md` | the dated log of every modelling choice and its rationale |

## The suite

Seven block-exogenous VARs spanning prior family (Minnesota vs Villani
steady-state), size (8-var SOE core vs 13-var medium), tightness
(marginal-likelihood-selected λ vs fixed tight/loose), lag length and
volatility (constant vs stochastic), plus four univariate anchors (random
walk, AR(4), UCSV, unconditional mean). Combination schemes: equal weights,
recursive log-score weights with forgetting, the Hall–Mitchell/Geweke–Amisano
optimal pool (all shrunk toward equal, per variable × horizon bucket), and
BMA reported as a diagnostic.

## Correctness properties (enforced, not aspirational)

- **Block exogeneity binds** — max |posterior-mean domestic coefficient in a
  foreign equation| is logged per fit and gated; engines with structural
  zeros report exactly 0.
- **No look-ahead** — forecasts at origin *t* are produced from `td[1:t, ]`
  alone (including data-driven hyperparameter selection); an automated test
  corrupts all post-origin data and requires bit-identical forecasts.
  Combination weights at *t* use only scores realized by *t*.
- **Reproducibility** — one master seed, deterministically derived per task;
  same seed ⇒ identical draws, regardless of parallel scheduling.

## Layout

```
config/config.yml   every knob: variables, transforms, suite, MCMC, evaluation, combination
R/                  data_sources, transforms, priors, engines, forecast,
                    benchmarks, evaluate, combine, report, utils
_targets.R          pipeline DAG (targets); parallel over origins via future/furrr
tests/testthat/     unit tests (priors, block exogeneity, IW sampler, recursion,
                    scoring, pooling, optimiser, no-look-ahead, reproducibility)
cache/              per-(member, origin) OOS results keyed by config hash
```

Run the tests with `Rscript -e 'testthat::test_dir("tests/testthat")'`.
