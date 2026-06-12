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
Rscript run_all.R
```

That restores pinned dependencies (`renv`), then runs `targets::tar_make()`.
**No network or API keys are required**: the default configuration
(`data.source: synthetic` in `config/config.yml`) simulates a block-exogenous
SOE data-generating process and runs the entire pipeline — estimation,
recursive evaluation, combination, figures, report — offline. Expect roughly
10–20 minutes on a laptop for the default settings.

### Real Australian data

Set in `config/config.yml`:

```yaml
data:
  source: real
```

Series are pulled key-free from the RBA (`readrba`), ABS (`readabs`) and
DBnomics (`rdbnomics`), cached under `data/raw/` with a manifest. If the
environment variable `FRED_API_KEY` is set, FRED is used as a fallback for the
US series; without it the pipeline still runs (DBnomics first) and any
irrecoverable download failure stops with a clear message suggesting the
synthetic path.

Two practical notes: the key-free US activity series (IMF IFS via DBnomics)
lags by about a year and the balanced panel is trimmed to the stalest series,
so a key-free real run ends roughly a year behind; setting `FRED_API_KEY`
recovers the missing quarters. A complete real-data run (June 2026, all §9
diagnostics green) is archived as `reports/report_real_data.html`.

## What it produces

| Where | What |
|---|---|
| `output/tables/` | scores by horizon (log score, CRPS, RMSE), DM tests, diagnostics, transform spec, combination weights |
| `output/figures/` | fan charts, PIT calibration histograms, weight-evolution plots, CRPS-by-horizon plots |
| `output/forecasts/` | combined + per-member point/interval forecast table |
| `reports/report.html` | rendered Quarto report (methodology + results narrative) |
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
