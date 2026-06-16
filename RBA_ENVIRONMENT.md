# Running the SOE-BVAR suite in a locked-down environment (e.g. the RBA)

This note records which R packages the suite depends on, which ones are
typically **missing** from a controlled corporate environment such as the RBA
work setup, what breaks when they are, and the changes made so the suite runs
with the smallest possible set of new approvals.

It was written after auditing the project against the RBA's approved package
list (≈400 packages). The findings and the remediation below are current as of
the `run_all_plain.R` / logging-facade change.

---

## TL;DR

- **Before remediation:** the project would not even load. Five dependencies
  were missing, four of them `library()`-ed at the top of source files, so the
  moment `_targets.R` sourced the `R/` directory, sourcing died.
- **After remediation:** the only hard requirements that may need to be added
  to the approved mirror are **`scoringRules`, `stochvol`, `coda`**. Everything
  else the suite needs is already on the RBA list.
- Run it with **`Rscript run_all_plain.R`** (no `targets` required).

---

## 1. Dependency audit

Every package the code actually references, found by scanning `library()` /
`requireNamespace()` calls and `pkg::` usages across `R/`, `_targets.R`,
`run_all.R` and `reports/report.qmd`.

| Package | On RBA list? | Role |
|---|---|---|
| `scoringRules` | ❌ missing | CRPS + log-score scoring (core) |
| `stochvol` | ❌ missing | stochastic-volatility estimator (SV members, UCSV) |
| `coda` | ❌ missing | MCMC effective sample size (convergence gate) |
| `targets` | ❌ missing | pipeline DAG orchestration |
| `logger` | ❌ missing | logging |
| `rdbnomics` | ❌ missing | DBnomics pulls (`alt_foreign` set only) |
| `quarto` (R pkg) | ❌ missing | renders the optional HTML report |
| `readrba` | ✅ | RBA data pulls |
| `readabs` | ✅ | ABS data pulls |
| `fredr` | ✅ | FRED (foreign block) data pulls |
| `furrr` | ✅ | parallel map over forecast origins |
| `future` | ✅ | parallel backend |
| `ggplot2` | ✅ | figures |
| `yaml` | ✅ | parse `config/config.yml` |
| `digest` | ✅ | cache hashing + deterministic seeds |
| `glue` | ✅ | log-message interpolation |

---

## 2. The missing packages and their impact (pre-remediation)

Four of the five blockers are `library()`-ed at the *top* of a source file, so
they fail at **load time** — this is a hard stop, not a "some features
degrade" situation.

| Package | Powers | Where | Severity |
|---|---|---|---|
| `targets` | the entire pipeline orchestration (`tar_make`, `tar_read`) | `_targets.R`, `run_all.R` | 🔴 blocks the documented entry point |
| `scoringRules` | CRPS + log score — the whole scoring / evaluation / combination / scorecard layer | top of `evaluate.R`, `combine.R` | 🔴 core methodology; file won't source |
| `stochvol` | SV engine (`small_sv`, `medium_minn`) + `ucsv` benchmark | top of `engines.R`, `benchmarks.R` | 🔴 hard load-blocker |
| `coda` | MCMC ESS feeding the §9 convergence gate | top of `engines.R`, `benchmarks.R` | 🔴 hard load-blocker (same line as stochvol) |
| `logger` | all logging | top of `utils.R`, ~15 call sites | 🔴 hard load-blocker, but trivially replaceable |
| `rdbnomics` | DBnomics pulls — **only** the `alt_foreign` variable set | `data_sources.R` (lazy `::`) | 🟡 default small+medium run never calls it; doesn't block loading |
| `quarto` (R pkg) | renders `report.qmd` → `report.html` | `_targets.R` | 🟢 already wrapped in `tryCatch`; logs "render failed" and continues |

---

## 3. What the *required* packages do

After remediation the hard requirements are `scoringRules`, `stochvol`, `coda`
plus the data/compute stack already present.

### Statistically load-bearing (not easily replaceable)

- **`scoringRules`** — the forecast scoring engine. `crps_sample()` and
  `logs_sample()` grade every density forecast against the realized value
  (`evaluate.R`), and again for the combination pools (`combine.R`). The
  scorecard, the Diebold–Mariano tests and the combination weights are all
  built on these.
- **`stochvol`** — the stochastic-volatility estimator. `specify_priors()` +
  the `sv_*()` prior constructors and the fast C++ sampler
  `svsample_fast_cpp()` drive the SV members (`small_sv`, `medium_minn`) and
  the `ucsv` benchmark (`engines.R`, `benchmarks.R`). The
  `sv_infinity()`/`sv_exponential()` switch is the Gaussian-vs-t-error COVID
  treatment. (The constant-volatility engines roll their own linear algebra and
  don't need it.)
- **`coda`** — MCMC convergence diagnostics. `effectiveSize()` computes the
  posterior effective sample size (`engines.R`, `benchmarks.R`); `ess_min`
  feeds the §9 diagnostics table and `assert_diagnostics()`, which hard-stops
  the pipeline on inadequate convergence.

### Data acquisition (only for `data.source: real`)

- **`readrba`** (`read_rba()`), **`readabs`** (`read_abs()`), **`fredr`**
  (`fredr()` / `fredr_set_key()`) — the three providers behind the real panel,
  dispatched per-series in `data_sources.R`: RBA (cash rate, commodities, TWI),
  ABS (GDP, CPI, unemployment), FRED (US foreign block). The synthetic DGP path
  uses none of these.

### Parallel compute

- **`future`** + **`furrr`** — `future::plan(multisession)` spawns workers and
  `furrr::future_map(..., seed = TRUE)` fans the per-origin forecasts across
  them with reproducible RNG (`evaluate.R`). Removing them would only force a
  slow serial run.

### Plumbing & output

- **`yaml`** — parses `config/config.yml`, the single source of truth.
- **`digest`** — content hashing for the on-disk OOS cache key
  (`config_hash()`) and the deterministic per-task seed (`derive_seed()`).
- **`ggplot2`** — all figures in `report.R`.
- **`glue`** — interpolates `{var}` log templates (only in the logging
  fallback; degrades gracefully if absent).

---

## 4. Remediation applied

Two changes removed the avoidable blockers (`targets`, `logger`) so only the
genuine statistical dependencies remain.

### `logger` made optional — `R/aaa_logging.R`

A logging facade. Bare `log_info()` / `log_warn()` / … calls **forward to
`logger` when it is installed** (preserving the layout and namespace behaviour,
with the caller's environment forwarded via `.topenv` so glue templates still
resolve) and **fall back to base-R `message()`** in the same
`HH:MM:SS [LEVEL] msg` layout — with glue interpolation and a severity
threshold — when it is not. The `logger::` qualifiers were stripped across
`R/*.R` and `library(logger)` removed from `utils.R`. The file sorts first in
`R/` so the facade exists before any caller (including parallel workers). It is
deliberately **not** in `config_hash()`'s estimation-file list, because logging
is a pure side effect and must not invalidate the OOS cache.

> Tip: `options(soe.no_logger = TRUE)` forces the fallback even on a machine
> that has `logger`, which is how the fallback path was tested.

### `targets` made unnecessary — `run_all_plain.R`

A drop-in alternative to `run_all.R` that runs the exact `_targets.R` DAG
sequentially, calling the same project functions in dependency order. It fails
fast with a clear message if a hard dependency is missing, and the HTML-report
step is optional.

- **Kept:** the expensive recursive OOS step is still disk-cached by
  `run_oos_member()` (`cache/oos_<hash>/`), so unchanged estimation is skipped
  on re-runs.
- **Lost:** `targets`' object-level caching of the cheap downstream steps
  (scoring, combination, DM, scorecard, figures) — these recompute each run
  (seconds) — and `tar_read()` / DAG introspection.

---

## 5. How to run at the RBA

```sh
# one-off: install the three statistical packages from the approved mirror
#   install.packages(c("scoringRules", "stochvol", "coda"))
# (stochvol builds against Rcpp + RcppArmadillo, both already on the RBA list)

# real data needs a FRED key in a gitignored .Renviron:
echo 'FRED_API_KEY=your_key_here' > .Renviron

Rscript run_all_plain.R
```

For a fully offline run (no network, no FRED key) set `data.source: synthetic`
in `config/config.yml` — this exercises the machinery on a block-exogenous DGP
but not Australian dynamics.

The scorecard PDF is rendered separately and needs the **quarto CLI** (a
standalone binary, not the `quarto` R package):

```sh
bash reports/render_scorecard_pdf.sh
```

---

## 6. What still needs an install / approval

| Need | Required for | If unavailable |
|---|---|---|
| `scoringRules`, `stochvol`, `coda` | the core pipeline | the suite cannot run — these are not optional |
| `quarto` CLI | the scorecard `.pdf` | the scorecard `.md` is still written; render the PDF elsewhere |
| `quarto` R package + CLI | the HTML narrative report (`report.html`) | step is skipped with a warning; all other outputs produced |
| `rdbnomics` | the `alt_foreign` variable set only | the default small + medium suite does not use it |

`logger` and `targets` are **no longer required**.
