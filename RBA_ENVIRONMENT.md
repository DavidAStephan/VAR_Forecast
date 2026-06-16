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
- **After remediation:** the suite **runs even with `scoringRules`, `stochvol`
  and `coda` all absent** — it degrades gracefully (see §4.3) rather than
  failing. Everything strictly needed to run is already on the RBA list. You
  only need those three to get the *full* output (scoring/scorecard, the SV
  members, and the convergence gate respectively).
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

After remediation `scoringRules`, `stochvol` and `coda` are **optional** — the
suite degrades gracefully when they are absent (§4.3) — leaving only the
data/compute stack (already present) as strictly required. You still want all
three installed for the full output. The table below is the quick
reference; the prose after it adds detail, and §3.1 spells out the
estimation-vs-evaluation distinction for the three statistical packages.

| Package | Category | Affects the estimates/forecasts? | What it does — and what breaks without it |
|---|---|---|---|
| `stochvol` | estimation | **Yes** (SV members only) | Samples the time-varying volatility path for `small_sv`, `medium_minn`, `ucsv`. Without it those members can't be fit (drop them and the rest of the suite still runs). |
| `coda` | estimation-time **diagnostic** | No | Computes MCMC effective sample size for the §9 convergence gate. Every VAR fit calls it, but it never changes a number — purely a sampler quality check. Stub-able. |
| `scoringRules` | evaluation + combination | No (post-estimation) | Scores forecasts (`crps_sample`/`logs_sample`). Drives the scorecard, DM tests, and the data-driven pools. Without it: still get member forecasts + `combo_equal`, but no scores/DM/scorecard/weighted pools. |
| `readrba`,`readabs`,`fredr` | data (real only) | n/a | Fetch the real RBA/ABS/FRED panel. Unused on the synthetic path. |
| `future`,`furrr` | compute | No | Parallelise the OOS loop. Without them it runs serially (slower, same results). |
| `yaml` | plumbing | n/a | Parses `config/config.yml`. |
| `digest` | plumbing | Indirectly | Hashes the config (OOS cache key) and derives per-task seeds (reproducibility). |
| `ggplot2` | output | No | All figures. |
| `glue` | logging | No | Interpolates `{var}` in log messages (fallback path only). |

### 3.1 The three statistical packages: what is *estimation* vs *evaluation*

A common point of confusion is whether all three are "needed for estimation".
They are not — only `stochvol` changes what is actually estimated:

- **`stochvol` — genuine estimation, SV members only.** `specify_priors()` +
  the `sv_*()` prior constructors and the fast C++ sampler
  `svsample_fast_cpp()` estimate the stochastic-volatility path for `small_sv`,
  `medium_minn` and the `ucsv` benchmark (`engines.R`, `benchmarks.R`). The
  `sv_infinity()`/`sv_exponential()` switch is the Gaussian-vs-t-error COVID
  treatment. The constant-volatility engines (conjugate, Gibbs-NIW,
  steady-state) and `rw`/`ar4`/`ucmean` roll their own linear algebra and do
  **not** use it. This is the only one of the three that affects the numbers.

- **`coda` — estimation-*time*, but only a diagnostic.** `effectiveSize()` is
  called inside `mcmc_diagnostics()`, which **every VAR engine runs while
  fitting** (`engines.R` `conj_br`/`gibbs`/`ss`/`sv`, plus `ucsv`). So it sits
  on the estimation path and the engine code errors without it — but all it does
  is compute the posterior effective sample size for the §9 convergence gate
  (`converged = ess_min > 50` → `assert_diagnostics()`). It does **not** touch
  any estimate or forecast; replace it with `NA`/`Inf` and every output number
  is identical, you just lose the convergence assertion. Mechanically required,
  statistically inert, trivially stub-able.

- **`scoringRules` — entirely post-estimation.** It appears only in
  `evaluate.R` and `combine.R` — never in any model-fitting code. `crps_sample()`
  / `logs_sample()` *grade* forecasts that already exist, which feeds (a) the
  OOS scoring, DM tests and scorecard performance tables, and (b) the
  **data-driven combination weights** (`combo_logscore`, `combo_pool`,
  `combo_bma` derive their weights from the log scores). So: evaluation + the
  weighted pools, not estimation.

**Minimum-install intuition:** for *forecasts alone* (no scoring, no convergence
gate, no SV members) you'd need none of the three. Add `stochvol` to restore
the SV members, `coda` to restore the convergence gate, `scoringRules` to
restore scoring + the scorecard + the weighted pools.

### 3.2 The other required packages (already on the RBA list)

**Data acquisition — only for `data.source: real`:**

- **`readrba`** (`read_rba()`), **`readabs`** (`read_abs()`), **`fredr`**
  (`fredr()` / `fredr_set_key()`) — the three providers behind the real panel,
  dispatched per-series in `data_sources.R`: RBA (cash rate, commodities, TWI),
  ABS (GDP, CPI, unemployment), FRED (US foreign block). The synthetic DGP path
  uses none of these, so an offline `data.source: synthetic` run needs none of
  the three.

**Parallel compute:**

- **`future`** + **`furrr`** — `future::plan(multisession)` spawns the worker
  processes and `furrr::future_map(..., seed = TRUE)` fans the per-origin
  forecasts across them with reproducible RNG (`evaluate.R`). `furrr` is just
  the purrr-style API on top of `future`. Dropping them would only force a slow
  serial run; results are unchanged.

**Plumbing & output:**

- **`yaml`** — `read_yaml()` parses `config/config.yml`, the single source of
  truth for variables, suite, MCMC and evaluation settings (`utils.R`).
- **`digest`** — content hashing for two jobs: `config_hash()` keys the on-disk
  OOS cache so stale results aren't served, and `derive_seed()` turns the master
  seed + a task string into a deterministic per-task seed (`utils.R`) — the
  basis of run-to-run reproducibility regardless of worker scheduling.
- **`ggplot2`** — every figure (fan charts, PIT histograms, weight-evolution,
  CRPS-by-horizon) in `report.R`.
- **`glue`** — interpolates `{var}` templates in log messages, used only by the
  logging fallback in `aaa_logging.R`; degrades to the raw template if absent.

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

### `scoringRules` / `stochvol` / `coda` made optional — `R/aaa_capabilities.R`

Capability flags (`has_scoringrules()`, `has_stochvol()`, `has_coda()`) gate the
points that use each package, so a missing one degrades instead of erroring:

- **`stochvol` absent** — `all_members()` drops the SV members (`small_sv`,
  `medium_minn`) and the `ucsv` benchmark (with a logged warning); the rest of
  the suite runs normally.
- **`coda` absent** — `mcmc_diagnostics()` returns `ess_min = NA` and lets the
  §9 convergence gate pass; estimates and forecasts are unchanged.
- **`scoringRules` absent** — scoring (`safe_crps`/`safe_logs`) returns `NA`,
  every combination scheme falls back to **equal weights**, and the scorecard
  prints a one-paragraph note in place of the performance/DM sections. Member
  forecasts and the equal-weighted pool are still produced.

`run_all_plain.R` no longer lists these three as hard requirements: it prints a
NOTE for each one that is missing, explaining what will be reduced. Each can
also be force-disabled for testing, e.g. `options(soe.disable_scoringRules =
TRUE)`.

---

## 5. How to run at the RBA

```sh
# recommended (for the FULL output) — the run also works without these:
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
| `scoringRules` | scoring, DM tests, scorecard performance tables, weighted pools | forecasts + equal-weighted pool still produced; performance sections replaced by a note |
| `stochvol` | the SV members (`small_sv`, `medium_minn`) + `ucsv` | those members are dropped; the rest of the suite runs |
| `coda` | the §9 MCMC convergence diagnostic | `ess_min = NA`, gate auto-passes; estimates unchanged |
| `quarto` CLI | the scorecard `.pdf` | the scorecard `.md` is still written; render the PDF elsewhere |
| `quarto` R package + CLI | the HTML narrative report (`report.html`) | step is skipped with a warning; all other outputs produced |
| `rdbnomics` | the `alt_foreign` variable set only | the default small + medium suite does not use it |

`logger` and `targets` are **no longer required**.
