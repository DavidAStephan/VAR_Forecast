# aaa_capabilities.R -- optional-package capability flags -----------------------
#
# The three statistical packages scoringRules, stochvol and coda are OPTIONAL:
# the suite degrades gracefully when one is absent instead of failing to load or
# crashing. What is lost in each case:
#   * stochvol missing -> the SV members (small_sv, medium_minn) and the ucsv
#     benchmark are dropped from the suite (all_members()); everything else runs.
#   * coda missing -> the MCMC effective-sample-size diagnostic is skipped
#     (ess_min = NA, the convergence gate auto-passes); estimates are unchanged.
#   * scoringRules missing -> forecasts are still produced and pooled with EQUAL
#     weights, but CRPS / log-score / RMSE scoring, the DM tests and the
#     scorecard performance tables are skipped (crps/logdens recorded as NA).
#
# Each capability can also be force-disabled (e.g. to test a degraded path or
# match a target environment) with an option:
#     options(soe.disable_scoringRules = TRUE)   # also stochvol / coda
# This file sorts first in R/ so the flags exist before any caller (including
# parallel workers via ensure_project_loaded()). It is intentionally NOT in
# config_hash()'s estimation-file list: it changes which members run and which
# diagnostics are computed, but never the numerical estimate of a member that
# does run, so it must not invalidate the OOS cache.

#' TRUE if `pkg` is installed AND not force-disabled via options(soe.disable_<pkg>).
has_pkg <- function(pkg)
  !isTRUE(getOption(paste0("soe.disable_", pkg), FALSE)) &&
    requireNamespace(pkg, quietly = TRUE)

has_scoringrules <- function() has_pkg("scoringRules")
has_stochvol     <- function() has_pkg("stochvol")
has_coda         <- function() has_pkg("coda")

# ---- safe wrappers (used in place of the bare pkg:: calls) -------------------

#' MCMC effective sample size; NA when coda is unavailable (diagnostic only,
#' so a NA never changes an estimate or forecast).
safe_ess <- function(x)
  if (has_coda()) as.numeric(coda::effectiveSize(x)) else NA_real_

#' CRPS of a sample forecast; NA_real_ when scoringRules is unavailable.
safe_crps <- function(y, dat)
  if (has_scoringrules()) scoringRules::crps_sample(y, dat) else NA_real_

#' Log predictive score of a sample forecast; NA_real_ when scoringRules is
#' unavailable. (score_member negates this to store logdens.)
safe_logs <- function(y, dat)
  if (has_scoringrules()) scoringRules::logs_sample(y, dat) else NA_real_
