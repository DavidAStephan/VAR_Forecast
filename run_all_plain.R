#!/usr/bin/env Rscript
# run_all_plain.R -- run the FULL pipeline WITHOUT the `targets` package.
#
# This is a drop-in alternative to `Rscript run_all.R` for environments where
# `targets` is unavailable (e.g. the RBA work setup). It runs exactly the same
# steps as _targets.R, in dependency order, calling the same project functions.
#
# What you DO and DON'T get vs the targets pipeline:
#   * The expensive recursive OOS step is still disk-cached by run_oos_member()
#     (cache/oos_<hash>/), so re-runs skip estimation that has not changed.
#   * You LOSE targets' object-level caching of the cheap downstream steps
#     (scoring, combination, DM, scorecard, figures) -- they recompute each run
#     (seconds), and there is no `tar_read()` / DAG introspection.
#
# Hard package requirements: scoringRules, stochvol, coda (plus the data/compute
# stack you already have: readrba, readabs, fredr, furrr, future, ggplot2, yaml,
# digest, glue). `logger` is optional (facade in R/aaa_logging.R); `targets` is
# not used; `rdbnomics` is only needed for the alt_foreign variable set; the
# `quarto` R package is only needed for the optional HTML report.
#
# Real data needs FRED_API_KEY in a gitignored .Renviron; for a fully offline
# run set data.source: synthetic in config/config.yml.

if (file.exists(".Renviron")) readRenviron(".Renviron")

if (file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)) {
  message("Restoring renv library (first run may take a while)...")
  options(renv.config.install.verbose = FALSE)
  try(renv::restore(prompt = FALSE), silent = TRUE)
}

# fail early if a TRULY required package is missing (the pipeline cannot run)
.need <- c("yaml", "digest", "furrr", "future", "ggplot2")
.missing <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing))
  stop("missing required packages: ", paste(.missing, collapse = ", "),
       "\nInstall them (e.g. via your Artifactory CRAN mirror) and re-run.")

# the three statistical packages are OPTIONAL: the run degrades instead of
# failing (see R/aaa_capabilities.R). Report what will be reduced.
.optional <- c(
  scoringRules = "scoring + DM tests + scorecard performance tables (forecasts still pooled with EQUAL weights)",
  stochvol     = "the SV members (small_sv, medium_minn) and the ucsv benchmark are skipped",
  coda         = "the MCMC convergence diagnostic (estimates unchanged)")
for (p in names(.optional))
  if (!requireNamespace(p, quietly = TRUE))
    message(sprintf("NOTE: '%s' not installed -- degraded: %s", p, .optional[[p]]))

# project functions; aaa_logging.R sorts first so the logging facade exists
# before anything calls it (also re-sourced in workers by ensure_project_loaded)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

main <- function() {
  setup_logging()

  # ---- config, spec, data ----
  cfg  <- load_config("config/config.yml")
  spec <- build_transform_spec(cfg)
  dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
  write.csv(spec, "output/tables/transform_spec.csv", row.names = FALSE)

  raw_data    <- get_raw_data(cfg, spec)
  td          <- transform_data(raw_data, spec)
  data_checks <- check_data(td, spec)

  members <- all_members(cfg)

  # ---- recursive OOS (parallel over origins; disk-cached) ----
  future::plan(future::multisession, workers = cfg$parallel$workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  oos_all <- setNames(
    lapply(members, function(m)
      timed(paste0("OOS ", m$name), run_oos_member(m, td, spec, cfg))),
    vapply(members, `[[`, "", "name"))
  future::plan(future::sequential)

  # ---- scoring + density combination ----
  scores <- do.call(rbind, lapply(names(oos_all), function(m)
    score_member(m, oos_all[[m]], td, spec, cfg)))
  draws_env <- collect_draws(oos_all, cfg)
  combos    <- timed("combination", combine_all(scores, draws_env, td, spec, cfg))
  allscores <- rbind(scores, combos$scores)

  # ---- Diebold-Mariano tests ----
  dm_table <- rbind(dm_vs_reference(allscores, "ar4"),
                    dm_vs_reference(allscores, "rw"))

  # ---- section-9 self-checks + diagnostics gate ----
  check_nla   <- test_no_lookahead(td, spec, cfg)
  check_repro <- test_reproducibility(td, spec, cfg)
  calibration <- check_calibration(allscores)
  diag_tab <- diagnostics_table(oos_all,
    checks = list(no_lookahead = check_nla, reproducible = check_repro))
  assert_diagnostics(diag_tab)

  # ---- final full-sample forecasts + outputs ----
  set.seed(derive_seed(cfg$master_seed, "final"))
  final_fc   <- timed("final forecasts", final_forecasts(td, spec, cfg, scores))
  fc_table    <- forecast_tables(final_fc, td, spec, cfg)
  eval_tables <- evaluation_tables(allscores, dm_table, cfg)

  qstr <- function(d) paste0(format(d, "%Y"), "Q",
                             (as.integer(format(d, "%m")) - 1) %/% 3 + 1)
  panel <- paste0(qstr(min(td$date)), "-", qstr(max(td$date)))
  write_model_scorecard(allscores, spec, dm_table, diag_tab, cfg,
                        src = cfg$data$source, panel = panel,
                        out_path = "reports/model_scorecard.md")

  write.csv(final_fc$weights, "output/tables/final_weights.csv", row.names = FALSE)
  write.csv(combos$weights, "output/tables/weights_over_time.csv", row.names = FALSE)
  make_figures(allscores, combos$weights, final_fc, td, spec, cfg)

  # ---- optional Quarto HTML report (skipped if quarto unavailable) ----
  if (requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))) {
    tryCatch(quarto::quarto_render("reports/report.qmd", quiet = TRUE),
             error = function(e) log_warn("report render failed: {conditionMessage(e)}"))
  } else {
    log_warn("quarto R package / CLI unavailable; skipping HTML report")
  }

  log_info("pipeline complete")
  invisible(TRUE)
}

main()

cat("\nDone. Outputs:\n",
    "  reports/model_scorecard.md  (render the PDF with reports/render_scorecard_pdf.sh)\n",
    "  output/tables/     evaluation + diagnostics tables\n",
    "  output/figures/    fan charts, PIT, weights, score plots\n",
    "  output/forecasts/  forecast tables\n",
    "  reports/report.html (only if the quarto R package + CLI are available)\n")
