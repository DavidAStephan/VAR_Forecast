# _targets.R -- pipeline definition ---------------------------------------------
library(targets)

tar_option_set(
  packages = c("yaml", "logger", "digest", "stochvol", "coda", "scoringRules",
               "furrr", "future", "ggplot2"),
  format = "rds"
)

# project functions (also re-sourced inside parallel workers)
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

setup_logging()

list(
  tar_target(cfg_file, "config/config.yml", format = "file"),
  tar_target(cfg, load_config(cfg_file)),
  tar_target(spec, build_transform_spec(cfg)),
  tar_target(spec_csv, {
    dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
    write.csv(spec, "output/tables/transform_spec.csv", row.names = FALSE)
    "output/tables/transform_spec.csv"
  }, format = "file"),

  tar_target(raw_data, get_raw_data(cfg, spec)),
  tar_target(td, transform_data(raw_data, spec)),
  tar_target(data_checks, check_data(td, spec)),

  tar_target(members, all_members(cfg)),

  # recursive OOS, one branch per member; parallel over origins inside
  tar_target(oos_member, {
    data_checks
    future::plan(future::multisession, workers = cfg$parallel$workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    timed(paste0("OOS ", members[[1]]$name),
          run_oos_member(members[[1]], td, spec, cfg))
  }, pattern = map(members), iteration = "list"),

  tar_target(oos_all, {
    x <- oos_member
    names(x) <- vapply(members, `[[`, "", "name")
    x
  }),

  tar_target(scores, {
    do.call(rbind, lapply(names(oos_all), function(m)
      score_member(m, oos_all[[m]], td, spec, cfg)))
  }),

  tar_target(combos, {
    draws_env <- collect_draws(oos_all, cfg)
    timed("combination", combine_all(scores, draws_env, td, spec, cfg))
  }),

  tar_target(allscores, rbind(scores, combos$scores)),

  tar_target(dm_table, {
    rbind(dm_vs_reference(allscores, "ar4"),
          dm_vs_reference(allscores, "rw"))
  }),

  # section 9 correctness checks
  tar_target(check_nla, test_no_lookahead(td, spec, cfg)),
  tar_target(check_repro, test_reproducibility(td, spec, cfg)),
  tar_target(calibration, check_calibration(allscores)),

  tar_target(diag_tab, {
    tab <- diagnostics_table(oos_all,
      checks = list(no_lookahead = check_nla, reproducible = check_repro))
    assert_diagnostics(tab)
    tab
  }),

  # final full-sample forecasts and outputs
  tar_target(final_fc, {
    set.seed(derive_seed(cfg$master_seed, "final"))
    timed("final forecasts", final_forecasts(td, spec, cfg, scores))
  }),
  tar_target(fc_table, forecast_tables(final_fc, td, spec, cfg)),
  tar_target(eval_tables, evaluation_tables(allscores, dm_table, cfg)),
  tar_target(weight_table, {
    dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
    write.csv(final_fc$weights, "output/tables/final_weights.csv", row.names = FALSE)
    write.csv(combos$weights, "output/tables/weights_over_time.csv", row.names = FALSE)
    final_fc$weights
  }),
  tar_target(figures, make_figures(allscores, combos$weights, final_fc, td, spec, cfg)),

  tar_target(report, {
    figures; fc_table; eval_tables; diag_tab; calibration; weight_table
    if (nzchar(Sys.which("quarto"))) {
      tryCatch({
        quarto::quarto_render("reports/report.qmd", quiet = TRUE)
        "reports/report.html"
      }, error = function(e) {
        logger::log_warn("report render failed: {conditionMessage(e)}")
        "render-failed"
      })
    } else {
      logger::log_warn("quarto CLI not found; skipping report render")
      "quarto-missing"
    }
  })
)
