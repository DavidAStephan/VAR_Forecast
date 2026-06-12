# report.R -- final forecasts, tables and figures -------------------------------

suppressPackageStartupMessages({ library(ggplot2) })

#' Fit every member on the FULL sample and produce h=1..12 predictive draws,
#' combined with weights estimated from the complete OOS score history.
final_forecasts <- function(td, spec, cfg, scores) {
  H <- cfg$horizons
  members <- all_members(cfg)
  T_n <- nrow(td)
  draws <- list()
  for (m in members) {
    set.seed(derive_seed(cfg$master_seed, paste0("final-", m$name)))
    fc <- forecast_at_origin(m, td, spec, cfg)
    draws[[m$name]] <- fc$draws
  }
  tgt <- spec$variable[spec$target]
  buckets <- cfg$combination$horizon_buckets
  schemes <- unlist(cfg$combination$schemes)
  mnames <- names(draws)
  pooled <- list(); wtab <- list()
  for (scheme in schemes) for (v in tgt) {
    arr <- array(NA_real_, c(dim(draws[[1]])[1], H, length(mnames)),
                 dimnames = list(NULL, NULL, mnames))
    for (m in mnames) arr[, , m] <- draws[[m]][, , v]
    pool_dr <- matrix(NA_real_, 600, H)
    for (bn in names(buckets)) {
      hs <- unlist(buckets[[bn]])
      w <- combo_weights(scheme, scores, v, hs, T_n + 1, mnames, cfg)
      wtab[[length(wtab) + 1]] <- data.frame(
        scheme = scheme, variable = v, bucket = bn,
        member = names(w), weight = as.numeric(w))
      set.seed(derive_seed(cfg$master_seed, paste("finalmix", scheme, v, bn)))
      cnt <- drop(rmultinom(1, 600, w))
      for (h in hs) {
        mix <- unlist(lapply(seq_along(mnames), function(i) {
          if (cnt[i] == 0) return(numeric(0))
          sample(arr[, h, mnames[i]], cnt[i], replace = TRUE)
        }))
        pool_dr[, h] <- mix
      }
    }
    pooled[[paste(scheme, v, sep = "|")]] <- pool_dr
  }
  list(member_draws = draws, pooled = pooled, weights = do.call(rbind, wtab),
       origin = T_n, origin_date = td$date[T_n])
}

#' Forecast table: point + interval for combined and per-member forecasts.
forecast_tables <- function(ff, td, spec, cfg, out_dir = "output/forecasts") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tgt <- spec$variable[spec$target]
  qs <- c(0.05, 0.5, 0.95)
  rows <- list()
  for (key in names(ff$pooled)) {
    sv <- strsplit(key, "|", fixed = TRUE)[[1]]
    qm <- t(apply(ff$pooled[[key]], 2, quantile, probs = qs))
    rows[[length(rows) + 1]] <- data.frame(
      member = paste0("combo_", sv[1]), variable = sv[2], h = seq_len(nrow(qm)),
      point = qm[, 2], lo90 = qm[, 1], hi90 = qm[, 3])
  }
  for (m in names(ff$member_draws)) for (v in tgt) {
    dr <- ff$member_draws[[m]][, , v]
    qm <- t(apply(dr, 2, quantile, probs = qs))
    rows[[length(rows) + 1]] <- data.frame(
      member = m, variable = v, h = seq_len(nrow(qm)),
      point = qm[, 2], lo90 = qm[, 1], hi90 = qm[, 3])
  }
  tab <- do.call(rbind, rows)
  tab$origin_date <- ff$origin_date
  write.csv(tab, file.path(out_dir, "forecast_table.csv"), row.names = FALSE)
  tab
}

#' Evaluation tables: scores by horizon with DM significance markers.
evaluation_tables <- function(allscores, dm, out_dir = "output/tables") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sm <- summarise_scores(allscores)
  write.csv(sm, file.path(out_dir, "scores_by_horizon.csv"), row.names = FALSE)
  if (!is.null(dm) && nrow(dm)) {
    dm$sig_crps <- cut(dm$p_crps, c(0, 0.01, 0.05, 0.1, 1),
                       labels = c("***", "**", "*", ""), include.lowest = TRUE)
    write.csv(dm, file.path(out_dir, "dm_tests.csv"), row.names = FALSE)
  }
  sm
}

# ---- figures ----------------------------------------------------------------------

.theme <- function() theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

#' Fan chart for one variable from pooled draws.
plot_fan <- function(pool_dr, td, v, spec, cfg, scheme,
                     out_dir = "output/figures", hist_quarters = 24) {
  H <- ncol(pool_dr)
  probs <- unlist(cfg$report$fan_quantiles)
  qm <- apply(pool_dr, 2, quantile, probs = probs)
  fdates <- seq(td$date[nrow(td)], by = "quarter", length.out = H + 1)[-1]
  hist <- tail(data.frame(date = td$date, value = td[[v]]), hist_quarters)
  bands <- data.frame()
  np <- length(probs)
  for (i in seq_len(floor(np / 2))) {
    bands <- rbind(bands, data.frame(
      date = fdates, lo = qm[i, ], hi = qm[np + 1 - i, ],
      band = factor(i)))
  }
  med <- data.frame(date = fdates, value = qm[(np + 1) / 2, ])
  lbl <- spec$label[spec$variable == v]
  g <- ggplot() +
    geom_ribbon(data = bands, aes(date, ymin = lo, ymax = hi, group = band),
                fill = "steelblue", alpha = 0.18) +
    geom_line(data = hist, aes(date, value), linewidth = 0.5) +
    geom_line(data = med, aes(date, value), color = "steelblue4",
              linewidth = 0.7, linetype = "21") +
    labs(title = paste0(lbl, " — ", scheme, " pool"), x = NULL, y = NULL) +
    .theme()
  f <- file.path(out_dir, sprintf("fan_%s_%s.png", scheme, v))
  ggsave(f, g, width = 7, height = 4, dpi = 130)
  f
}

#' PIT histograms by member at a given horizon.
plot_pits <- function(allscores, h = 1, out_dir = "output/figures") {
  d <- allscores[allscores$measure == "q" & allscores$h == h, ]
  g <- ggplot(d, aes(pit)) +
    geom_histogram(breaks = seq(0, 1, 0.1), fill = "steelblue", color = "white") +
    facet_wrap(~member, scales = "free_y") +
    labs(title = sprintf("PIT histograms, h = %d (uniform = calibrated)", h),
         x = "PIT", y = NULL) + .theme()
  f <- file.path(out_dir, sprintf("pit_h%d.png", h))
  ggsave(f, g, width = 9, height = 6, dpi = 130)
  f
}

#' Combination weights over origins (stacked area), per variable/bucket.
plot_weights <- function(weights, td, scheme = "pool", v, bucket = "near",
                         out_dir = "output/figures") {
  d <- weights[weights$scheme == scheme & weights$variable == v &
               weights$bucket == bucket, ]
  if (!nrow(d)) return(NULL)
  d$date <- td$date[d$origin]
  g <- ggplot(d, aes(date, weight, fill = member)) +
    geom_area(position = "stack") +
    labs(title = sprintf("%s weights over time — %s (%s horizons)", scheme, v, bucket),
         x = NULL, y = "weight") + .theme()
  f <- file.path(out_dir, sprintf("weights_%s_%s_%s.png", scheme, v, bucket))
  ggsave(f, g, width = 8, height = 4.5, dpi = 130)
  f
}

#' Score-by-horizon plot (mean CRPS by h), members vs combinations.
plot_scores_by_h <- function(allscores, v, out_dir = "output/figures") {
  sm <- summarise_scores(allscores[allscores$variable == v &
                                   allscores$measure == "q", ])
  sm$type <- ifelse(grepl("^combo_", sm$member), "combination", "member")
  g <- ggplot(sm, aes(h, crps, color = member, linetype = type)) +
    geom_line() + geom_point(size = 0.8) +
    labs(title = paste0("Mean CRPS by horizon — ", v), x = "horizon (quarters)",
         y = "CRPS (lower = better)") + .theme()
  f <- file.path(out_dir, sprintf("crps_by_h_%s.png", v))
  ggsave(f, g, width = 8, height = 5, dpi = 130)
  f
}

#' Produce every figure; returns the file list.
make_figures <- function(allscores, cmb_weights, ff, td, spec, cfg,
                         out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tgt <- spec$variable[spec$target]
  files <- c()
  for (v in tgt) {
    key <- paste("pool", v, sep = "|")
    if (key %in% names(ff$pooled))
      files <- c(files, plot_fan(ff$pooled[[key]], td, v, spec, cfg, "pool"))
    files <- c(files, plot_scores_by_h(allscores, v))
    for (b in names(cfg$combination$horizon_buckets))
      files <- c(files, plot_weights(cmb_weights, td, "pool", v, b))
  }
  files <- c(files, plot_pits(allscores, 1), plot_pits(allscores, 4))
  files[!vapply(files, is.null, logical(1))]
}

#' Section 9 diagnostics summary table, written to output/tables.
diagnostics_table <- function(oos_all, checks, out_dir = "output/tables") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(names(oos_all), function(m) {
    dg <- lapply(oos_all[[m]], function(r) r$diagnostics)
    sn <- lapply(oos_all[[m]], function(r) r$sanity)
    data.frame(
      member = m,
      ess_min = min(vapply(dg, function(d) d$ess_min, numeric(1)), na.rm = TRUE),
      stable_share_min = suppressWarnings(min(vapply(dg, function(d)
        ifelse(is.null(d$stable_share) || is.na(d$stable_share), 1,
               d$stable_share), numeric(1)))),
      block_exog_max = max(vapply(dg, function(d)
        ifelse(is.null(d$block_exog_max), 0, d$block_exog_max), numeric(1))),
      converged_all = all(vapply(dg, function(d) isTRUE(d$converged), logical(1))),
      sanity_all = all(vapply(sn, function(s) isTRUE(s$ok), logical(1))))
  })
  tab <- do.call(rbind, rows)
  tab$no_lookahead <- checks$no_lookahead
  tab$reproducible <- checks$reproducible
  write.csv(tab, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)
  tab
}

#' Hard gate: stop the pipeline if a section-9 correctness property fails.
assert_diagnostics <- function(diag_tab) {
  stopifnot(
    "block exogeneity violated"   = all(diag_tab$block_exog_max < 1e-2),
    "MCMC convergence failure"    = all(diag_tab$converged_all),
    "explosive forecasts"         = all(diag_tab$sanity_all),
    "look-ahead detected"         = all(diag_tab$no_lookahead),
    "not reproducible"            = all(diag_tab$reproducible))
  logger::log_info("all section-9 diagnostics passed")
  invisible(TRUE)
}
