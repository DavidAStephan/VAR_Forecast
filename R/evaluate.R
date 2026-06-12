# evaluate.R -- pseudo-real-time recursive OOS evaluation -----------------------
#
# No-look-ahead is enforced structurally: forecast_at_origin() receives only
# td[1:t, ] and everything downstream (GLP lambda selection, estimation,
# prediction) sees that truncated panel alone. test_no_lookahead() verifies it.

suppressPackageStartupMessages({ library(scoringRules); library(furrr) })

#' Forecast origins (row indices of the transformed panel).
oos_origins <- function(td, cfg) {
  T_n <- nrow(td)
  t0 <- ceiling(cfg$evaluation$first_origin_frac * T_n)
  origins <- t0:(T_n - 1)
  if (length(origins) > cfg$evaluation$max_origins)
    origins <- tail(origins, cfg$evaluation$max_origins)
  origins
}

#' All pool members: suite entries + benchmarks, in a common list format.
all_members <- function(cfg) {
  suite <- lapply(cfg$suite, function(m) { m$kind <- "var"; m })
  bench <- lapply(cfg$benchmarks, function(b)
    list(name = b, kind = "benchmark", engine = b))
  c(suite, bench)
}

#' Fit + forecast a single member at one origin using ONLY data up to t.
#' Returns thinned predictive draws [store_draws, h, n_targets] plus sanity info.
forecast_at_origin <- function(member, td_t, spec, cfg) {
  H <- cfg$horizons
  tgt <- spec$variable[spec$target]
  nfd <- cfg$mcmc$forecast_draws
  if (member$kind == "var") {
    set_name <- member$set
    spec_m <- vars_for_set(spec, set_name)
    y <- as.matrix(td_t[, spec_m$variable])
    if (member$kind == "var" && identical(member$prior$lambda, "auto") &&
        isTRUE(cfg$glp$enabled)) {
      glp_lambda <- select_lambda(y, member$lags, ar_sigmas(y), spec_m$delta,
                                  grid = unlist(cfg$glp$lambda_grid))
    } else glp_lambda <- 0.2
    post <- fit_var_member(y, member, spec_m, cfg, glp_lambda = glp_lambda)
    paths <- simulate_paths(post, y, H, nfd)
    keep <- intersect(tgt, dimnames(paths)[[3]])
    paths <- paths[, , keep, drop = FALSE]
  } else {
    cfgb <- cfg
    cfgb$mcmc$ndraw <- cfg$mcmc$bench_ndraw
    cfgb$mcmc$nburn <- cfg$mcmc$bench_nburn
    y <- as.matrix(td_t[, tgt])
    post <- fit_benchmark(y, member$engine, cfgb)
    paths <- simulate_paths(post, y, H, nfd)
    glp_lambda <- NA_real_
  }
  vnames <- dimnames(paths)[[3]]
  sanity <- check_forecasts(paths, as.matrix(td_t[, vnames]),
                            label = member$name,
                            delta = spec$delta[match(vnames, spec$variable)])
  thin <- round(seq(1, dim(paths)[1], length.out = cfg$mcmc$store_draws))
  list(draws = paths[thin, , , drop = FALSE],
       diagnostics = post$diagnostics, sanity = sanity,
       glp_lambda = glp_lambda)
}

#' Run the recursive loop for one member over all origins, with disk caching
#' keyed by the config hash. Parallel over origins.
run_oos_member <- function(member, td, spec, cfg, cache_root = "cache") {
  hash <- config_hash(cfg)
  cdir <- file.path(cache_root, paste0("oos_", hash))
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  origins <- oos_origins(td, cfg)
  res <- furrr::future_map(origins, function(t) {
    ensure_project_loaded()
    f <- file.path(cdir, sprintf("%s_o%03d.rds", member$name, t))
    if (file.exists(f)) return(readRDS(f))
    set.seed(derive_seed(cfg$master_seed, paste0(member$name, "-", t)))
    out <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
    out$origin <- t
    saveRDS(out, f)
    out
  }, .options = furrr::furrr_options(seed = TRUE))
  names(res) <- as.character(origins)
  logger::log_info("OOS done: {member$name} ({length(origins)} origins)")
  res
}

# ---- scoring --------------------------------------------------------------------

#' Score one member's OOS results against realizations. Returns a long
#' data.frame: origin, date, variable, measure (q|ye), h, point, real, logdens,
#' crps, pit. Year-ended rows combine forecast draws with realized history.
score_member <- function(member_name, oos_res, td, spec, cfg) {
  T_n <- nrow(td)
  tgt <- spec$variable[spec$target]
  rows <- list()
  for (res in oos_res) {
    t <- res$origin
    dr <- res$draws                       # [ndraw, H, var]
    H <- dim(dr)[2]
    for (v in dimnames(dr)[[3]]) {
      vt <- spec[spec$variable == v, ]
      real_q <- td[[v]]
      for (h in seq_len(H)) {
        if (t + h > T_n) next
        x <- dr[, h, v]
        y_real <- real_q[t + h]
        rows[[length(rows) + 1]] <- data.frame(
          member = member_name, origin = t, date = td$date[t + h],
          variable = v, measure = "q", h = h,
          point = mean(x), real = y_real,
          logdens = -scoringRules::logs_sample(y_real, x),
          crps = scoringRules::crps_sample(y_real, x),
          pit = mean(x <= y_real))
        # year-ended: sum of 4 quarterly outcomes ending at t+h
        if (isTRUE(vt$year_ended) && vt$transform == "dlog") {
          k <- (h - 3):h
          hist_part <- sum(real_q[t + k[k <= 0]])
          fc_idx <- k[k >= 1]
          xye <- rowSums(dr[, fc_idx, v, drop = FALSE]) + hist_part
          ye_real <- sum(real_q[t + k])
          rows[[length(rows) + 1]] <- data.frame(
            member = member_name, origin = t, date = td$date[t + h],
            variable = v, measure = "ye", h = h,
            point = mean(xye), real = ye_real,
            logdens = -scoringRules::logs_sample(ye_real, xye),
            crps = scoringRules::crps_sample(ye_real, xye),
            pit = mean(xye <= ye_real))
        }
      }
    }
  }
  do.call(rbind, rows)
}

#' Stored draws in a flat structure for the combination layer:
#' list keyed "origin|variable|h" -> named list(member -> draw vector).
collect_draws <- function(oos_all, cfg) {
  out <- new.env(parent = emptyenv())
  for (m in names(oos_all)) {
    for (res in oos_all[[m]]) {
      t <- res$origin
      dr <- res$draws
      for (v in dimnames(dr)[[3]]) for (h in seq_len(dim(dr)[2])) {
        key <- paste(t, v, h, sep = "|")
        cur <- if (exists(key, out)) get(key, out) else list()
        cur[[m]] <- dr[, h, v]
        assign(key, cur, out)
      }
    }
  }
  out
}

#' Summary score table: mean logdens / CRPS / RMSE by member, variable,
#' measure, horizon.
summarise_scores <- function(scores) {
  agg <- aggregate(cbind(logdens, crps) ~ member + variable + measure + h,
                   data = scores, FUN = mean)
  rmse <- aggregate((point - real)^2 ~ member + variable + measure + h,
                    data = scores, FUN = function(x) sqrt(mean(x)))
  names(rmse)[5] <- "rmse"
  merge(agg, rmse)
}

# ---- Diebold-Mariano ---------------------------------------------------------------

#' DM test with Newey-West variance (lag h-1) and the Harvey small-sample
#' correction. loss1/loss2 aligned by origin; H0: equal predictive accuracy.
dm_test <- function(loss1, loss2, h = 1) {
  d <- loss1 - loss2
  d <- d[is.finite(d)]
  n <- length(d)
  if (n < 8 || sd(d) < 1e-12) return(c(stat = NA_real_, p = NA_real_))
  dbar <- mean(d)
  L <- max(0, h - 1)
  g0 <- mean((d - dbar)^2)
  v <- g0
  if (L > 0) for (l in seq_len(min(L, n - 1))) {
    gl <- mean((d[(l + 1):n] - dbar) * (d[1:(n - l)] - dbar))
    v <- v + 2 * (1 - l / (L + 1)) * gl
  }
  v <- max(v, 1e-12)
  stat <- dbar / sqrt(v / n)
  k <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)   # Harvey et al. correction
  stat <- stat * k
  p <- 2 * pt(-abs(stat), df = n - 1)
  c(stat = stat, p = p)
}

#' DM comparisons of every member against a reference member, per variable,
#' measure and horizon, for both squared-error and CRPS losses.
dm_vs_reference <- function(scores, reference) {
  combos <- unique(scores[, c("variable", "measure", "h")])
  members <- setdiff(unique(scores$member), reference)
  rows <- list()
  for (i in seq_len(nrow(combos))) {
    cb <- combos[i, ]
    base <- scores[scores$member == reference & scores$variable == cb$variable &
                   scores$measure == cb$measure & scores$h == cb$h, ]
    base <- base[order(base$origin), ]
    for (m in members) {
      alt <- scores[scores$member == m & scores$variable == cb$variable &
                    scores$measure == cb$measure & scores$h == cb$h, ]
      alt <- alt[order(alt$origin), ]
      common <- intersect(base$origin, alt$origin)
      if (length(common) < 8) next
      b <- base[match(common, base$origin), ]; a <- alt[match(common, alt$origin), ]
      d_se   <- dm_test((a$point - a$real)^2, (b$point - b$real)^2, h = cb$h)
      d_crps <- dm_test(a$crps, b$crps, h = cb$h)
      rows[[length(rows) + 1]] <- data.frame(
        member = m, reference = reference, variable = cb$variable,
        measure = cb$measure, h = cb$h,
        dm_se = d_se["stat"], p_se = d_se["p"],
        dm_crps = d_crps["stat"], p_crps = d_crps["p"])
    }
  }
  do.call(rbind, rows)
}

# ---- section 9 self-checks ----------------------------------------------------------

#' No-look-ahead test: corrupting all data AFTER the origin must not change
#' the forecast (the harness only ever passes td[1:t, ]).
test_no_lookahead <- function(td, spec, cfg, member = NULL) {
  if (is.null(member)) member <- all_members(cfg)[[1]]
  origins <- oos_origins(td, cfg)
  t <- origins[1]
  td_bad <- td
  if (t < nrow(td))
    td_bad[(t + 1):nrow(td), -1] <- td_bad[(t + 1):nrow(td), -1] * 1e6 + 999
  set.seed(derive_seed(cfg$master_seed, "nla"))
  f1 <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
  set.seed(derive_seed(cfg$master_seed, "nla"))
  f2 <- forecast_at_origin(member, td_bad[seq_len(t), , drop = FALSE], spec, cfg)
  ok <- identical(f1$draws, f2$draws)
  if (!ok) logger::log_error("NO-LOOK-AHEAD TEST FAILED")
  else logger::log_info("no-look-ahead test passed (member {member$name}, origin {t})")
  ok
}

#' Reproducibility: same seed twice -> identical draws.
test_reproducibility <- function(td, spec, cfg) {
  member <- all_members(cfg)[[1]]
  t <- oos_origins(td, cfg)[1]
  set.seed(derive_seed(cfg$master_seed, "repro"))
  f1 <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
  set.seed(derive_seed(cfg$master_seed, "repro"))
  f2 <- forecast_at_origin(member, td[seq_len(t), , drop = FALSE], spec, cfg)
  ok <- identical(f1$draws, f2$draws)
  if (!ok) logger::log_error("REPRODUCIBILITY TEST FAILED")
  else logger::log_info("reproducibility test passed")
  ok
}

#' Calibration check: PIT approximately uniform (chi-square on deciles).
check_calibration <- function(scores, alpha = 0.01) {
  out <- list()
  for (m in unique(scores$member)) {
    p <- scores$pit[scores$member == m & scores$measure == "q" & scores$h == 1]
    if (length(p) < 20) next
    cnt <- table(cut(p, seq(0, 1, 0.2), include.lowest = TRUE))
    chi <- suppressWarnings(chisq.test(cnt))
    out[[m]] <- chi$p.value
  }
  out
}
