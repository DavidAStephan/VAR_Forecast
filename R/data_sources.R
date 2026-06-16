# data_sources.R -- synthetic generator and real downloads with caching --------
#
# Contract: both paths return a list(data = data.frame(date, <one column per
# variable in transform_spec, in RAW LEVEL units>), manifest = data.frame).
# transforms.R turns raw levels into model units. The synthetic path also
# attaches the true DGP (attr "dgp") for the block-exogeneity ground-truth check.

# ---- synthetic ----------------------------------------------------------------

#' Simulate a block-exogenous SOE DGP in transformed units, then integrate back
#' to raw levels so the transform layer is exercised identically to real data.
generate_synthetic_data <- function(cfg, spec) {
  seed <- derive_seed(cfg$master_seed, "synthetic-data")
  set.seed(seed)
  n   <- cfg$synthetic$n_quarters
  p   <- cfg$synthetic$dgp_lags
  M   <- nrow(spec)
  nf  <- sum(spec$block == "foreign")
  stopifnot(all(spec$block[seq_len(nf)] == "foreign"))

  # Steady states in transformed units; fill the NA anchors (100*log indices ~ 460).
  mu <- ifelse(is.na(spec$ss_mean), 460, spec$ss_mean)

  # Own-lag persistence by variable type: levels persistent, growth mildly so.
  rho <- ifelse(spec$delta == 1, 0.92, 0.30)

  # Lag matrices A1, A2 (M x M) with block exogeneity: foreign rows have zero
  # coefficients on domestic columns. Off-diagonals are a fixed, economically
  # structured pattern (random draws compound badly under high persistence).
  A1 <- matrix(0, M, M, dimnames = list(spec$variable, spec$variable))
  A2 <- matrix(0, M, M)
  diag(A1) <- rho * 1.3
  diag(A2) <- -rho * 0.3
  put <- function(eq, on, coef) {
    if (eq %in% rownames(A1) && on %in% rownames(A1)) A1[eq, on] <<- coef
  }
  # foreign block internal dynamics
  put("f_act",    "f_comm", 0.02)
  put("f_rate",   "f_act",  0.06)
  # domestic real activity loads on the foreign block (price-taker channel)
  for (g in c("gdp_growth", "emp_growth", "cons_growth", "wpi_growth")) {
    put(g, "f_act", 0.15); put(g, "f_comm", 0.04)
  }
  put("cpi_inflation", "f_comm",        0.015)
  put("cpi_inflation", "gdp_growth",    0.04)
  put("unemp_rate",    "gdp_growth",   -0.06)
  put("cash_rate",     "cpi_inflation", 0.10)   # Taylor-rule-ish
  put("cash_rate",     "gdp_growth",    0.03)
  put("cash_rate",     "f_rate",        0.05)
  put("bond10y",       "cash_rate",     0.10)
  put("bond10y",       "f_rate",        0.08)
  put("rtwi",          "f_comm",        0.06)
  put("tot",           "f_comm",        0.10)
  put("gdp_growth",    "cash_rate",    -0.05)   # monetary transmission
  put("cpi_inflation", "rtwi",         -0.004)
  stopifnot(all(A1[seq_len(nf), (nf + 1):M] == 0),
            all(A2[seq_len(nf), (nf + 1):M] == 0))

  # shrink until stable
  scale_stable <- function(A1, A2) {
    for (k in 1:50) {
      comp <- rbind(cbind(A1, A2), cbind(diag(M), matrix(0, M, M)))
      if (max(Mod(eigen(comp, only.values = TRUE)$values)) < 0.985) break
      A1 <- A1 * 0.95; A2 <- A2 * 0.95
    }
    list(A1 = A1, A2 = A2)
  }
  As <- scale_stable(A1, A2); A1 <- As$A1; A2 <- As$A2

  # Structural errors: domestic loads contemporaneously on foreign shocks,
  # never the reverse (price-taker).
  sd_e <- ifelse(spec$delta == 1, 0.22, 0.55)
  sd_e[spec$variable == "f_comm"] <- 3.0       # commodity prices are wild
  sd_e[spec$variable == "cpi_inflation"] <- 0.25
  sd_e[spec$variable %in% c("cash_rate", "f_rate")] <- 0.18
  G <- matrix(0, M, M); diag(G) <- 1
  for (i in (nf + 1):M) for (j in seq_len(nf)) G[i, j] <- runif(1, 0.0, 0.25)

  cvec <- drop((diag(M) - A1 - A2) %*% mu)     # intercept consistent with mu

  y <- matrix(0, n + 20, M)                     # 20 burn-in quarters
  y[1:2, ] <- matrix(mu, 2, M, byrow = TRUE)
  for (t in 3:(n + 20)) {
    e <- rnorm(M) * sd_e
    u <- drop(G %*% e)
    y[t, ] <- cvec + drop(A1 %*% y[t - 1, ]) + drop(A2 %*% y[t - 2, ]) + u
  }
  y <- y[-(1:20), , drop = FALSE]

  # COVID-style outlier: hit growth/level variables hard, partial rebound.
  if (isTRUE(cfg$synthetic$covid_break)) {
    tc <- cfg$synthetic$covid_t
    if (tc > 2 && tc < n - 1) {
      shock <- rep(0, M)
      shock[spec$variable == "gdp_growth"]  <- -7
      shock[spec$variable == "cons_growth"] <- -10
      shock[spec$variable == "emp_growth"]  <- -5
      shock[spec$variable == "unemp_rate"]  <- 3
      shock[spec$variable == "f_act"]       <- -8
      shock[spec$variable == "f_act_tw"]    <- -8
      shock[spec$variable == "cash_rate"]   <- -1.5
      y[tc, ] <- y[tc, ] + shock
      y[tc + 1, ] <- y[tc + 1, ] - 0.7 * shock   # rebound
    }
  }

  colnames(y) <- spec$variable

  # Integrate transformed series back to raw levels.
  raw <- as.data.frame(y)
  for (i in seq_len(M)) {
    v <- spec$variable[i]
    if (spec$transform[i] == "dlog") {
      raw[[v]] <- 100 * exp(cumsum(y[, i] / 100))      # index, base 100
    } else if (spec$transform[i] == "loglevel") {
      raw[[v]] <- exp(y[, i] / 100)
    }                                                   # level: as-is
  }
  dates <- quarter_seq(cfg$data$start, n)
  out <- list(
    data = cbind(data.frame(date = dates), raw),
    manifest = data.frame(variable = spec$variable, provider = "synthetic",
                          series_id = "synthetic", pulled = as.character(Sys.Date()))
  )
  attr(out, "dgp") <- list(A1 = A1, A2 = A2, G = G, mu = mu, nf = nf, seed = seed)
  log_info("Synthetic data generated: {n} quarters x {M} variables (seed {seed})")
  out
}

# ---- real ----------------------------------------------------------------------
#
# The foreign block is FRED-primary (GDPC1, FEDFUNDS) and Australian data comes
# from RBA/ABS. There is deliberately NO silent series substitution: a download
# failure errors loudly rather than swapping in a different concept.

#' Pull one series from its provider; returns data.frame(date, value) at native
#' frequency. Errors are caught by the caller.
fetch_series <- function(v, src, cfg) {
  provider <- src$provider
  if (provider == "rba") {
    x <- readrba::read_rba(series_id = src$id)
    data.frame(date = as.Date(x$date), value = x$value)
  } else if (provider == "abs") {
    x <- readabs::read_abs(series_id = src$id)
    data.frame(date = as.Date(x$date), value = x$value)
  } else if (provider == "fred") {
    # key is validated and set once in download_real_data
    x <- fredr::fredr(series_id = src$id)
    x <- x[!is.na(x$value), ]
    if (nrow(x) == 0) stop("FRED series '", src$id, "' returned no observations")
    data.frame(date = as.Date(x$date), value = x$value)
  } else if (provider == "dbnomics") {
    x <- rdbnomics::rdb(ids = src$id)
    x <- x[!is.na(x$value), ]
    data.frame(date = as.Date(x$period), value = x$value)
  } else stop("unknown provider: ", provider)
}

#' Quarterly aggregation: average of observations within the quarter, but
#' COVERAGE-AWARE -- a quarter with materially fewer source observations than a
#' full quarter (a partial trailing/leading quarter at the data frontier) is
#' returned as NA rather than a 1-2-month/half-quarter average that would
#' silently enter the panel. NA partial quarters are then dropped (if trailing)
#' or error (if interior) downstream, never silently averaged.
to_quarterly <- function(df, min_coverage = 0.8) {
  df <- df[is.finite(df$value), , drop = FALSE]
  q <- as.Date(cut(df$date, "quarter"))
  m <- tapply(df$value, q, mean)
  n <- tapply(df$value, q, length)
  dates <- as.Date(names(m))
  expected <- max(n)                       # fullest quarter (1 qtrly, 3 mthly, ~63 daily)
  value <- as.numeric(m)
  value[as.integer(n) < min_coverage * expected] <- NA_real_
  data.frame(date = dates, value = value)
}

#' Download all series with on-disk caching; fall back per-series to FRED if a
#' key is available. Series already published as quarterly % change (trimmed
#' mean CPI) are cumulated into a level index so the transform layer is uniform.
download_real_data <- function(cfg, spec, raw_dir = "data/raw") {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(raw_dir, "manifest.csv")
  rows <- list(); series <- list()
  if (any(spec$provider == "fred")) {
    if (!requireNamespace("fredr", quietly = TRUE))
      stop("package 'fredr' is required for FRED series; run renv::restore().")
    key <- Sys.getenv(cfg$data$fred_api_key_env, "")
    if (!nzchar(key))
      stop("FRED series configured but ", cfg$data$fred_api_key_env,
           " is not set. Put it in .Renviron (gitignored).")
    fredr::fredr_set_key(key)
  }

  for (i in seq_len(nrow(spec))) {
    v <- spec$variable[i]
    src <- list(provider = spec$provider[i], id = spec$series_id[i])
    cache_file <- file.path(raw_dir, paste0(v, ".rds"))
    # cache is keyed by variable name, so it must also record the provider+id
    # it was pulled with -- otherwise changing a series ID silently serves
    # stale data (a production footgun). A cache hit requires fresh AND a
    # matching (provider, id).
    cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
    cache_ok <- !is.null(cached) && is.list(cached) &&
      identical(cached$provider, src$provider) &&
      identical(cached$id, src$id) &&
      difftime(Sys.time(), file.mtime(cache_file), units = "days") <
        cfg$data$cache_max_age_days
    if (cache_ok) {
      log_info("cache hit: {v} ({src$provider}:{src$id})")
      series[[v]] <- cached$data
    } else {
      df <- tryCatch(fetch_series(v, src, cfg), error = function(e) {
        log_warn("download failed for {v} ({src$provider}:{src$id}): {conditionMessage(e)}")
        NULL
      })
      if (is.null(df)) stop("Could not obtain series '", v, "' (", src$provider, ":",
                            src$id, "). Fix the source or set data.source: synthetic to run offline.")
      saveRDS(list(provider = src$provider, id = src$id, data = df,
                   pulled = as.character(Sys.Date())), cache_file)
      series[[v]] <- df
    }
    rows[[v]] <- data.frame(variable = v, provider = src$provider,
                            series_id = src$id,
                            last_obs = as.character(max(series[[v]]$date)),
                            pulled = as.character(Sys.Date()))
  }

  # quarterly panel on a common date index
  qs <- lapply(names(series), function(v) {
    q <- to_quarterly(series[[v]])
    # series published as a qtr % change (e.g. trimmed-mean CPI, G20 GDP
    # growth) are cumulated to an index so the uniform dlog transform applies;
    # cumprod treats the published number as an arithmetic % change (exact),
    # the subsequent dlog yields log growth consistent with the other series.
    # cumprod propagates any interior NA across the whole tail, which would
    # silently shorten the panel -- so require a contiguous finite support and
    # error loudly on an interior gap (never silently truncate).
    if (identical(spec$pre[spec$variable == v], "pct_change")) {
      stopifnot(identical(spec$transform[spec$variable == v], "dlog"))
      fin <- which(is.finite(q$value))
      if (length(fin)) {
        span <- min(fin):max(fin)
        if (anyNA(q$value[span]))
          stop("interior NA in pct_change series '", v, "' before cumulation")
        q$value[span] <- 100 * cumprod(1 + q$value[span] / 100)
      }
    }
    names(q)[2] <- v
    q
  })
  panel <- Reduce(function(a, b) merge(a, b, by = "date", all = TRUE), qs)
  panel <- panel[panel$date >= as.Date(cfg$data$start) &
                 panel$date <= as.Date(cfg$data$end), ]
  manifest <- do.call(rbind, rows)
  write.csv(manifest, manifest_path, row.names = FALSE)
  log_info("real data panel: {nrow(panel)} quarters x {ncol(panel)-1} series")
  list(data = panel, manifest = manifest)
}

#' Entry point used by the pipeline.
get_raw_data <- function(cfg, spec) {
  if (identical(cfg$data$source, "real")) {
    download_real_data(cfg, spec)
  } else {
    generate_synthetic_data(cfg, spec)
  }
}
