# utils.R -- config, logging, seeds, transform_spec, shared helpers ------------

suppressPackageStartupMessages({
  library(yaml)
  library(logger)
  library(digest)
})

#' Load the project configuration.
load_config <- function(path = "config/config.yml") {
  cfg <- yaml::read_yaml(path)$default
  stopifnot(!is.null(cfg$master_seed), !is.null(cfg$variables))
  cfg
}

#' Hash of the config sections that affect estimation, used to key the OOS cache.
config_hash <- function(cfg) {
  digest::digest(cfg[c("master_seed", "data", "synthetic", "variables",
                       "horizons", "mcmc", "glp", "suite", "benchmarks",
                       "evaluation")], algo = "xxhash64")
}

#' Set up the logger once per session.
setup_logging <- function(level = "INFO") {
  logger::log_threshold(level)
  logger::log_layout(logger::layout_glue_generator(
    format = "{format(time, '%H:%M:%S')} [{level}] {msg}"))
  invisible(NULL)
}

#' Deterministic per-task seed derived from the master seed and a string key,
#' so parallel workers are reproducible regardless of scheduling.
derive_seed <- function(master_seed, key) {
  h <- digest::digest(paste0(master_seed, "::", key), algo = "xxhash32")
  (strtoi(substr(h, 1, 7), base = 16L) + master_seed) %% .Machine$integer.max
}

#' Time an expression and log its duration.
timed <- function(label, expr) {
  t0 <- Sys.time()
  res <- force(expr)
  logger::log_info("{label} took {round(as.numeric(difftime(Sys.time(), t0, units = 'secs')), 1)}s")
  res
}

# transform_spec ---------------------------------------------------------------

#' Build the transform_spec table (single source of truth) from config.
#' Rows are ordered foreign block first, then domestic, preserving config order.
build_transform_spec <- function(cfg) {
  rows <- lapply(names(cfg$variables), function(v) {
    x <- cfg$variables[[v]]
    data.frame(
      variable   = v,
      label      = x$label,
      block      = x$block,
      transform  = x$transform,
      delta      = x$delta,
      ss_mean    = if (is.null(x$ss_mean)) NA_real_ else x$ss_mean,
      ss_sd      = x$ss_sd,
      sets       = paste(unlist(x$sets), collapse = ","),
      provider   = x$source$provider,
      series_id  = x$source$id,
      target     = isTRUE(x$target),
      year_ended = isTRUE(x$year_ended),
      stringsAsFactors = FALSE
    )
  })
  spec <- do.call(rbind, rows)
  spec <- spec[order(match(spec$block, c("foreign", "domestic"))), ]
  rownames(spec) <- NULL
  stopifnot(all(spec$block %in% c("foreign", "domestic")),
            all(spec$transform %in% c("dlog", "level", "loglevel")))
  spec
}

#' Variables belonging to a model set, foreign first.
vars_for_set <- function(spec, set) {
  keep <- vapply(strsplit(spec$sets, ","),
                 function(s) set %in% s, logical(1))
  spec[keep, , drop = FALSE]
}

# linear algebra helpers --------------------------------------------------------

#' Stable log-sum-exp.
logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Draw from N(mu, Sigma) given a precomputed upper Cholesky of Sigma.
rmvnorm_chol <- function(mu, chol_sigma) {
  drop(mu + crossprod(chol_sigma, rnorm(length(mu))))
}

#' Draw from an inverse-Wishart IW(S, nu): if X ~ W(S^-1, nu) then X^-1 ~ IW(S, nu).
riwish <- function(nu, S) {
  m <- nrow(S)
  stopifnot(nu > m - 1)
  # Bartlett decomposition of the Wishart W(S^{-1}, nu)
  L <- t(chol(chol2inv(chol(S))))   # S^{-1} = L L'
  A <- matrix(0, m, m)
  diag(A) <- sqrt(rchisq(m, df = nu - seq_len(m) + 1))
  A[lower.tri(A)] <- rnorm(m * (m - 1) / 2)
  W <- L %*% A %*% t(A) %*% t(L)
  chol2inv(chol(W))
}

#' Companion matrix of a VAR coefficient matrix B (K x M, rows: intercept then lags).
companion_matrix <- function(B, M, p) {
  A <- t(B[-1, , drop = FALSE])               # M x (M*p), drop intercept row
  if (p == 1) return(A)
  rbind(A, cbind(diag(M * (p - 1)), matrix(0, M * (p - 1), M)))
}

#' Max eigenvalue modulus of the companion matrix.
max_eig_mod <- function(B, M, p) {
  max(Mod(eigen(companion_matrix(B, M, p), only.values = TRUE)$values))
}

#' Build lagged regressor matrix: X[t,] = (1, y_{t-1}, ..., y_{t-p}); returns
#' list(Y = y[(p+1):T_n, ], X) so that Y = X B + E.
build_XY <- function(y, p, intercept = TRUE) {
  T_n <- nrow(y); M <- ncol(y)
  stopifnot(T_n > p)
  X <- matrix(NA_real_, T_n - p, M * p)
  for (l in seq_len(p)) {
    X[, ((l - 1) * M + 1):(l * M)] <- y[(p + 1 - l):(T_n - l), , drop = FALSE]
  }
  if (intercept) X <- cbind(1, X)
  list(Y = y[(p + 1):T_n, , drop = FALSE], X = X)
}

#' Per-variable residual sd from a univariate AR(plag) fit -- the sigma_i used
#' to scale the Minnesota prior.
ar_sigmas <- function(y, plag = 4) {
  apply(y, 2, function(z) {
    z <- z[is.finite(z)]
    p <- min(plag, floor(length(z) / 4))
    fit <- tryCatch(stats::ar(z, order.max = p, aic = FALSE, method = "ols"),
                    error = function(e) NULL)
    if (is.null(fit) || !is.finite(fit$var.pred) || fit$var.pred <= 0)
      stats::sd(z) else sqrt(fit$var.pred)
  })
}

#' Quarterly date sequence helper.
quarter_seq <- function(start, n) seq(as.Date(start), by = "quarter", length.out = n)
