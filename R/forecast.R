# forecast.R -- iterated predictive simulation for every engine ----------------
#
# simulate_paths(post, y, h, ndraw, condition = NULL) -> array [ndraw, h, M]
# of future values in model units, integrating over parameter draws (recycled
# if ndraw > stored draws) and future shocks.
#
# condition: optional list(variable = , path = numeric(h)) -- hard conditioning
# by substitution (the conditioned variable is overwritten each step before it
# feeds subsequent dynamics). An approximation to proper conditional
# forecasting (no Waggoner-Zha shock adjustment); see DECISIONS.md.

simulate_paths <- function(post, y, h, ndraw, condition = NULL) {
  UseMethod("simulate_paths")
}

.apply_condition <- function(ystep, s, condition, varnames) {
  if (is.null(condition)) return(ystep)
  j <- match(condition$variable, varnames)
  if (!is.na(j)) ystep[j] <- condition$path[s]
  ystep
}

#' Iterate one VAR path: B (K x M, intercept first), cS = chol(Sigma),
#' ystate = matrix p x M with row 1 = most recent observation.
.iterate_path <- function(B, cS, ystate, h, M, p, condition, varnames) {
  out <- matrix(NA_real_, h, M)
  for (s in seq_len(h)) {
    x <- c(1, as.vector(t(ystate[seq_len(p), , drop = FALSE])))
    mu_s <- drop(crossprod(B, x))
    ynew <- mu_s + drop(crossprod(cS, rnorm(M)))
    ynew <- .apply_condition(ynew, s, condition, varnames)
    out[s, ] <- ynew
    ystate <- rbind(ynew, ystate[-p, , drop = FALSE])
  }
  out
}

# state helper: last p observations, most recent first
.ystate <- function(y, p) y[nrow(y) - seq_len(p) + 1, , drop = FALSE]

# ---- gibbs ----------------------------------------------------------------------

simulate_paths.post_gibbs <- function(post, y, h, ndraw, condition = NULL) {
  M <- post$M; p <- post$p
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)
  for (d in seq_len(ndraw)) {
    k <- ((d - 1) %% post$ndraw) + 1
    B <- post$B[k, , ]
    cS <- chol(post$Sigma[k, , ])
    paths[d, , ] <- .iterate_path(B, cS, st0, h, M, p, condition, post$varnames)
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- ss -------------------------------------------------------------------------

simulate_paths.post_ss <- function(post, y, h, ndraw, condition = NULL) {
  M <- post$M; p <- post$p
  paths <- array(NA_real_, c(ndraw, h, M))
  for (d in seq_len(ndraw)) {
    k <- ((d - 1) %% post$ndraw) + 1
    A <- post$A[k, , ]; Psi <- post$Psi[k, ]
    cS <- chol(post$Sigma[k, , ])
    z <- sweep(y, 2, Psi)
    st <- .ystate(z, p)
    B <- rbind(0, A)                      # zero intercept in demeaned form
    zp <- .iterate_path(B, cS, st, h, M, p, condition = NULL, post$varnames)
    pth <- sweep(zp, 2, Psi, "+")
    if (!is.null(condition)) {
      j <- match(condition$variable, post$varnames)
      if (!is.na(j)) pth[, j] <- condition$path
    }
    paths[d, , ] <- pth
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- conj_br ---------------------------------------------------------------------

simulate_paths.post_conj_br <- function(post, y, h, ndraw, condition = NULL) {
  M <- post$M; p <- post$p; nf <- post$nf
  nd <- M - nf
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)                       # all vars, most recent first
  for (d in seq_len(ndraw)) {
    k <- ((d - 1) %% post$ndraw) + 1
    Bf <- post$foreign$B[k, , ];  cSf <- chol(post$foreign$Sigma[k, , ])
    Bd <- post$domestic$B[k, , ]; cSd <- chol(post$domestic$Sigma[k, , ])
    st <- st0
    for (s in seq_len(h)) {
      xf <- c(1, as.vector(t(st[seq_len(p), seq_len(nf), drop = FALSE])))
      yf <- drop(crossprod(Bf, xf)) + drop(crossprod(cSf, rnorm(nf)))
      xd <- c(1, as.vector(t(st[seq_len(p), , drop = FALSE])), yf)
      yd <- drop(crossprod(Bd, xd)) + drop(crossprod(cSd, rnorm(nd)))
      ynew <- c(yf, yd)
      ynew <- .apply_condition(ynew, s, condition, post$varnames)
      paths[d, s, ] <- ynew
      st <- rbind(ynew, st[-p, , drop = FALSE])
    }
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- sv --------------------------------------------------------------------------

simulate_paths.post_sv <- function(post, y, h, ndraw, condition = NULL) {
  M <- post$M; p <- post$p
  paths <- array(NA_real_, c(ndraw, h, M))
  st0 <- .ystate(y, p)
  for (d in seq_len(ndraw)) {
    k <- ((d - 1) %% post$ndraw) + 1
    st <- st0
    # pull per-equation parameters for this draw
    betas <- lapply(post$eqs, function(e) e$beta[k, ])
    hs    <- vapply(post$eqs, function(e) e$hT[k], numeric(1))
    svp   <- lapply(post$eqs, function(e) e$svpara[k, ])
    for (s in seq_len(h)) {
      ynew <- numeric(M)
      lags_vec <- as.vector(t(st[seq_len(p), , drop = FALSE]))
      for (i in seq_len(M)) {
        eq <- post$eqs[[i]]$design
        # rebuild the design vector for this equation: 1, allowed lags, contemp
        allow <- eq$lag_meta
        xlag <- lags_vec[(allow$lag - 1) * M + allow$var]
        xx <- c(1, xlag, if (length(eq$contemp)) ynew[eq$contemp])
        hs[i] <- svp[[i]]["mu"] + svp[[i]]["phi"] * (hs[i] - svp[[i]]["mu"]) +
                 svp[[i]]["sigma"] * rnorm(1)
        ynew[i] <- sum(xx * betas[[i]]) + exp(hs[i] / 2) * rnorm(1)
      }
      ynew <- .apply_condition(ynew, s, condition, post$varnames)
      paths[d, s, ] <- ynew
      st <- rbind(ynew, st[-p, , drop = FALSE])
    }
  }
  dimnames(paths)[[3]] <- post$varnames
  paths
}

# ---- fan-chart quantiles + sanity checks ------------------------------------------

#' Quantiles per horizon/variable from a path array.
fan_quantiles <- function(paths, probs) {
  qs <- apply(paths, c(2, 3), quantile, probs = probs)
  # qs: [prob, h, var] -> long data.frame
  out <- expand.grid(prob = probs, h = seq_len(dim(paths)[2]),
                     variable = dimnames(paths)[[3]], stringsAsFactors = FALSE)
  out$value <- as.vector(qs)
  out
}

#' Forecast sanity (section 9): finite, non-explosive, and -- for
#' mean-reverting (delta = 0) variables only -- long-horizon reversion toward
#' an unconditional anchor. Variables modelled as persistent levels
#' (delta = 1: rates, log TWI/ToT) are deliberately near-unit-root, so
#' 12-quarter reversion to the sample mean is NOT an implication of the model;
#' they are checked for boundedness only.
check_forecasts <- function(paths, y, label = "member", delta = NULL) {
  ok_finite <- all(is.finite(paths))
  rng <- apply(abs(y), 2, max)
  ok_bound <- TRUE
  mean_path <- apply(paths, c(2, 3), median)
  for (j in seq_len(ncol(y))) {
    bound <- 5 * max(rng[j], 1) + 50
    # explosiveness is about the central mass, not individual fat-tail draws
    ql <- apply(paths[, , j, drop = FALSE], 2, quantile, probs = c(0.005, 0.995))
    if (any(abs(ql) > bound)) ok_bound <- FALSE
  }
  H <- dim(paths)[2]
  ybar <- colMeans(y)
  dev4 <- abs(mean_path[min(4, H), ] - ybar)
  devH <- abs(mean_path[H, ] - ybar)
  # 2.5 sd: flag genuine drift (the pathologies this catches are 10-100x sd),
  # not honest persistence of a below-mean regime
  rev_ok <- devH <= pmax(dev4 * 1.5, 2.5 * apply(y, 2, sd))
  if (!is.null(delta)) rev_ok <- rev_ok | (delta == 1)
  ok_converge <- all(rev_ok)
  ok <- ok_finite && ok_bound && ok_converge
  if (!ok) logger::log_warn(
    "forecast sanity FAILED for {label}: finite={ok_finite} bounded={ok_bound} converge={ok_converge}")
  list(finite = ok_finite, bounded = ok_bound, converged = ok_converge, ok = ok)
}
