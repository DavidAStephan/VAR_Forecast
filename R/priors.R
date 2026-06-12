# priors.R -- Minnesota, block-exogeneity, dummy observations, steady-state ----
#
# Conventions: the VAR in regression form is Y = X B + E with
# X[t,] = (1, y_{t-1}, ..., y_{t-p}), B is K x M, K = 1 + M*p.
# Coefficient (row k>1, col i) is the loading of equation i on lag l of
# variable j, where k = 1 + (l-1)*M + j.

#' Index helper: for each non-intercept row of B, which (lag, variable) it is.
coef_index <- function(M, p) {
  data.frame(row = 2:(1 + M * p),
             lag = rep(seq_len(p), each = M),
             var = rep(seq_len(M), times = p))
}

#' Independent-prior Minnesota moments (mean b0 and sd s0, both K x M), with
#' block exogeneity imposed: near-zero mean AND near-zero sd on domestic-lag
#' coefficients inside foreign equations.
#'
#' sigma: M-vector of AR residual sds. delta: M-vector of own-lag prior means.
#' blocks: character M-vector ("foreign"/"domestic").
minnesota_prior <- function(M, p, sigma, delta, blocks,
                            lambda = 0.2, cross = 0.5, lag_decay = 1,
                            intercept_scale = 10,
                            block_exog_sd = 1e-5) {
  K <- 1 + M * p
  b0 <- matrix(0, K, M)
  s0 <- matrix(0, K, M)
  idx <- coef_index(M, p)
  b0[1, ] <- 0
  s0[1, ] <- intercept_scale * sigma
  for (r in seq_len(nrow(idx))) {
    k <- idx$row[r]; l <- idx$lag[r]; j <- idx$var[r]
    for (i in seq_len(M)) {
      own <- (i == j)
      if (own && l == 1) b0[k, i] <- delta[i]
      s0[k, i] <- lambda * (if (own) 1 else cross) *
        (sigma[i] / sigma[j]) / l^lag_decay
    }
  }
  # block exogeneity: foreign equations must not load on domestic lags
  for_eq <- which(blocks == "foreign")
  dom_var <- which(blocks == "domestic")
  if (length(for_eq) && length(dom_var)) {
    rows <- idx$row[idx$var %in% dom_var]
    b0[rows, for_eq] <- 0
    s0[rows, for_eq] <- block_exog_sd
  }
  list(b0 = b0, s0 = s0)
}

#' Conjugate (Kronecker) Minnesota: prior B|Sigma ~ MN(B0, Sigma x Omega0),
#' Sigma ~ IW(S0, nu0). Omega0 is K x K diagonal; cross-shrinkage is forced to
#' 1 by the Kronecker structure (which is why block exogeneity needs the
#' independent-prior engines).
conjugate_prior <- function(M, p, sigma, delta,
                            lambda = 0.2, lag_decay = 1, intercept_scale = 10) {
  K <- 1 + M * p
  B0 <- matrix(0, K, M)
  B0[1 + seq_len(M), ] <- diag(delta, M)        # own first lag
  omega <- numeric(K)
  omega[1] <- intercept_scale^2
  idx <- coef_index(M, p)
  omega[idx$row] <- (lambda / (sigma[idx$var] * idx$lag^lag_decay))^2
  nu0 <- M + 2
  S0  <- diag(sigma^2, M) * (nu0 - M - 1)        # prior mean of Sigma = diag(sigma^2)
  list(B0 = B0, Omega0_diag = omega, S0 = S0, nu0 = nu0)
}

#' Sum-of-coefficients (Doan-Litterman-Sims) and dummy-initial-observation
#' (Sims-Zha) artificial observations, appended to (Y, X).
#' ybar: M-vector of pre-sample means (mean of the first p observations).
soc_dio_dummies <- function(M, p, ybar, delta,
                            soc = TRUE, soc_mu = 1, dio = TRUE, dio_delta = 1) {
  Yd <- NULL; Xd <- NULL
  if (soc && soc_mu > 0) {
    # one dummy row per variable with delta=1: y*_j = ybar_j / mu on both sides
    keep <- which(delta == 1)
    if (length(keep)) {
      Ys <- matrix(0, length(keep), M)
      Xs <- matrix(0, length(keep), 1 + M * p)
      for (r in seq_along(keep)) {
        j <- keep[r]
        Ys[r, j] <- ybar[j] / soc_mu
        Xs[r, 1 + j + (seq_len(p) - 1) * M] <- ybar[j] / soc_mu
      }
      Yd <- rbind(Yd, Ys); Xd <- rbind(Xd, Xs)
    }
  }
  if (dio && dio_delta > 0) {
    Yi <- matrix(ybar / dio_delta, 1, M)
    Xi <- matrix(0, 1, 1 + M * p)
    Xi[1, 1] <- 1 / dio_delta
    Xi[1, -1] <- rep(ybar, p) / dio_delta
    Yd <- rbind(Yd, Yi); Xd <- rbind(Xd, Xi)
  }
  list(Y = Yd, X = Xd)
}

#' Steady-state prior moments from the transform_spec: N(psi0, diag(psi_sd^2)).
#' NA anchors fall back to the estimation-window sample mean with the spec sd.
steadystate_prior <- function(spec_m, y) {
  psi0 <- ifelse(is.na(spec_m$ss_mean), colMeans(y), spec_m$ss_mean)
  list(psi0 = psi0, psi_sd = spec_m$ss_sd)
}

# ---- marginal likelihood (conjugate NIW) for GLP-style lambda selection -------

lmvgamma <- function(a, m) {
  m * (m - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(m)) / 2))
}

#' Closed-form log marginal likelihood of the conjugate NIW BVAR.
log_marginal_likelihood <- function(Y, X, prior) {
  T_n <- nrow(Y); M <- ncol(Y)
  Oi  <- diag(1 / prior$Omega0_diag)
  XtX <- crossprod(X)
  P1  <- Oi + XtX
  cP1 <- chol(P1)
  B1  <- backsolve(cP1, forwardsolve(t(cP1), Oi %*% prior$B0 + crossprod(X, Y)))
  S1  <- prior$S0 + crossprod(Y) + t(prior$B0) %*% Oi %*% prior$B0 - t(B1) %*% P1 %*% B1
  S1  <- (S1 + t(S1)) / 2
  nu1 <- prior$nu0 + T_n
  ldet <- function(A) as.numeric(determinant(A, logarithm = TRUE)$modulus)
  -(M * T_n / 2) * log(pi) +
    lmvgamma(nu1 / 2, M) - lmvgamma(prior$nu0 / 2, M) +
    (prior$nu0 / 2) * ldet(prior$S0) - (nu1 / 2) * ldet(S1) +
    (M / 2) * (-2 * sum(log(diag(cP1))) - sum(log(prior$Omega0_diag)))
}

#' GLP-style data-driven overall tightness: grid-search the conjugate marginal
#' likelihood on the given data set. Returns the selected lambda. Optional row
#' weights apply the LP COVID rescaling (the s_t Jacobian is constant in
#' lambda, so it drops out of the comparison).
select_lambda <- function(y, p, sigma, delta, grid, lag_decay = 1,
                          weights = NULL) {
  xy <- build_XY(y, p)
  if (!is.null(weights)) {
    w <- weights[(p + 1):nrow(y)]
    xy$Y <- xy$Y * w; xy$X <- xy$X * w
  }
  M <- ncol(y)
  mls <- vapply(grid, function(lam) {
    pr <- conjugate_prior(M, p, sigma, delta, lambda = lam, lag_decay = lag_decay)
    log_marginal_likelihood(xy$Y, xy$X, pr)
  }, numeric(1))
  lam <- grid[which.max(mls)]
  logger::log_debug("GLP lambda grid: {paste(round(mls,1), collapse=' ')} -> lambda={lam}")
  lam
}
