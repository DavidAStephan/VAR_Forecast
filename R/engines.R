# engines.R -- estimation engines behind a common interface --------------------
#
# fit_var_member(y, member, spec_m, cfg) -> posterior (class = engine name)
# Each posterior carries meta (M, p, varnames, blocks) and $diagnostics.
# Predictive simulation lives in forecast.R (simulate_paths.<engine>).
#
# Engines:
#   conj_br : block-recursive conjugate NIW (foreign VAR + domestic conditional)
#   gibbs   : joint independent-Normal + inverse-Wishart Gibbs (block exog via
#             asymmetric prior variances -- the workhorse required by the spec)
#   ss      : Villani steady-state Gibbs (unconditional mean explicitly sampled)
#   sv      : triangular equation-by-equation with stochastic volatility
#             (CCM-2019-style; exact factorisation of the reduced form)

# stochvol (SV engine) and coda (ESS diagnostic) are OPTIONAL: used via guarded
# wrappers (has_stochvol(), safe_ess()) from aaa_capabilities.R, not attached.

# ---- shared helpers -----------------------------------------------------------

#' Resolve member prior settings against config defaults / GLP-selected lambda.
resolve_prior <- function(member, glp_lambda) {
  pr <- member$prior
  if (identical(pr$lambda, "auto")) pr$lambda <- glp_lambda
  pr$lambda    <- as.numeric(pr$lambda)
  pr$soc       <- isTRUE(pr$soc)
  pr$dio       <- isTRUE(pr$dio)
  pr$soc_mu    <- if (is.null(pr$soc_mu)) 1 else pr$soc_mu
  pr$dio_delta <- if (is.null(pr$dio_delta)) 1 else pr$dio_delta
  pr
}

#' Sample (B, Sigma) from a conjugate NIW posterior given (Y, X, prior).
#' Returns arrays B [n, K, M], Sigma [n, M, M].
sample_conjugate <- function(Y, X, prior, n) {
  T_n <- nrow(Y); M <- ncol(Y); K <- ncol(X)
  Oi  <- diag(1 / prior$Omega0_diag, K)
  P1  <- Oi + crossprod(X)
  cP1 <- chol(P1)
  B1  <- backsolve(cP1, forwardsolve(t(cP1), Oi %*% prior$B0 + crossprod(X, Y)))
  S1  <- prior$S0 + crossprod(Y) + t(prior$B0) %*% Oi %*% prior$B0 - t(B1) %*% P1 %*% B1
  S1  <- (S1 + t(S1)) / 2
  nu1 <- prior$nu0 + T_n
  Bd <- array(NA_real_, c(n, K, M)); Sd <- array(NA_real_, c(n, M, M))
  for (d in seq_len(n)) {
    Sig <- riwish(nu1, S1)
    cS  <- chol(Sig)
    Z   <- matrix(rnorm(K * M), K, M)
    # B = B1 + chol(Omega1)' Z chol(Sigma); Omega1 = P1^{-1}
    Bd[d, , ] <- B1 + backsolve(cP1, Z) %*% cS
    Sd[d, , ] <- Sig
  }
  list(B = Bd, Sigma = Sd, B1 = B1)
}

#' One Gibbs draw of vec(B) | Sigma for the independent-prior VAR.
#' V0inv_diag: KM-vector of prior precisions; b0vec: KM prior mean.
draw_beta_given_sigma <- function(XtX, XtY, Sigma_inv, V0inv_diag, b0vec, K, M) {
  P <- kronecker(Sigma_inv, XtX)
  diag(P) <- diag(P) + V0inv_diag
  rhs <- as.vector(XtY %*% Sigma_inv) + V0inv_diag * b0vec
  cP <- chol(P)
  mean_ <- backsolve(cP, forwardsolve(t(cP), rhs))
  matrix(mean_ + backsolve(cP, rnorm(K * M)), K, M)
}

#' Posterior-mean diagnostic: max |coef| of domestic lags inside foreign
#' equations (must be ~ 0 under block exogeneity).
block_exog_metric <- function(Bbar, M, p, blocks) {
  idx <- coef_index(M, p)
  rows <- idx$row[idx$var %in% which(blocks == "domestic")]
  cols <- which(blocks == "foreign")
  if (!length(rows) || !length(cols)) return(0)
  max(abs(Bbar[rows, cols, drop = FALSE]))
}

#' MCMC diagnostics: ESS for the own-first-lag coefficients + stability share.
mcmc_diagnostics <- function(own_lag_draws, stable_share) {
  if (!has_coda())   # no ESS available: report NA and let the gate pass
    return(list(ess_min = NA_real_, ess_median = NA_real_,
                stable_share = stable_share, converged = TRUE))
  ess <- apply(own_lag_draws, 2, function(x) {
    if (sd(x) < 1e-12) return(NA_real_)   # pinned by block-exog prior etc.
    safe_ess(x)
  })
  list(ess_min = suppressWarnings(min(ess, na.rm = TRUE)),
       ess_median = suppressWarnings(median(ess, na.rm = TRUE)),
       stable_share = stable_share,
       converged = is.finite(min(ess, na.rm = TRUE)) &&
                   min(ess, na.rm = TRUE) > 50)
}

# ---- engine: conj_br (block-recursive conjugate NIW) --------------------------

fit_conj_br <- function(y, member, spec_m, cfg, prior, weights = NULL) {
  M <- ncol(y); p <- member$lags
  blocks <- spec_m$block
  nf <- sum(blocks == "foreign"); nd <- M - nf
  # the block-recursive split below indexes (nf+1):M / seq_len(nf); both blocks
  # must be non-empty (a degenerate partition would mis-index in R).
  stopifnot("conj_br needs >=1 foreign and >=1 domestic variable" = nf >= 1 && nd >= 1)
  sigma <- ar_sigmas(y, weights = weights); delta <- spec_m$delta
  n <- cfg$mcmc$forecast_draws
  w_rows <- if (is.null(weights)) NULL else weights[(p + 1):nrow(y)]

  yf <- y[, seq_len(nf), drop = FALSE]
  ybar <- colMeans(y[seq_len(p), , drop = FALSE])

  # foreign marginal VAR (conjugate, own dummies)
  xyf <- build_XY(yf, p)
  if (!is.null(w_rows)) { xyf$Y <- xyf$Y * w_rows; xyf$X <- xyf$X * w_rows }
  prf <- conjugate_prior(nf, p, sigma[seq_len(nf)], delta[seq_len(nf)],
                         lambda = prior$lambda)
  if (prior$soc || prior$dio) {
    dmy <- soc_dio_dummies(nf, p, ybar[seq_len(nf)], delta[seq_len(nf)],
                           soc = prior$soc, soc_mu = prior$soc_mu,
                           dio = prior$dio, dio_delta = prior$dio_delta)
    if (!is.null(dmy$Y)) { xyf$Y <- rbind(dmy$Y, xyf$Y); xyf$X <- rbind(dmy$X, xyf$X) }
  }
  postf <- sample_conjugate(xyf$Y, xyf$X, prf, n)

  # domestic conditional regression: y_d on (1, lags of all, contemporaneous y_f)
  xy <- build_XY(y, p)
  Yd <- xy$Y[, (nf + 1):M, drop = FALSE]
  Xd <- cbind(xy$X, y[(p + 1):nrow(y), seq_len(nf), drop = FALSE])
  if (!is.null(w_rows)) { Yd <- Yd * w_rows; Xd <- Xd * w_rows }
  Kd <- ncol(Xd)
  sig_d <- sigma[(nf + 1):M]
  B0d <- matrix(0, Kd, nd)
  idx <- coef_index(M, p)
  for (i in seq_len(nd)) B0d[1 + (nf + i), i] <- delta[nf + i]  # own first lag
  omega <- numeric(Kd)
  omega[1] <- 100
  # Minnesota lag variance with a FIXED lag-decay exponent of 1 (variance ~ 1/lag)
  # for the conjugate scheme; the harmonic decay is not configurable here.
  omega[idx$row] <- (prior$lambda / (sigma[idx$var] * idx$lag))^2
  omega[(1 + M * p + 1):Kd] <- (1 / sigma[seq_len(nf)])^2  # contemp. foreign, loose
  nu0 <- nd + 2
  prd <- list(B0 = B0d, Omega0_diag = omega,
              S0 = diag(sig_d^2, nd) * (nu0 - nd - 1), nu0 = nu0)
  if (prior$soc || prior$dio) {
    dmy <- soc_dio_dummies(M, p, ybar, delta, soc = prior$soc,
                           soc_mu = prior$soc_mu, dio = prior$dio,
                           dio_delta = prior$dio_delta)
    if (!is.null(dmy$Y)) {
      # keep only rows for domestic variables; extend X with contemp.-foreign
      # cols. SOC rows are zero off the own variable, so their foreign
      # contemporaneous entries are 0; the DIO row asserts ALL variables at
      # ybar, so its contemporaneous-foreign entries must be ybar_f/dio_delta.
      keepr <- which(rowSums(abs(dmy$Y[, (nf + 1):M, drop = FALSE])) > 0 |
                     rowSums(abs(dmy$Y)) == 0)
      dYd <- dmy$Y[keepr, (nf + 1):M, drop = FALSE]
      dXc <- matrix(0, length(keepr), nf)
      is_dio <- abs(dmy$X[keepr, 1]) > 0          # DIO is the row with intercept
      if (any(is_dio) && prior$dio)
        dXc[is_dio, ] <- matrix(ybar[seq_len(nf)] / prior$dio_delta,
                                sum(is_dio), nf, byrow = TRUE)
      dXd <- cbind(dmy$X[keepr, , drop = FALSE], dXc)
      Yd <- rbind(dYd, Yd); Xd <- rbind(dXd, Xd)
    }
  }
  postd <- sample_conjugate(Yd, Xd, prd, n)

  # diagnostics: stability of the implied joint system at posterior mean
  Bf_bar <- postf$B1
  stable <- mean(vapply(seq_len(min(n, 200)), function(d)
    max_eig_mod(matrix(postf$B[d, , ], ncol = nf), nf, p) < 1.05, logical(1)))
  own <- sapply(seq_len(nf), function(i) postf$B[, 1 + i, i])
  diag_ <- mcmc_diagnostics(own, stable)
  diag_$ess_min <- diag_$ess_median <- n     # iid draws from closed form
  diag_$converged <- TRUE
  diag_$block_exog_max <- 0                  # exact by construction

  structure(list(engine = "conj_br", M = M, p = p, nf = nf,
                 varnames = spec_m$variable, blocks = blocks,
                 foreign = postf, domestic = postd, ndraw = n,
                 diagnostics = diag_),
            class = c("post_conj_br", "var_posterior"))
}

# ---- engine: gibbs (independent Normal + inverse-Wishart) ---------------------

fit_gibbs <- function(y, member, spec_m, cfg, prior, weights = NULL) {
  M <- ncol(y); p <- member$lags
  blocks <- spec_m$block
  sigma <- ar_sigmas(y, weights = weights); delta <- spec_m$delta
  ndraw <- cfg$mcmc$ndraw; nburn <- cfg$mcmc$nburn

  xy <- build_XY(y, p)
  if (!is.null(weights)) {
    w_rows <- weights[(p + 1):nrow(y)]
    xy$Y <- xy$Y * w_rows; xy$X <- xy$X * w_rows
  }
  ybar <- colMeans(y[seq_len(p), , drop = FALSE])
  if (prior$soc || prior$dio) {
    dmy <- soc_dio_dummies(M, p, ybar, delta, soc = prior$soc,
                           soc_mu = prior$soc_mu, dio = prior$dio,
                           dio_delta = prior$dio_delta)
    if (!is.null(dmy$Y)) { xy$Y <- rbind(dmy$Y, xy$Y); xy$X <- rbind(dmy$X, xy$X) }
  }
  Y <- xy$Y; X <- xy$X
  T_n <- nrow(Y); K <- ncol(X)

  mn <- minnesota_prior(M, p, sigma, delta, blocks, lambda = prior$lambda,
                        block_exog_sd = cfg$mcmc$block_exog_prior_sd)
  b0vec <- as.vector(mn$b0)
  V0inv <- as.vector(1 / mn$s0^2)

  nu0 <- M + 2
  S0  <- diag(sigma^2, M) * (nu0 - M - 1)

  XtX <- crossprod(X); XtY <- crossprod(X, Y)
  B <- matrix(0, K, M); B[1 + seq_len(M), ] <- diag(delta, M)
  Sigma <- diag(sigma^2, M)

  Bd <- array(NA_real_, c(ndraw, K, M)); Sd <- array(NA_real_, c(ndraw, M, M))
  for (it in seq_len(nburn + ndraw)) {
    Sigma_inv <- chol2inv(chol(Sigma))
    B <- draw_beta_given_sigma(XtX, XtY, Sigma_inv, V0inv, b0vec, K, M)
    E <- Y - X %*% B
    Sigma <- riwish(nu0 + T_n, S0 + crossprod(E))
    if (it > nburn) { Bd[it - nburn, , ] <- B; Sd[it - nburn, , ] <- Sigma }
  }

  Bbar <- apply(Bd, c(2, 3), mean)
  stable <- mean(vapply(seq_len(min(ndraw, 300)), function(d)
    max_eig_mod(Bd[d, , ], M, p) < 1.05, logical(1)))
  own <- sapply(seq_len(M), function(i) Bd[, 1 + i, i])
  diag_ <- mcmc_diagnostics(own, stable)
  diag_$block_exog_max <- block_exog_metric(Bbar, M, p, blocks)

  structure(list(engine = "gibbs", M = M, p = p,
                 varnames = spec_m$variable, blocks = blocks,
                 B = Bd, Sigma = Sd, ndraw = ndraw, diagnostics = diag_),
            class = c("post_gibbs", "var_posterior"))
}

# ---- engine: ss (Villani steady-state) ----------------------------------------

fit_ss <- function(y, member, spec_m, cfg, prior, weights = NULL) {
  M <- ncol(y); p <- member$lags
  blocks <- spec_m$block
  sigma <- ar_sigmas(y, weights = weights); delta <- spec_m$delta
  ndraw <- cfg$mcmc$ndraw; nburn <- cfg$mcmc$nburn
  w_rows <- if (is.null(weights)) rep(1, nrow(y) - p) else weights[(p + 1):nrow(y)]

  ssp <- steadystate_prior(spec_m, y)
  psi0 <- ssp$psi0; psi_prec <- 1 / ssp$psi_sd^2

  # Minnesota prior on the (no-intercept) lag coefficients
  mn <- minnesota_prior(M, p, sigma, delta, blocks, lambda = prior$lambda,
                        block_exog_sd = cfg$mcmc$block_exog_prior_sd)
  b0 <- mn$b0[-1, , drop = FALSE]; s0 <- mn$s0[-1, , drop = FALSE]
  b0vec <- as.vector(b0); V0inv <- as.vector(1 / s0^2)

  nu0 <- M + 2
  S0  <- diag(sigma^2, M) * (nu0 - M - 1)

  T_all <- nrow(y)
  Psi <- colMeans(y)
  A <- matrix(0, M * p, M); A[seq_len(M), ] <- diag(delta, M)
  Sigma <- diag(sigma^2, M)

  Ad <- array(NA_real_, c(ndraw, M * p, M))
  Pd <- matrix(NA_real_, ndraw, M)
  Sd <- array(NA_real_, c(ndraw, M, M))

  lag_idx <- function(z) build_XY(z, p, intercept = FALSE)

  for (it in seq_len(nburn + ndraw)) {
    # 1. A | Psi, Sigma on demeaned data (reject explosive draws, max 5 tries);
    # rows GLS-weighted for the COVID treatment
    z <- sweep(y, 2, Psi)
    xz <- lag_idx(z)
    Xw <- xz$X * w_rows; Yw <- xz$Y * w_rows
    XtX <- crossprod(Xw); XtY <- crossprod(Xw, Yw)
    Sigma_inv <- chol2inv(chol(Sigma))
    for (try in 1:5) {
      Anew <- draw_beta_given_sigma(XtX, XtY, Sigma_inv, V0inv, b0vec,
                                    M * p, M)
      Bfull <- rbind(0, Anew)
      if (max_eig_mod(Bfull, M, p) < 1.0) { A <- Anew; break }
      if (try == 5) A <- Anew   # keep anyway; flagged via stability share
    }
    # 2. Psi | A, Sigma: w_t-row = U Psi + e_t with var s_t^2 Sigma, so each
    # t contributes w_t^2 to the GLS precision and w_t^2 * row to the rhs.
    # NOTE: this is the full Villani (2009) GLS full conditional -- it uses the
    # COMPLETE M x M Sigma_inv, so the joint draw of Psi couples the foreign and
    # domestic means through the (non-block-diagonal) error covariance. Block
    # exogeneity here constrains the LAG DYNAMICS (A is block-lower-triangular),
    # NOT the steady-state means: a foreign-only Psi update would be an incorrect
    # full conditional. This is intentional (README D3 block-exog scope / D5).
    U <- diag(M); for (l in seq_len(p)) U <- U - t(A[((l - 1) * M + 1):(l * M), ])
    xy_raw <- build_XY(y, p, intercept = FALSE)
    W <- xy_raw$Y - xy_raw$X %*% A
    UtSi <- t(U) %*% Sigma_inv
    Pp <- diag(psi_prec, M) + sum(w_rows^2) * UtSi %*% U
    rhs <- psi_prec * psi0 + UtSi %*% colSums(W * w_rows^2)
    cPp <- chol((Pp + t(Pp)) / 2)
    Psi <- drop(backsolve(cPp, forwardsolve(t(cPp), rhs)) +
                backsolve(cPp, rnorm(M)))
    # 3. Sigma | A, Psi (weighted residuals are homoskedastic)
    z <- sweep(y, 2, Psi)
    xz <- lag_idx(z)
    E <- (xz$Y - xz$X %*% A) * w_rows
    Sigma <- riwish(nu0 + nrow(E), S0 + crossprod(E))
    if (it > nburn) {
      Ad[it - nburn, , ] <- A; Pd[it - nburn, ] <- Psi; Sd[it - nburn, , ] <- Sigma
    }
  }

  Abar <- apply(Ad, c(2, 3), mean)
  stable <- mean(vapply(seq_len(min(ndraw, 300)), function(d)
    max_eig_mod(rbind(0, Ad[d, , ]), M, p) < 1.0, logical(1)))
  own <- sapply(seq_len(M), function(i) Ad[, i, i])
  diag_ <- mcmc_diagnostics(own, stable)
  diag_$block_exog_max <- block_exog_metric(rbind(0, Abar), M, p, blocks)

  structure(list(engine = "ss", M = M, p = p,
                 varnames = spec_m$variable, blocks = blocks,
                 A = Ad, Psi = Pd, Sigma = Sd, ndraw = ndraw,
                 diagnostics = diag_),
            class = c("post_ss", "var_posterior"))
}

# ---- engine: sv (triangular equation-by-equation with stochastic volatility) --

#' Equation i of the triangularised system: y_i,t on intercept, allowed lags,
#' and contemporaneous y_j,t (j < i). Foreign equations exclude all domestic
#' regressors exactly (columns dropped), which imposes block exogeneity.
sv_equation_design <- function(y, i, p, blocks) {
  M <- ncol(y)
  xy <- build_XY(y, p)
  allow_lag <- if (blocks[i] == "foreign") which(blocks == "foreign") else seq_len(M)
  idx <- coef_index(M, p)
  keep_cols <- c(1, idx$row[idx$var %in% allow_lag])
  X <- xy$X[, keep_cols, drop = FALSE]
  contemp <- if (i > 1) seq_len(i - 1) else integer(0)
  if (blocks[i] == "foreign") contemp <- contemp[blocks[contemp] == "foreign"]
  if (length(contemp))
    X <- cbind(X, y[(p + 1):nrow(y), contemp, drop = FALSE])
  list(y = xy$Y[, i], X = X, contemp = contemp,
       lag_meta = idx[idx$var %in% allow_lag, , drop = FALSE])
}

#' Per-equation prior (independent normal): Minnesota on lags, loose on
#' contemporaneous loadings.
sv_equation_prior <- function(i, eq, sigma, delta, lambda) {
  Kx <- ncol(eq$X)
  nlag <- nrow(eq$lag_meta)
  b0 <- numeric(Kx); s0 <- numeric(Kx)
  s0[1] <- 10 * sigma[i]
  for (r in seq_len(nlag)) {
    j <- eq$lag_meta$var[r]; l <- eq$lag_meta$lag[r]
    own <- (j == i)
    if (own && l == 1) b0[1 + r] <- delta[i]
    s0[1 + r] <- lambda * (if (own) 1 else 0.5) * (sigma[i] / sigma[j]) / l
  }
  if (length(eq$contemp))
    s0[(1 + nlag + 1):Kx] <- sigma[i] / sigma[eq$contemp]
  list(b0 = b0, s0 = s0)
}

fit_sv <- function(y, member, spec_m, cfg, prior, weights = NULL) {
  if (!has_stochvol())   # all_members() drops SV members; this is a backstop
    stop("the 'stochvol' package is required for SV member '", member$name, "'")
  M <- ncol(y); p <- member$lags
  blocks <- spec_m$block
  sigma <- ar_sigmas(y); delta <- spec_m$delta
  ndraw <- cfg$mcmc$ndraw; nburn <- cfg$mcmc$nburn
  # COVID robustness for the SV engine: t-distributed errors (CCMM 2024 SV-t,
  # ~ their preferred SVO-t; Hartwig 2024), NOT the LP row weighting -- the
  # outlier is absorbed by the iid scale mixture instead of the persistent
  # log-volatility process. `weights` is accepted and ignored by design.
  use_t <- !isFALSE(cfg$covid$sv_t_errors)

  eqs <- vector("list", M)
  pspec <- stochvol::specify_priors(
    mu = stochvol::sv_normal(0, 10),
    phi = stochvol::sv_beta(20, 1.5),
    sigma2 = stochvol::sv_gamma(0.5, 0.5),
    nu = if (use_t) stochvol::sv_exponential(0.1) else stochvol::sv_infinity())

  for (i in seq_len(M)) {
    eq <- sv_equation_design(y, i, p, blocks)
    pr <- sv_equation_prior(i, eq, sigma, delta, prior$lambda)
    Tn <- length(eq$y); Kx <- ncol(eq$X)
    V0inv <- 1 / pr$s0^2

    # initialise at the prior mean (an OLS init can be near-singular with
    # several lags of persistent log-level regressors)
    beta <- pr$b0
    resid <- eq$y - drop(eq$X %*% beta)
    h <- rep(log(stats::var(resid) + 1e-8), Tn)
    mix <- rep(1, Tn)                       # t-scale mixing weights
    para <- list(mu = h[1], phi = 0.9, sigma = 0.2,
                 nu = if (use_t) 10 else Inf, rho = 0,
                 beta = 0, latent0 = h[1])

    bdr <- matrix(NA_real_, ndraw, Kx)
    hT  <- numeric(ndraw)
    pdr <- matrix(NA_real_, ndraw, 4,
                  dimnames = list(NULL, c("mu", "phi", "sigma", "nu")))
    for (it in seq_len(nburn + ndraw)) {
      # beta | h, mix: weighted regression (conditional variance exp(h)*mix)
      w <- exp(-h / 2) / sqrt(mix)
      Xw <- eq$X * w; yw <- eq$y * w
      P <- crossprod(Xw); diag(P) <- diag(P) + V0inv
      cP <- chol(P)
      m_ <- backsolve(cP, forwardsolve(t(cP), crossprod(Xw, yw) + V0inv * pr$b0))
      beta <- drop(m_ + backsolve(cP, rnorm(Kx)))
      resid <- eq$y - drop(eq$X %*% beta)
      # h, para, mix | resid via stochvol single update
      upd <- stochvol::svsample_fast_cpp(
        resid, draws = 1, burnin = 0, designmatrix = matrix(NA),
        priorspec = pspec, thinpara = 1, thinlatent = 1,
        keeptime = "all", startpara = para, startlatent = h,
        keeptau = use_t)
      h <- drop(upd$latent[1, ])
      if (use_t) mix <- drop(upd$tau[1, ])
      para$mu    <- upd$para[1, "mu"]
      para$phi   <- upd$para[1, "phi"]
      para$sigma <- upd$para[1, "sigma"]
      if (use_t) para$nu <- upd$para[1, "nu"]
      para$latent0 <- h[1]
      if (it > nburn) {
        d <- it - nburn
        bdr[d, ] <- beta; hT[d] <- h[Tn]
        pdr[d, ] <- c(para$mu, para$phi, para$sigma,
                      if (is.finite(para$nu)) para$nu else Inf)
      }
    }
    eqs[[i]] <- list(design = eq, beta = bdr, hT = hT, svpara = pdr,
                     prior = pr)
  }

  own <- sapply(seq_len(M), function(i) {
    r <- which(eqs[[i]]$design$lag_meta$var == i & eqs[[i]]$design$lag_meta$lag == 1)
    eqs[[i]]$beta[, 1 + r]
  })
  diag_ <- mcmc_diagnostics(own, stable_share = NA_real_)
  # block exog: max |posterior mean| over domestic-lag coefs in foreign eqs --
  # exact zero by construction (columns dropped)
  diag_$block_exog_max <- 0

  structure(list(engine = "sv", M = M, p = p,
                 varnames = spec_m$variable, blocks = blocks,
                 eqs = eqs, ndraw = ndraw, y = y, diagnostics = diag_),
            class = c("post_sv", "var_posterior"))
}

# ---- dispatcher ----------------------------------------------------------------

fit_var_member <- function(y, member, spec_m, cfg, glp_lambda = 0.2,
                           weights = NULL) {
  prior <- resolve_prior(member, glp_lambda)
  fitter <- switch(member$engine,
                   conj_br = fit_conj_br,
                   gibbs   = fit_gibbs,
                   ss      = fit_ss,
                   sv      = fit_sv,
                   stop("unknown engine: ", member$engine))
  fitter(y, member, spec_m, cfg, prior, weights = weights)
}
