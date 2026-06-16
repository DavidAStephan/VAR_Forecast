# benchmarks.R -- univariate benchmark members: rw, ar4, ucsv, ucmean ----------
#
# Benchmarks implement the same interface as VAR members: fit_benchmark()
# returns an object with a simulate_paths() method producing [ndraw, h, nvar]
# arrays over the *target* variables. They are pool members and the bar every
# VAR must clear.

#' Random walk on the modelled (transformed) series with Gaussian increments.
#' Optional weights downweight COVID-quarter increments in the sd estimate.
fit_rw <- function(y, cfg, weights = NULL) {
  w <- if (is.null(weights)) rep(1, nrow(y) - 1) else pmin(weights[-1], weights[-nrow(y)])
  sdd <- apply(diff(y) * w, 2, sd)
  structure(list(engine = "rw", varnames = colnames(y),
                 yT = y[nrow(y), ], sdd = sdd,
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_rw", "var_posterior"))
}

simulate_paths.post_rw <- function(post, y, h, ndraw, condition = NULL,
                                   shock_scale = NULL) {
  nv <- length(post$yT)
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    inc <- matrix(rnorm(ndraw * h, 0, post$sdd[j]), ndraw, h) *
           matrix(ss_, ndraw, h, byrow = TRUE)
    paths[, , j] <- post$yT[j] + t(apply(inc, 1, cumsum))
  }
  paths
}

#' Bayesian AR(p) per variable: conjugate normal-inverse-gamma with
#' Minnesota-style lag shrinkage (sd = 0.5/lag); iterated density forecasts.
#' A flat ridge lets near-unit oscillatory coefficient draws resonate off
#' COVID-sized outliers in the initial conditions.
fit_ar <- function(y, cfg, p = 4, weights = NULL) {
  fits <- lapply(seq_len(ncol(y)), function(j) {
    z <- y[, j]
    xy <- build_XY(matrix(z, ncol = 1), p)
    X <- xy$X; Y <- drop(xy$Y)
    if (!is.null(weights)) {
      w <- weights[(p + 1):length(z)]
      X <- X * w; Y <- Y * w
    }
    K <- ncol(X)
    V0inv <- diag(c(1e-4, (seq_len(p) / 0.5)^2))  # loose intercept, sd 0.5/l on lag l
    P <- crossprod(X) + V0inv
    cP <- chol(P)
    bhat <- backsolve(cP, forwardsolve(t(cP), crossprod(X, Y)))
    resid <- Y - drop(X %*% bhat)
    s2 <- sum(resid^2) / (length(Y) - K)
    list(bhat = bhat, cP = cP, s2 = s2, df = length(Y) - K, p = p)
  })
  structure(list(engine = "ar4", varnames = colnames(y), fits = fits, p = p,
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ar", "var_posterior"))
}

simulate_paths.post_ar <- function(post, y, h, ndraw, condition = NULL,
                                   shock_scale = NULL) {
  nv <- length(post$fits); p <- post$p
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    f <- post$fits[[j]]
    z <- y[, j]
    for (d in seq_len(ndraw)) {
      sig2 <- f$s2 * f$df / rchisq(1, f$df)
      # stationarity-truncated posterior: an explosive AR draw compounds over
      # 12 iterated steps (seen at origins straddling the COVID outlier).
      # Stationarity requires ALL roots of 1 - b1 z - ... - bp z^p OUTSIDE the
      # unit circle, i.e. min |root| > 1 (max companion eigenvalue < 1).
      for (try in 1:20) {
        beta <- f$bhat + sqrt(sig2) * backsolve(f$cP, rnorm(p + 1))
        rmin <- min(Mod(polyroot(c(1, -beta[-1]))))
        if (rmin > 1.0) break
        if (try == 20) {
          # deflate lag-l coefficient by (0.95*rmin)^l: scales the max
          # companion eigenvalue 1/rmin down to ~0.95/1
          beta[-1] <- beta[-1] * (0.95 * rmin)^seq_len(p)
        }
      }
      st <- z[length(z) - seq_len(p) + 1]      # most recent first
      for (s in seq_len(h)) {
        ynew <- sum(c(1, st) * beta) + ss_[s] * sqrt(sig2) * rnorm(1)
        paths[d, s, j] <- ynew
        st <- c(ynew, st[-p])
      }
    }
  }
  paths
}

# ---- UCSV (Stock-Watson unobserved components with twin SV) --------------------

#' FFBS for the local level model with time-varying variances:
#' y_t = tau_t + e_t, e_t ~ N(0, s2e_t); tau_t = tau_{t-1} + u_t, u_t ~ N(0, s2u_t).
ffbs_local_level <- function(y, s2e, s2u, tau0 = y[1], P0 = 10 * var(y)) {
  Tn <- length(y)
  af <- numeric(Tn); Pf <- numeric(Tn)
  a <- tau0; P <- P0
  for (t in seq_len(Tn)) {
    P <- P + s2u[t]
    Kk <- P / (P + s2e[t])
    a <- a + Kk * (y[t] - a)
    P <- (1 - Kk) * P
    af[t] <- a; Pf[t] <- P
  }
  tau <- numeric(Tn)
  tau[Tn] <- rnorm(1, af[Tn], sqrt(Pf[Tn]))
  for (t in (Tn - 1):1) {
    Pp <- Pf[t] + s2u[t + 1]
    J <- Pf[t] / Pp
    m <- af[t] + J * (tau[t + 1] - af[t])
    V <- Pf[t] * (1 - J)
    tau[t] <- rnorm(1, m, sqrt(max(V, 1e-12)))
  }
  tau
}

#' UC-SV per variable via Gibbs: FFBS for the trend; stochastic volatility with
#' t-distributed errors on the transitory component (outlier-robust,
#' Lenza-Primiceri spirit); CONSTANT trend-shock variance with a conjugate
#' inverse-gamma update. The canonical twin-SV UCSV mixes pathologically when a
#' COVID-sized outlier sits at the sample endpoint (trend-vs-noise attribution
#' is bimodal and the trend-vol process piles up near zero); a direct IG draw
#' for the trend variance removes that dimension entirely. See README.md D9.
fit_ucsv <- function(y, cfg) {
  ndraw <- cfg$mcmc$ndraw; nburn <- cfg$mcmc$nburn
  thin <- if (is.null(cfg$mcmc$bench_thin)) 1 else cfg$mcmc$bench_thin
  # tight vol-of-vol prior (mean 0.04): Stock-Watson fix gamma ~ 0.2; a loose
  # prior lets an endpoint outlier inflate sigma_eta and the 12-step variance
  # forecast explodes through exp(h/2) compounding
  pspec_e <- stochvol::specify_priors(
    mu = stochvol::sv_normal(0, 10),
    phi = stochvol::sv_beta(20, 1.5),
    sigma2 = stochvol::sv_gamma(2, 50),
    nu = stochvol::sv_exponential(0.1))
  run_chain <- function(z, thin_j) {
    Tn <- length(z)
    tau <- stats::filter(z, rep(1 / 8, 8), sides = 1)
    tau[is.na(tau)] <- z[is.na(tau)]
    tau <- as.numeric(tau)
    he <- rep(log(var(z) / 2 + 1e-8), Tn)
    mix_e <- rep(1, Tn)                       # t-mixing weights, transitory eq
    # informative IG prior keeps the trend-shock variance off zero:
    # mean ~ 2% of series variance, ~10 prior "observations"
    a0 <- 5; b0 <- (a0 - 1) * 0.02 * var(z)
    s2u <- b0 / (a0 - 1)
    pe <- list(mu = he[1], phi = 0.9, sigma = 0.2, nu = 10, rho = 0, beta = 0,
               latent0 = he[1])
    tauT <- numeric(ndraw); heT <- numeric(ndraw); s2uT <- numeric(ndraw)
    pe_d <- matrix(NA_real_, ndraw, 4)
    for (it in seq_len(nburn + ndraw * thin_j)) {
      tau <- ffbs_local_level(z, exp(he) * mix_e, rep(s2u, Tn))
      eres <- z - tau
      ures <- diff(tau)
      # conjugate IG update for the constant trend-shock variance
      s2u <- 1 / rgamma(1, a0 + length(ures) / 2, b0 + sum(ures^2) / 2)
      ue <- stochvol::svsample_fast_cpp(eres, draws = 1, burnin = 0,
        priorspec = pspec_e, startpara = pe, startlatent = he, keeptau = TRUE)
      he <- drop(ue$latent[1, ])
      mix_e <- drop(ue$tau[1, ])
      pe$mu <- ue$para[1, "mu"]; pe$phi <- ue$para[1, "phi"]
      pe$sigma <- ue$para[1, "sigma"]; pe$nu <- ue$para[1, "nu"]
      pe$latent0 <- he[1]
      if (it > nburn && (it - nburn) %% thin_j == 0) {
        d <- (it - nburn) %/% thin_j
        tauT[d] <- tau[Tn]; heT[d] <- he[Tn]; s2uT[d] <- s2u
        pe_d[d, ] <- c(pe$mu, pe$phi, pe$sigma, pe$nu)
      }
    }
    list(tauT = tauT, heT = heT, s2uT = s2uT, pe = pe_d)
  }
  # Convergence gate: the trend/noise SPLIT of a near-white series is weakly
  # identified in any UC model (classic pile-up), so raw parameter ESS is the
  # wrong criterion for a FORECASTING benchmark -- the predictive sum is what
  # must be Monte-Carlo-precise. Gate: MCSE of the trend endpoint must be
  # small relative to the one-step predictive sd. ESS values are still
  # reported. README.md D9.
  tau_ess <- function(f) as.numeric(coda::effectiveSize(f$tauT))
  mcse_ok <- function(f) {
    pred_sd <- sqrt(var(f$tauT) + mean(exp(f$heT)) + mean(f$s2uT))
    mcse <- sd(f$tauT) / sqrt(max(tau_ess(f), 1))
    mcse < 0.15 * pred_sd
  }
  fits <- lapply(seq_len(ncol(y)), function(j) {
    f <- run_chain(y[, j], thin)
    if (!mcse_ok(f)) f <- run_chain(y[, j], thin * 3)   # adaptive retry
    f
  })
  ess_tau <- min(vapply(fits, tau_ess, numeric(1)))
  conv <- all(vapply(fits, mcse_ok, logical(1)))
  if (!conv)
    log_warn("ucsv: trend-endpoint MCSE exceeds 15% of predictive sd")
  structure(list(engine = "ucsv", varnames = colnames(y), fits = fits,
                 ndraw = ndraw,
                 diagnostics = list(converged = conv, ess_min = ess_tau,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ucsv", "var_posterior"))
}

simulate_paths.post_ucsv <- function(post, y, h, ndraw, condition = NULL,
                                     shock_scale = NULL) {
  # shock_scale ignored: UCSV is already outlier-robust (t errors)
  nv <- length(post$fits)
  paths <- array(NA_real_, c(ndraw, h, nv), dimnames = list(NULL, NULL, post$varnames))
  for (j in seq_len(nv)) {
    f <- post$fits[[j]]
    for (d in seq_len(ndraw)) {
      k <- ((d - 1) %% post$ndraw) + 1
      tau <- f$tauT[k]; he <- f$heT[k]; s2u <- f$s2uT[k]
      pe <- f$pe[k, ]
      nu <- f$pe[k, 4]
      for (s in seq_len(h)) {
        he <- pe[1] + pe[2] * (he - pe[1]) + pe[3] * rnorm(1)
        tau <- tau + sqrt(s2u) * rnorm(1)
        eps <- if (is.finite(nu)) rt(1, df = nu) else rnorm(1)
        paths[d, s, j] <- tau + exp(he / 2) * eps
      }
    }
  }
  paths
}

#' Unconditional mean with Gaussian predictive (expanding moments,
#' COVID-weighted when weights supplied).
fit_ucmean <- function(y, cfg, weights = NULL) {
  w <- if (is.null(weights)) rep(1, nrow(y)) else weights^2
  mu <- colSums(y * w) / sum(w)
  sdv <- sqrt(colSums(sweep(y, 2, mu)^2 * w) / sum(w))
  structure(list(engine = "ucmean", varnames = colnames(y),
                 mu = mu, sd = sdv,
                 diagnostics = list(converged = TRUE, ess_min = Inf,
                                    stable_share = 1, block_exog_max = 0)),
            class = c("post_ucmean", "var_posterior"))
}

simulate_paths.post_ucmean <- function(post, y, h, ndraw, condition = NULL,
                                       shock_scale = NULL) {
  ss_ <- if (is.null(shock_scale)) rep(1, h) else shock_scale
  nv <- length(post$mu)
  paths <- array(rnorm(ndraw * h * nv), c(ndraw, h, nv),
                 dimnames = list(NULL, NULL, post$varnames))
  scl <- matrix(ss_, ndraw, h, byrow = TRUE)
  for (j in seq_len(nv)) paths[, , j] <- post$mu[j] + paths[, , j] * post$sd[j] * scl
  paths
}

#' Dispatcher mirroring fit_var_member.
fit_benchmark <- function(y_targets, name, cfg, weights = NULL) {
  switch(name,
         rw     = fit_rw(y_targets, cfg, weights = weights),
         ar4    = fit_ar(y_targets, cfg, p = 4, weights = weights),
         ucsv   = fit_ucsv(y_targets, cfg),     # already outlier-robust (t)
         ucmean = fit_ucmean(y_targets, cfg, weights = weights),
         stop("unknown benchmark: ", name))
}

#' Hook for externally supplied forecasts (e.g. published RBA forecasts):
#' if a CSV with columns (origin_date, variable, h, q05..q95 or point, sd)
#' exists at `path`, it is read and returned for comparison in the report.
read_external_forecasts <- function(path = "data/external_forecasts.csv") {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}
