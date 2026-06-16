# covid.R -- COVID-period volatility treatment (Lenza-Primiceri 2022) ----------
#
# Model: y_t = c + B(L)y_{t-1} + s_t * eps_t, eps_t ~ N(0, Sigma), with s_t = 1
# outside the pandemic, a free scale factor for each configured COVID quarter,
# and geometric decay of the EXCESS volatility afterwards:
#   s_{T0+j} = 1 + (s_last - 1) * rho^j.
# Dividing the observation equation at time t by s_t restores a homoskedastic
# VAR in transformed data (rows of Y and X scaled by 1/s_t), so every
# constant-volatility engine -- including the Kronecker conjugate -- applies
# unchanged (LP 2022, JAE). The scales are estimated by maximizing the
# closed-form conjugate marginal likelihood of the rescaled data INCLUDING the
# Jacobian -M*sum(log s_t) (LP eq. A.1), via coordinate descent on a grid --
# the GLP-style treatment of (s, rho) as hyperparameters.
#
# treatment = "dummy" instead sets s_t -> large at the COVID quarters
# (equivalent to dummying-out/dropping those rows from the likelihood, cf.
# Schorfheide-Song 2024; Cascaldi-Garcia's Pandemic Priors with phi -> 0) and
# implies s = 1 for the forecast path: the documented cost is predictive bands
# stuck at pre-COVID widths (CCMM 2024).
#
# The SV engine is NOT weighted: it gets t-distributed errors instead (see
# engines.R / README.md D17), per CCMM (SV-t ~ SVO-t) and Hartwig (2024).

#' Quarter index (year*4 + quarter): all date comparisons in this module are
#' done at quarter granularity, because the real panel stamps quarters at
#' their start (2020-01-01) while the synthetic panel stamps them mid-quarter
#' (2020-03-01) -- exact-date matching would silently disable the treatment.
.qidx <- function(d) {
  d <- as.Date(d)
  as.integer(format(d, "%Y")) * 4L + (as.integer(format(d, "%m")) - 1L) %/% 3L
}

#' Which configured COVID quarters fall inside the estimation window.
covid_quarters_in <- function(dates, cfg) {
  if (is.null(cfg$covid) || identical(cfg$covid$treatment, "none"))
    return(as.Date(character(0)))
  q <- as.Date(unlist(cfg$covid$quarters))
  q[.qidx(q) <= max(.qidx(dates))]
}

#' Volatility path s_t over the sample dates given per-quarter scales and the
#' decay rho applied to quarters after the last scaled quarter.
covid_s_path <- function(dates, cq, scales, rho) {
  s <- rep(1, length(dates))
  if (!length(cq)) return(s)
  di <- .qidx(dates); qi <- .qidx(cq)
  for (i in seq_along(qi)) s[di == qi[i]] <- scales[i]
  after <- which(di > max(qi))
  if (length(after) && scales[length(scales)] > 1) {
    j <- di[after] - max(qi)
    s[after] <- 1 + (scales[length(scales)] - 1) * rho^j
  }
  pmax(s, 1)
}

#' Future volatility path for horizons 1..H beyond the origin (used to scale
#' predictive shocks at COVID-era origins -- the LP forecasting payoff).
covid_s_future <- function(dates, H, cq, scales, rho) {
  if (!length(cq) || scales[length(scales)] <= 1) return(rep(1, H))
  j0 <- max(.qidx(dates)) - max(.qidx(cq))   # quarters elapsed past T0
  j0 <- max(j0, 0)
  1 + (scales[length(scales)] - 1) * rho^(j0 + seq_len(H))
}

#' Log marginal likelihood of the LP-scaled conjugate model: closed-form
#' conjugate ML on the rescaled rows plus the Jacobian -M * sum(log s_t).
.lp_logml <- function(y, p, sigma, delta, lambda, s_rows) {
  xy <- build_XY(y, p)
  w <- 1 / s_rows
  Yw <- xy$Y * w
  Xw <- xy$X * w
  pr <- conjugate_prior(ncol(y), p, sigma, delta, lambda = lambda)
  log_marginal_likelihood(Yw, Xw, pr) - ncol(y) * sum(log(s_rows))
}

#' Estimate the COVID scale factors and decay by coordinate descent on grids,
#' maximizing the LP marginal likelihood. Uses only data in `y`/`dates`
#' (i.e. up to the forecast origin), so the recursive no-look-ahead property
#' holds. Returns list(cq, scales, rho); scales is empty when no COVID quarter
#' is in sample or treatment is "none".
estimate_covid_scales <- function(y, dates, p, sigma, delta, cfg, lambda = 0.2) {
  cq <- covid_quarters_in(dates, cfg)
  if (!length(cq)) return(list(cq = cq, scales = numeric(0), rho = 0))
  if (identical(cfg$covid$treatment, "dummy")) {
    return(list(cq = cq, scales = rep(1e3, length(cq)), rho = 0))
  }
  sgrid <- as.numeric(unlist(cfg$covid$scale_grid))
  rgrid <- as.numeric(unlist(cfg$covid$decay_grid))
  scales <- rep(sgrid[ceiling(length(sgrid) / 2)], length(cq))
  rho <- rgrid[1]
  rows <- dates[(p + 1):length(dates)]     # dates aligned with regression rows
  eval_ml <- function(scales, rho) {
    s_rows <- covid_s_path(rows, cq, scales, rho)
    .lp_logml(y, p, sigma, delta, lambda, s_rows)
  }
  for (sweep in 1:2) {
    for (i in seq_along(cq)) {
      mls <- vapply(sgrid, function(s) {
        sc <- scales; sc[i] <- s; eval_ml(sc, rho)
      }, numeric(1))
      scales[i] <- sgrid[which.max(mls)]
    }
    mls <- vapply(rgrid, function(r) eval_ml(scales, r), numeric(1))
    rho <- rgrid[which.max(mls)]
  }
  log_debug("covid scales at {max(dates)}: {paste(scales, collapse='/')} rho={rho}")
  list(cq = cq, scales = scales, rho = rho)
}

#' One-call helper for the evaluation harness: returns the row weights
#' (length T, multiply rows of y's regression representation; NULL when
#' inactive) and the future shock-scale path (length H).
covid_treatment <- function(y, dates, p, sigma, delta, cfg, H, lambda = 0.2) {
  none <- list(weights = NULL, s_future = NULL, scales = numeric(0), rho = 0)
  if (is.null(cfg$covid) || identical(cfg$covid$treatment, "none")) return(none)
  est <- estimate_covid_scales(y, dates, p, sigma, delta, cfg, lambda)
  if (!length(est$cq)) return(none)
  s_all <- covid_s_path(dates, est$cq, est$scales, est$rho)
  s_fut <- if (identical(cfg$covid$treatment, "dummy")) rep(1, H)
           else covid_s_future(dates, H, est$cq, est$scales, est$rho)
  list(weights = 1 / s_all, s_future = s_fut,
       scales = est$scales, rho = est$rho)
}
