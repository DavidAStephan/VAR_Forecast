source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()

cfg0 <- load_config(testthat::test_path("..", "..", "config", "config.yml"))

test_that("AR draws are stationarity-truncated (audit finding 1 regression)", {
  set.seed(4)
  # explosive-ish data: AR(1) coef 1.01
  z <- as.numeric(stats::filter(rnorm(120), 1.01, method = "recursive"))
  y <- matrix(z, ncol = 1, dimnames = list(NULL, "x"))
  post <- fit_ar(y, cfg0, p = 4)
  paths <- simulate_paths(post, y, 12, 400)
  # every retained coefficient draw must imply a stable companion matrix:
  # the 12-step path central mass cannot explode
  q99 <- apply(abs(paths[, , 1]), 2, quantile, probs = 0.99)
  expect_lt(max(q99) / max(abs(z)), 8)
})

test_that("covid_s_path and covid_s_future implement the LP decay", {
  dates <- seq(as.Date("2019-03-01"), by = "quarter", length.out = 10)
  cq <- as.Date(c("2020-03-01", "2020-06-01"))
  s <- covid_s_path(dates, cq, scales = c(10, 5), rho = 0.5)
  expect_equal(s[dates < cq[1]], rep(1, sum(dates < cq[1])))
  expect_equal(s[dates == cq[1]], 10)
  expect_equal(s[dates == cq[2]], 5)
  # first quarter after: 1 + (5-1)*0.5 = 3
  expect_equal(s[dates == as.Date("2020-09-01")], 3)
  expect_equal(s[dates == as.Date("2020-12-01")], 2)
  sf <- covid_s_future(dates, 3, cq, scales = c(10, 5), rho = 0.5)
  # 6 quarters already elapsed past 2020Q2 within dates (2020Q3..2021Q4)
  j0 <- sum(dates > max(cq))
  expect_equal(sf, 1 + 4 * 0.5^(j0 + 1:3))
})

test_that("LP marginal likelihood detects the synthetic COVID outlier", {
  cfg <- cfg0
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  spec_s <- vars_for_set(spec, "small")
  y <- as.matrix(td[1:126, spec_s$variable])
  dates <- td$date[1:126]
  est <- estimate_covid_scales(y, dates, 4, ar_sigmas(y), spec_s$delta, cfg)
  expect_length(est$scales, 3)
  # the outlier quarters demand scales well above 1
  expect_gte(max(est$scales), 8)
  # and the ML must prefer them to no scaling
  rows <- dates[5:length(dates)]
  ml_scaled <- .lp_logml(y, 4, ar_sigmas(y), spec_s$delta, 0.2,
                         covid_s_path(rows, est$cq, est$scales, est$rho))
  ml_unit <- .lp_logml(y, 4, ar_sigmas(y), spec_s$delta, 0.2,
                       rep(1, length(rows)))
  expect_gt(ml_scaled, ml_unit)
})

test_that("LP weighting protects coefficient estimates (parameter recovery)", {
  cfg <- cfg0
  cfg$mcmc$ndraw <- 300; cfg$mcmc$nburn <- 80
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  spec_s <- vars_for_set(spec, "small")
  member <- list(name = "t", kind = "var", engine = "gibbs", set = "small",
                 lags = 2, prior = list(lambda = 0.2, soc = FALSE, dio = FALSE))
  y_pre <- as.matrix(td[1:118, spec_s$variable])     # ends before the outlier
  y_post <- as.matrix(td[1:124, spec_s$variable])    # includes outlier+rebound
  dates_post <- td$date[1:124]
  fitB <- function(y, w) {
    set.seed(99)
    post <- fit_var_member(y, member, spec_s, cfg, weights = w)
    apply(post$B, c(2, 3), mean)
  }
  B_pre <- fitB(y_pre, NULL)                          # clean reference
  B_raw <- fitB(y_post, NULL)                         # contaminated
  ct <- covid_treatment(y_post, dates_post, 2, ar_sigmas(y_post),
                        spec_s$delta, cfg, 12)
  B_lp <- fitB(y_post, ct$weights)                    # LP-weighted
  d_raw <- sqrt(mean((B_raw - B_pre)^2))
  d_lp <- sqrt(mean((B_lp - B_pre)^2))
  # the treated estimates must be substantially closer to the clean ones
  expect_lt(d_lp, 0.7 * d_raw)
})

test_that("dummy treatment is equivalent to dropping the COVID rows (coefficients)", {
  set.seed(11)
  Tn <- 90; M <- 2; p <- 1
  y <- matrix(rnorm(Tn * M), Tn, M)
  y[, 1] <- as.numeric(stats::filter(y[, 1], 0.6, method = "recursive"))
  y[60, ] <- y[60, ] + c(12, -9)               # huge outlier at row 60
  spec_m <- data.frame(variable = c("a", "b"), block = c("foreign", "domestic"),
                       delta = c(0, 0))
  cfg <- cfg0; cfg$mcmc$forecast_draws <- 200
  member <- list(name = "t", kind = "var", engine = "conj_br", set = "x",
                 lags = p, prior = list(lambda = 0.2, soc = FALSE, dio = FALSE))
  # Downweighting removes observation 60 from the LIKELIHOOD; y_60 stays in
  # the regressors of row 61 by design (the jump-off/conditioning channel,
  # CCMM 2024). The testable limit property: once the weight is tiny, the
  # posterior is insensitive to how tiny (the dropping limit is reached), and
  # it differs from the unweighted fit.
  fitB <- function(yy, w) {
    set.seed(5)
    apply(fit_var_member(yy, member, spec_m, cfg,
                         weights = w)$domestic$B, c(2, 3), mean)
  }
  w3 <- rep(1, Tn); w3[60] <- 1e-3
  w6 <- rep(1, Tn); w6[60] <- 1e-6
  B3 <- fitB(y, w3); B6 <- fitB(y, w6); B1 <- fitB(y, NULL)
  expect_equal(B3, B6, tolerance = 1e-3)        # limit reached
  expect_gt(sqrt(mean((B1 - B3)^2)), 0.01)      # and it matters vs untreated
})

test_that("SV engine t-errors: nu is estimated finite and draws use t tails", {
  cfg <- cfg0
  cfg$mcmc$ndraw <- 150; cfg$mcmc$nburn <- 50
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  spec_s <- vars_for_set(spec, "small")
  y <- as.matrix(td[1:126, spec_s$variable])   # includes the outlier
  member <- list(name = "t", kind = "var", engine = "sv", set = "small",
                 lags = 2, prior = list(lambda = 0.2))
  set.seed(2)
  post <- fit_var_member(y, member, spec_s, cfg)
  nus <- sapply(post$eqs, function(e) median(e$svpara[, "nu"]))
  expect_true(all(is.finite(nus)))
  expect_true(all(nus > 2))
})

test_that("no-look-ahead holds with the COVID treatment active", {
  cfg <- cfg0
  cfg$mcmc$ndraw <- 80; cfg$mcmc$nburn <- 20
  cfg$mcmc$forecast_draws <- 80; cfg$mcmc$store_draws <- 40
  cfg$evaluation$max_origins <- 4
  cfg$glp$enabled <- FALSE
  spec <- build_transform_spec(cfg)
  raw <- generate_synthetic_data(cfg, spec)
  td <- transform_data(raw, spec)
  expect_true(test_no_lookahead(td, spec, cfg))
})
