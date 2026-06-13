source_project <- function() {
  root <- testthat::test_path("..", "..")
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    sys.source(f, envir = globalenv())
}
source_project()
suppressPackageStartupMessages(library(scoringRules))

test_that("iterated forecast recursion matches a manual deterministic recursion", {
  M <- 2; p <- 2
  B <- rbind(c(0.1, -0.2),            # intercept
             c(0.5, 0.1), c(0.0, 0.4),   # lag 1
             c(-0.1, 0.0), c(0.05, 0.2)) # lag 2
  y <- matrix(c(1, 2, 0.5, -1), 2, 2, byrow = TRUE)  # y_{T-1}; y_T rows
  ystate <- y[2:1, ]                   # most recent first
  cS <- matrix(0, M, M)                # zero shocks -> deterministic
  path <- .iterate_path(B, cS, ystate, h = 3, M, p, condition = NULL,
                        varnames = c("a", "b"))
  # manual recursion
  yy <- rbind(y, matrix(NA, 3, 2))
  for (s in 3:5) {
    x <- c(1, yy[s - 1, ], yy[s - 2, ])
    yy[s, ] <- drop(crossprod(B, x))
  }
  expect_equal(path, yy[3:5, ], ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("companion matrix and stability check are correct", {
  # univariate AR(1) with coef 0.9: max eigenvalue 0.9
  B <- rbind(0, 0.9)
  expect_equal(max_eig_mod(B, M = 1, p = 1), 0.9, tolerance = 1e-12)
  # explosive
  B2 <- rbind(0, 1.1)
  expect_gt(max_eig_mod(B2, 1, 1), 1)
})

test_that("logsumexp is stable and exact", {
  x <- c(-1000, -1001)
  expect_equal(logsumexp(x), -1000 + log(1 + exp(-1)))
  expect_equal(logsumexp(c(0, 0)), log(2))
})

test_that("log score and CRPS from draws approximate analytic values", {
  set.seed(3)
  draws <- rnorm(50000)
  y0 <- 0.5
  # analytic: log N(0.5; 0,1); crps closed form for standard normal
  expect_equal(-logs_sample(y0, draws), dnorm(y0, log = TRUE), tolerance = 0.02)
  crps_an <- y0 * (2 * pnorm(y0) - 1) + 2 * dnorm(y0) - 1 / sqrt(pi)
  expect_equal(crps_sample(y0, draws), crps_an, tolerance = 0.01)
})

test_that("linear pool density and PIT are the correct mixture quantities", {
  # two Gaussian members, known mixture
  ld1 <- dnorm(0.3, 0, 1, log = TRUE); ld2 <- dnorm(0.3, 2, 1, log = TRUE)
  w <- c(0.7, 0.3)
  pool_ld <- logsumexp(log(w) + c(ld1, ld2))
  expect_equal(exp(pool_ld), 0.7 * dnorm(0.3) + 0.3 * dnorm(0.3, 2), tolerance = 1e-12)
  pit <- sum(w * c(pnorm(0.3), pnorm(0.3, 2)))
  expect_equal(pit, 0.7 * pnorm(0.3) + 0.3 * pnorm(0.3, 2))
})

test_that("optimal pool concentrates on the dominant member and stays on the simplex", {
  set.seed(5)
  n <- 120
  # member A: correct density; member B: badly biased
  ldA <- dnorm(rnorm(n), log = TRUE)
  ldB <- dnorm(rnorm(n), mean = 3, log = TRUE)
  tr <- data.frame(origin = rep(1:n, 2), h = 1,
                   member = rep(c("A", "B"), each = n),
                   logdens = c(ldA, ldB))
  w <- optimal_pool_weights(tr, c("A", "B"), t = n + 2, forgetting = 1)
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_gt(w["A"], 0.9)
})

test_that("year_ended is the rolling 4-quarter sum for dlog series", {
  x <- 1:8
  ye <- year_ended(x, "dlog")
  expect_true(all(is.na(ye[1:3])))
  expect_equal(ye[4], sum(1:4))
  expect_equal(ye[8], sum(5:8))
  expect_equal(year_ended(x, "level"), x)
})

test_that("DM test flags an obviously inferior forecast", {
  set.seed(9)
  n <- 60
  l_good <- rchisq(n, 1)
  l_bad <- l_good + 1 + 0.2 * rnorm(n)
  r <- dm_test(l_bad, l_good, h = 1)
  expect_lt(r["p"], 0.01)
  expect_gt(r["stat"], 0)
})

test_that("score_member computes q, year-ended and cumulative-level measures", {
  set.seed(1)
  H <- 4; ndraw <- 200
  # one dlog target 'g' and one level target 'r'
  td <- data.frame(date = seq(as.Date("2000-01-01"), by = "quarter", length.out = 12),
                   g = rnorm(12, 0.6, 0.5), r = rnorm(12, 4, 0.3))
  spec <- data.frame(variable = c("g", "r"), transform = c("dlog", "level"),
                     target = c(TRUE, TRUE), year_ended = c(TRUE, FALSE),
                     stringsAsFactors = FALSE)
  cfg <- list()
  t <- 6
  dr <- array(rnorm(ndraw * H * 2, 0.5, 0.4), c(ndraw, H, 2),
              dimnames = list(NULL, NULL, c("g", "r")))
  oos <- list(list(origin = t, draws = dr))
  sc <- score_member("m", oos, td, spec, cfg)
  # cum exists for the dlog target only, h >= 2; never for the level target
  expect_true(all(sc$variable[sc$measure == "cum"] == "g"))
  expect_equal(sort(unique(sc$h[sc$measure == "cum"])), 2:4)
  expect_equal(nrow(sc[sc$measure == "cum" & sc$variable == "r", ]), 0)
  # the cumulative realized value at h=3 is the sum of the 3 quarterly outcomes
  cum3 <- sc[sc$measure == "cum" & sc$h == 3, ]
  expect_equal(cum3$real, sum(td$g[(t + 1):(t + 3)]))
  # and its CRPS matches scoring the summed draws directly
  expect_equal(cum3$crps,
               crps_sample(sum(td$g[(t + 1):(t + 3)]), rowSums(dr[, 1:3, "g"])),
               tolerance = 1e-9)
  # at h=1 cum is omitted (equals q); q still present at all horizons
  expect_equal(sort(unique(sc$h[sc$measure == "q" & sc$variable == "g"])), 1:4)
})
