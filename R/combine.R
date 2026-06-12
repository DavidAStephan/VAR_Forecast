# combine.R -- density combination schemes -------------------------------------
#
# Weights are estimated per target variable and per horizon bucket
# (near/medium/far), using ONLY scores whose realization was observable at the
# forecast origin (s + h <= t) -- strictly recursive, no look-ahead. All
# schemes except BMA are shrunk toward equal weights (combination puzzle).
# The linear pool's log score uses log-sum-exp; its CRPS uses draws resampled
# from the mixture; its PIT is the weighted member PIT (exact mixture CDF).

#' Training scores available at origin t for variable v and horizons hs.
#' Restricted to (origin, h) cells where EVERY member has a score, so that
#' sums of log densities are comparable across members (a member with fewer
#' negative terms would otherwise get spuriously large weight).
.train_scores <- function(scores, v, hs, t, members) {
  s <- scores[scores$variable == v & scores$measure == "q" &
              scores$h %in% hs & (scores$origin + scores$h) <= t &
              scores$member %in% members, ]
  if (!nrow(s)) return(s)
  cell <- paste(s$origin, s$h)
  complete <- names(which(table(cell) == length(members)))
  s[cell %in% complete, ]
}

#' Discounted sum of log predictive densities per member.
.disc_logdens <- function(tr, t, forgetting) {
  age <- t - (tr$origin + tr$h)
  w <- forgetting^age
  vapply(split(seq_len(nrow(tr)), tr$member),
         function(ix) sum(w[ix] * tr$logdens[ix]), numeric(1))
}

#' Number of distinct training origins.
.n_train_origins <- function(tr) length(unique(tr$origin))

.shrink <- function(w, kappa) {
  w <- kappa / length(w) + (1 - kappa) * w
  w / sum(w)
}

#' Optimal prediction pool (Hall-Mitchell / Geweke-Amisano): maximise the
#' discounted historical log score of the linear pool over the simplex,
#' softmax-parametrised, numerically via BFGS.
optimal_pool_weights <- function(tr, members, t, forgetting) {
  # matrix of log densities: rows = (origin, h) obs, cols = members
  key <- paste(tr$origin, tr$h)
  obs <- unique(key)
  L <- matrix(NA_real_, length(obs), length(members),
              dimnames = list(obs, members))
  L[cbind(match(key, obs), match(tr$member, members))] <- tr$logdens
  keep <- rowSums(is.na(L)) == 0
  L <- L[keep, , drop = FALSE]
  if (nrow(L) < 4) return(setNames(rep(1 / length(members), length(members)), members))
  age <- t - (tr$origin + tr$h)[match(rownames(L), key)]
  dsc <- forgetting^age
  negll <- function(theta) {
    lw <- c(0, theta); lw <- lw - logsumexp(lw)
    -sum(dsc * apply(L, 1, function(ld) logsumexp(lw + ld)))
  }
  opt <- optim(rep(0, length(members) - 1), negll, method = "BFGS",
               control = list(maxit = 200))
  lw <- c(0, opt$par); lw <- lw - logsumexp(lw)
  setNames(exp(lw), members)
}

#' Compute weights for one (scheme, variable, bucket) at origin t.
combo_weights <- function(scheme, scores, v, hs, t, members, cfg) {
  n <- length(members)
  eq <- setNames(rep(1 / n, n), members)
  if (scheme == "equal") return(eq)
  tr <- .train_scores(scores, v, hs, t, members)
  if (.n_train_origins(tr) < cfg$combination$min_train_origins) return(eq)
  if (scheme == "logscore") {
    ld <- .disc_logdens(tr, t, cfg$combination$forgetting)
    ld <- ld[members]; ld[is.na(ld)] <- -Inf
    # softmax of discounted log-score sums (comparable across members because
    # .train_scores keeps only complete cells)
    w <- exp(ld - max(ld, na.rm = TRUE))
    w[!is.finite(w)] <- 0
    if (sum(w) <= 0) return(eq)
    return(.shrink(w / sum(w), cfg$combination$shrink_kappa))
  }
  if (scheme == "bma") {
    # predictive-likelihood BMA: one-step-ahead log scores, no forgetting,
    # no shrinkage -- reported as a which-model-does-the-data-favour diagnostic
    tr1 <- tr[tr$h == min(tr$h), ]
    ld <- vapply(split(tr1$logdens, tr1$member), sum, numeric(1))[members]
    ld[is.na(ld)] <- -Inf
    w <- exp(ld - max(ld))
    w[!is.finite(w)] <- 0
    if (sum(w) <= 0) return(eq)
    return(w / sum(w))
  }
  if (scheme == "pool") {
    w <- optimal_pool_weights(tr, members, t, cfg$combination$forgetting)
    return(.shrink(w, cfg$combination$shrink_kappa))
  }
  stop("unknown scheme: ", scheme)
}

#' Reconstruct year-ended draws for variable v, horizon h at origin t from
#' quarterly draws + realized history (mirrors score_member).
.ye_from_q <- function(dr_list_vh, td, v, h, t) {
  k <- (h - 3):h
  hist_part <- sum(td[[v]][t + k[k <= 0]])
  lapply(dr_list_vh, function(mat) {
    # mat: [draw, h_index] for horizons k >= 1
    rowSums(mat) + hist_part
  })
}

#' Run all combination schemes. Returns list(scores = combo score rows in the
#' same long format, weights = long weight table).
combine_all <- function(scores, draws_env, td, spec, cfg) {
  members <- vapply(all_members(cfg), `[[`, "", "name")
  tgt <- spec$variable[spec$target]
  buckets <- cfg$combination$horizon_buckets
  schemes <- unlist(cfg$combination$schemes)
  origins <- sort(unique(scores$origin))
  T_n <- nrow(td)
  out_rows <- list(); w_rows <- list()

  for (scheme in schemes) for (v in tgt) {
    vt <- spec[spec$variable == v, ]
    for (t in origins) {
      for (bn in names(buckets)) {
        hs <- unlist(buckets[[bn]])
        w <- combo_weights(scheme, scores, v, hs, t, members, cfg)
        w_rows[[length(w_rows) + 1]] <- data.frame(
          scheme = scheme, variable = v, bucket = bn, origin = t,
          member = names(w), weight = as.numeric(w))
        for (h in hs) {
          if (t + h > T_n) next
          sub <- scores[scores$origin == t & scores$variable == v &
                        scores$h == h & scores$member %in% members, ]
          for (msr in unique(sub$measure)) {
            ss <- sub[sub$measure == msr, ]
            ss <- ss[match(members, ss$member), ]
            if (anyNA(ss$logdens)) next
            lw <- log(pmax(w, 1e-12))
            ld_pool <- logsumexp(lw + ss$logdens)
            pit_pool <- sum(w * ss$pit)
            pt_pool <- sum(w * ss$point)
            # mixture CRPS via resampled draws
            key <- paste(t, v, h, sep = "|")
            dl <- if (exists(key, draws_env)) get(key, draws_env) else NULL
            crps_pool <- NA_real_
            if (!is.null(dl) && all(members %in% names(dl))) {
              set.seed(derive_seed(cfg$master_seed,
                                   paste("mix", scheme, v, t, h, msr)))
              n_mix <- 500
              cnt <- drop(rmultinom(1, n_mix, w))
              if (msr == "ye") {
                kpos <- ((h - 3):h); kpos <- kpos[kpos >= 1]
                dmat <- lapply(members, function(m) {
                  do.call(cbind, lapply(kpos, function(hh)
                    get(paste(t, v, hh, sep = "|"), draws_env)[[m]]))
                })
                names(dmat) <- members
                dl_use <- .ye_from_q(dmat, td, v, h, t)
              } else dl_use <- dl
              mix <- unlist(lapply(seq_along(members), function(i) {
                if (cnt[i] == 0) return(numeric(0))
                sample(dl_use[[members[i]]], cnt[i], replace = TRUE)
              }))
              crps_pool <- scoringRules::crps_sample(ss$real[1], mix)
            }
            out_rows[[length(out_rows) + 1]] <- data.frame(
              member = paste0("combo_", scheme), origin = t,
              date = ss$date[1], variable = v, measure = msr, h = h,
              point = pt_pool, real = ss$real[1], logdens = ld_pool,
              crps = crps_pool, pit = pit_pool)
          }
        }
      }
    }
  }
  list(scores = do.call(rbind, out_rows), weights = do.call(rbind, w_rows))
}
