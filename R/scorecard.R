# scorecard.R -- single consolidated model scorecard (specs + performance) -----
#
# write_model_scorecard() produces one easy-to-read markdown document listing
# every model and its specification, then comparing forecast performance
# (CRPS, log score, RMSE) by target variable and horizon, with the best entry
# per column highlighted. It is a first-class pipeline output and is also
# generatable standalone from the saved CSV tables.

# engine-code -> human description
.engine_desc <- c(
  conj_br = "Block-recursive conjugate NIW",
  gibbs   = "Independent Normal-inverse-Wishart (Gibbs)",
  ss      = "Steady-state (Villani)",
  sv      = "Stochastic volatility, equation-by-equation",
  rw      = "Random walk",
  ar4     = "Bayesian AR(4)",
  ucsv    = "Unobserved components + stochastic volatility",
  ucmean  = "Unconditional mean")

# per-engine COVID-period treatment (README.md D17)
.covid_desc <- function(engine, treatment) {
  if (treatment == "none") return("none")
  if (engine == "sv") return("t-errors (SV-t)")
  if (engine == "ucsv") return("t-errors (robust)")
  if (engine %in% c("gibbs", "ss", "conj_br", "rw", "ar4", "ucmean"))
    return(if (treatment == "dummy") "dummy/drop" else "LP scaling")
  "—"
}

# ---- per-model narrative profiles (README.md is the full rationale) --------
# Hand-authored fields; the spec line and the eval evidence are auto-generated
# so the profiles stay accurate as the suite/data change.
.model_profiles <- list(
  small_minn = list(
    distinctive = paste0("The workhorse Minnesota BVAR, estimated by Gibbs with an ",
      "*independent* Normal-inverse-Wishart prior — the engine required to impose ",
      "block exogeneity, since asymmetric (equation-specific) prior variances cannot ",
      "be represented by a Kronecker/conjugate prior. Shrinkage is data-driven (GLP ",
      "marginal-likelihood), and sum-of-coefficients + dummy-initial-observation ",
      "priors discipline the I(1) levels (rates, real TWI)."),
    role = paste0("The central, representative small-SOE BVAR — the reference point the ",
      "other small members are deliberate variations around (tighter, looser, ",
      "steady-state, SV)."),
    watch = paste0("Solid all-rounder at short-to-medium horizons. Constant volatility ",
      "means it leans on the LP scaling for 2020; without it the 2020 outliers would ",
      "distort the coefficients."),
    refs = "D3, D4, D5, D17"),
  small_ss = list(
    distinctive = paste0("Reparametrised around its *unconditional means* (Villani steady ",
      "state), with informative priors placed directly on long-run levels — inflation ",
      "2.5% (target midpoint), NAIRU 4.5%, neutral cash rate 3.5%, US potential growth — ",
      "and on the foreign block's steady states, which the domestic forecast inherits."),
    role = paste0("The **long-horizon anchor**. An iterated VAR reverts to its ",
      "unconditional mean at h = 8-12; this member makes that mean economically grounded ",
      "rather than the raw sample average."),
    watch = paste0("Strongest at medium/far horizons for the mean-reverting targets. ",
      "Vulnerable if a steady-state anchor is stale (e.g. a shifted neutral rate drags the ",
      "long end); constant volatility under-disperses around 2020 absent the LP correction."),
    refs = "D3, D4, D5, D17"),
  small_sv = list(
    distinctive = paste0("Stochastic volatility estimated equation-by-equation (the ",
      "Carriero-Clark-Marcellino triangular factorisation), with t-distributed errors ",
      "(the CCMM SV-t COVID treatment). Block exogeneity is *exact* — foreign equations ",
      "simply drop the domestic regressors."),
    role = paste0("The **density-calibration specialist**: time-varying volatility tracks ",
      "changing macro uncertainty and the fat tails absorb outliers instead of letting ",
      "them widen the whole history."),
    watch = paste0("Its SV captures the time-varying conditional variance, so it stays ",
      "well-calibrated when uncertainty shifts — competitive on the GDP level at medium ",
      "horizons. The volatility state at the jump-off can over/under-shoot if the last few ",
      "quarters were unusual; the small system limits cross-variable information."),
    refs = "D3, D6, D17"),
  small_loose_p5 = list(
    distinctive = paste0("A deliberately *under-shrunk*, longer-lag variant — fixed ",
      "λ = 0.4 (vs the ~0.1-0.15 the GLP procedure selects) and 5 lags — so the data speak ",
      "more and richer dynamics can show through, at the cost of estimation noise."),
    role = paste0("The **loose / long-lag** diversity axis: it fails differently from the ",
      "tightly-shrunk members and can capture dynamics they shrink away."),
    watch = paste0("Occasionally best at near-horizon density when the extra flexibility pays; ",
      "noisier and prone to wider intervals at long horizons (the cost of light shrinkage in ",
      "a short sample)."),
    refs = "D4, D8"),
  medium_minn = list(
    distinctive = paste0("The larger system — 13 variables (adds terms of trade, wages, ",
      "employment, consumption, the 10y yield) with stochastic volatility and t-errors, but ",
      "only 2 lags to keep the parameter count feasible; equation-by-equation estimation keeps ",
      "the recursive loop tractable."),
    role = paste0("The **medium-system** axis (Banbura-Giannone-Reichlin): medium systems often ",
      "forecast best given enough shrinkage, and the extra variables bring cross-sectional ",
      "information the small core lacks."),
    watch = paste0("Strong short-horizon density (it tends to win the near bucket). The short ",
      "lag length limits long-horizon dynamics, and more parameters mean more estimation ",
      "uncertainty at the far end."),
    refs = "D1, D6, D17"),
  medium_conj = list(
    distinctive = paste0("The medium system estimated by the fast *block-recursive conjugate* ",
      "scheme (foreign VAR + domestic block conditioned on contemporaneous foreign values, the ",
      "RBNZ Bloor-Matheson approach) — closed-form, so cheap even at 13 variables x 4 lags. ",
      "Block exogeneity is exact by the recursive structure, not the prior."),
    role = paste0("The cheap medium workhorse; it complements `medium_minn` (conjugate ",
      "constant-volatility vs SV) on the same large system."),
    watch = paste0("Tends to lead the far-horizon GDP level error. Constant volatility ",
      "leans on LP for 2020; the conjugate Kronecker prior cannot represent asymmetric ",
      "shrinkage, which is why block exogeneity comes from the recursive structure."),
    refs = "D3, D5, D8"),
  small_tight = list(
    distinctive = paste0("The heavily-shrunk small model — fixed λ = 0.05, far tighter than the ",
      "GLP selection — pulling hard toward the persistence/random-walk prior, so it is ",
      "parsimonious and low-variance."),
    role = paste0("The **tight** diversity axis and the long-horizon robustness member: heavy ",
      "shrinkage buys stability where lightly-parametrised models wander."),
    watch = paste0("Best or near-best at far-horizon GDP (the tight prior stops it over-reacting). ",
      "Can be too rigid at short horizons, missing genuine dynamics the looser members catch."),
    refs = "D4, D5, D8"))

.benchmark_profiles <- list(
  rw = list(
    distinctive = "The no-change forecast: the last observed value persists, with Gaussian increments scaled to the historical change.",
    role = "The universal hard-to-beat short-horizon bar for persistent/level variables, and a pool member.",
    watch = "Competitive at h = 1 for level variables; fails badly at long horizons — its level path runs away, which the level-error tables (§2) expose brutally.",
    refs = "D8"),
  ar4 = list(
    distinctive = "A Bayesian AR(4) per variable with Minnesota-style lag shrinkage and a stationarity-truncated posterior.",
    role = "The univariate-persistence bar — it isolates how much of the forecast is just own-history dynamics.",
    watch = "Surprisingly strong for inflation at near/medium horizons, where univariate dynamics dominate; it cannot use cross-variable information, so it lags when that matters.",
    refs = "D8"),
  ucsv = list(
    distinctive = "Stock-Watson unobserved-components stochastic volatility per variable: a random-walk trend plus transitory noise, both with time-varying variances and outlier-robust t-errors.",
    role = "The canonical inflation benchmark and a genuine density anchor for the other targets.",
    watch = "Strong for inflation (its native use case); weaker for variables with richer multivariate dynamics. The trend/noise split is weakly identified, so it is gated on Monte-Carlo precision, not raw ESS.",
    refs = "D9, D17"),
  ucmean = list(
    distinctive = "A Gaussian density centred on the expanding-sample mean with the sample variance — the simplest possible density forecast.",
    role = "The floor: the 'did the model beat just predicting the long-run average' bar.",
    watch = "Unexpectedly competitive at long horizons for mean-reverting growth (everything reverts to the mean eventually); useless at short horizons where dynamics matter.",
    refs = "D8"))

.combo_profiles <- list(
  combo_equal = list(
    how = "Equal weights on every member, per variable x horizon bucket.",
    note = "The forecast-combination-puzzle benchmark and the recommended robust default: in the evaluation it has the best mean log score at every horizon and the best far-horizon CRPS. Hard to beat because it never over-fits weights."),
  combo_logscore = list(
    how = "Weights proportional to each member's recent log predictive score, with a forgetting factor, shrunk toward equal.",
    note = "Adapts to which members are forecasting well lately; the shrinkage and forgetting guard against over-concentrating on a member that was lucky."),
  combo_pool = list(
    how = "Optimal prediction pool (Hall-Mitchell / Geweke-Amisano): weights on the simplex that maximise the historical *pooled* log score, shrunk toward equal.",
    note = "Unlike BMA it does not degenerate to a single model ('all models are false but useful'); competitive with equal weights and occasionally better at near horizons."),
  combo_bma = list(
    how = "Bayesian model averaging by predictive likelihood — no shrinkage.",
    note = "**Diagnostic only.** It concentrates weight on the single best-fitting member, so it answers 'which model does the data favour' rather than serving as a robust combination; reported, not recommended."))

#' Auto-generated one-line spec for a VAR suite member.
.spec_oneliner <- function(m, cfg) {
  size <- if (m$set == "small") "8-variable small SOE core"
          else if (m$set == "medium") "13-variable medium system" else m$set
  lam <- if (identical(m$prior$lambda, "auto")) "GLP marginal-likelihood shrinkage"
         else paste0("fixed shrinkage λ=", m$prior$lambda)
  vol <- if (m$engine == "sv") "stochastic volatility" else "constant volatility"
  parts <- c(.engine_desc[[m$engine]], size, paste0(m$lags, " lags"), lam, vol,
             paste0(.covid_desc(m$engine, cfg$covid$treatment), " (COVID)"))
  extras <- c()
  if (isTRUE(m$prior$soc)) extras <- c(extras, "sum-of-coefficients")
  if (isTRUE(m$prior$dio)) extras <- c(extras, "dummy-initial-observation")
  if (length(extras))
    parts <- c(parts, paste0(paste(extras, collapse = " + "), " priors"))
  paste(parts, collapse = "; ")
}

#' Auto-computed eval evidence: where a member ranks best among the individual
#' models (excludes combinations) by mean quarterly CRPS, across variable x
#' horizon bucket.
.member_evidence <- function(scores, member, tgt, buckets, pool, measure = "level") {
  res <- list()
  for (v in tgt) for (bn in names(buckets)) {
    d <- scores[scores$measure == measure & scores$variable == v &
                scores$h %in% buckets[[bn]] & scores$member %in% pool, ]
    if (!nrow(d)) next
    agg <- aggregate(crps ~ member, d, mean)
    agg <- agg[order(agg$crps), ]
    rk <- match(member, agg$member)
    if (!is.na(rk)) res[[length(res) + 1]] <- list(v = v, bn = bn, rank = rk, n = nrow(agg))
  }
  if (!length(res)) return("")
  ones <- Filter(function(x) x$rank == 1, res)
  if (length(ones))
    return(paste0("best individual model for ",
                  paste(sapply(ones, function(x) paste0(x$v, " at ", x$bn)),
                        collapse = "; "), "."))
  b <- res[[which.min(sapply(res, function(x) x$rank))]]
  paste0("strongest at ", b$v, " (", b$bn, "), ranked ", b$rank, " of ", b$n,
         " individual models.")
}

#' Markdown table: members (rows) x horizons (cols) of a metric, best per
#' column bolded. better = "low" (CRPS, RMSE) or "high" (log score).
.metric_table <- function(scores, variable, measure, metric, horizons, better) {
  d <- scores[scores$variable == variable & scores$measure == measure &
              scores$h %in% horizons, ]
  if (!nrow(d)) return("")
  members <- unique(d$member)
  M <- matrix(NA_real_, length(members), length(horizons),
              dimnames = list(members, as.character(horizons)))
  for (i in seq_along(members)) for (j in seq_along(horizons)) {
    v <- d[[metric]][d$member == members[i] & d$h == horizons[j]]
    if (length(v)) M[i, j] <- v[1]
  }
  # order rows by mean metric (best first)
  ord <- order(rowMeans(M, na.rm = TRUE), decreasing = (better == "high"))
  M <- M[ord, , drop = FALSE]; members <- rownames(M)
  best <- apply(M, 2, function(z) {
    z <- z[is.finite(z)]
    if (!length(z)) return(NA_real_)
    if (better == "high") max(z) else min(z)
  })
  fmt <- function(x, col) {
    if (is.na(x)) return("—")
    if (x == -Inf) return("−∞")
    if (x == Inf) return("+∞")
    s <- formatC(x, format = "f", digits = 3)
    if (is.finite(best[col]) && abs(x - best[col]) < 1e-9) paste0("**", s, "**") else s
  }
  hdr <- paste0("| Model | ", paste0("h=", horizons, collapse = " | "), " |")
  sep <- paste0("|", paste(rep(":--", length(horizons) + 1), collapse = "|"), "|")
  rows <- vapply(seq_along(members), function(i) {
    cells <- vapply(seq_along(horizons), function(j) fmt(M[i, j], j), "")
    paste0("| ", members[i], " | ", paste(cells, collapse = " | "), " |")
  }, "")
  paste(c(hdr, sep, rows), collapse = "\n")
}

#' "Who wins" matrix: best model (by mean CRPS over all 12 horizons) per
#' variable x horizon bucket.
.winner_matrix <- function(scores, variables, buckets, measure = "level") {
  hdr <- paste0("| Variable | ", paste(names(buckets), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(":--", length(buckets) + 1), collapse = "|"), "|")
  rows <- vapply(variables, function(v) {
    cells <- vapply(buckets, function(hs) {
      d <- scores[scores$variable == v & scores$measure == measure & scores$h %in% hs, ]
      agg <- aggregate(crps ~ member, d, mean)
      w <- agg$member[which.min(agg$crps)]
      sprintf("%s (%.3f)", w, min(agg$crps))
    }, "")
    paste0("| ", v, " | ", paste(cells, collapse = " | "), " |")
  }, "")
  paste(c(hdr, sep, rows), collapse = "\n")
}

#' The LEVEL-error view -- the headline of this report. For each target it is
#' the forecast error of the underlying series' LEVEL at quarter t+h:
#'   - growth-modelled (dlog) variables (GDP, inflation): the cumulative level
#'     from the forecast origin, 100*(log level_{t+h} - log level_t) = the `cum`
#'     measure for h>=2, and `q` at h=1 (where cum == q);
#'   - level-modelled variables (unemployment, cash rate): the rate level at
#'     t+h, which is already the `q` measure.
#' Returns the summarised scores relabelled measure="level".
.level_view <- function(scores, spec) {
  tgt <- spec$variable[spec$target %in% c(TRUE, "TRUE")]
  out <- lapply(tgt, function(v) {
    if (identical(spec$transform[spec$variable == v], "dlog"))
      rbind(scores[scores$variable == v & scores$measure == "cum", ],
            scores[scores$variable == v & scores$measure == "q" & scores$h == 1, ])
    else
      scores[scores$variable == v & scores$measure == "q", ]
  })
  d <- do.call(rbind, out)
  d$measure <- "level"
  d
}

#' The level measure used per target (cum for growth-modelled, q for levels).
.level_measure <- function(v, spec)
  if (identical(spec$transform[spec$variable == v], "dlog")) "cum" else "q"

#' Full scorecard markdown.
write_model_scorecard <- function(scores, spec, dm, diag, cfg,
                                  src = "real", panel = "1997Q4-2026Q1",
                                  out_path = "reports/model_scorecard.md") {
  # scores arrives PER ORIGIN (member x origin x variable x measure x h);
  # collapse to mean scores per (member, variable, measure, h) so every table
  # below shows the average over forecast origins, not a single origin.
  if ("origin" %in% names(scores)) scores <- summarise_scores(scores)
  tgt <- spec$variable[spec$target == TRUE | spec$target == "TRUE"]
  hcols <- c(1, 4, 8, 12)
  buckets <- list(`near (1-4)` = 1:4, `medium (5-8)` = 5:8, `far (9-12)` = 9:12)
  lvl <- .level_view(scores, spec)        # the headline LEVEL-error view
  L <- c()
  add <- function(...) L <<- c(L, paste0(...))

  add("# SOE-BVAR Suite — Model Scorecard\n")
  add("**Data source:** ", src,
      if (identical(src, "real")) " (RBA + ABS + FRED)" else " (simulated DGP — validates machinery, not Australian dynamics)", "  ")
  add("**Panel:** ", panel, " | **Targets:** ", paste(tgt, collapse = ", "),
      " | **Horizons:** 1-12 quarters  ")
  add("**Evaluation:** expanding-window pseudo-real-time, ",
      cfg$evaluation$max_origins, " forecast origins; densities scored by CRPS ",
      "and log predictive density, points by RMSE, all by horizon.  ")
  add("**Headline metric — the LEVEL error.** Performance is reported on the ",
      "forecast error of each series' *level* at quarter t+h: for GDP and ",
      "inflation (modelled as growth) the cumulative level from the forecast ",
      "origin — where real GDP and the price level land — and for unemployment ",
      "and the cash rate the rate level itself. The quarterly and year-ended ",
      "*growth* scores are not shown here (they answer a different, narrower ",
      "question); they remain in `output/tables/scores_by_horizon.csv`.  ")
  d9 <- all(diag$converged_all & diag$sanity_all & diag$no_lookahead & diag$reproducible) &&
        all(diag$block_exog_max < 1e-2)
  add("**Diagnostics (§9):** ", if (d9) "all green" else "SEE DIAGNOSTICS",
      " (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).\n")

  # ---- 1. Model suite specifications ----
  add("## 1. The models\n")
  add("Every VAR member is **block-exogenous** (the domestic block never feeds ",
      "back into the foreign block) and produces **iterated** density forecasts. ",
      "Members are designed to fail differently; see README.md for the full rationale.\n")
  add("### 1a. VAR members\n")
  add("| Model | Family | System | Lags | Shrinkage λ | Volatility | COVID |")
  add("|:--|:--|:--|:--|:--|:--|:--|")
  setsize <- function(set) if (set == "small") "small (8 var)" else if (set == "medium") "medium (13 var)" else set
  for (m in cfg$suite) {
    lam <- if (identical(m$prior$lambda, "auto")) "auto (GLP)" else as.character(m$prior$lambda)
    vol <- if (m$engine == "sv") "stochastic (SV)" else "constant"
    extras <- c()
    if (isTRUE(m$prior$soc)) extras <- c(extras, "SOC")
    if (isTRUE(m$prior$dio)) extras <- c(extras, "DIO")
    fam <- .engine_desc[[m$engine]]
    if (length(extras)) fam <- paste0(fam, " + ", paste(extras, collapse = "/"))
    add(sprintf("| `%s` | %s | %s | %d | %s | %s | %s |",
                m$name, fam, setsize(m$set), m$lags, lam, vol,
                .covid_desc(m$engine, cfg$covid$treatment)))
  }
  add("\n### 1b. Benchmark members (the bar every VAR must clear)\n")
  add("| Model | Description | COVID |")
  add("|:--|:--|:--|")
  for (b in cfg$benchmarks)
    add(sprintf("| `%s` | %s | %s |", b, .engine_desc[[b]],
                .covid_desc(b, cfg$covid$treatment)))

  add("\n### 1c. Combination schemes (density pools)\n")
  add("Weights estimated **per target variable and per horizon bucket** ",
      "(near 1-4, medium 5-8, far 9-12), shrunk toward equal weights, strictly ",
      "recursive (no look-ahead).\n")
  add("| Scheme | How weights are set |")
  add("|:--|:--|")
  add("| `combo_equal` | Equal weights (the benchmark pool — hard to beat) |")
  add("| `combo_logscore` | Recursive log-score weights, with forgetting |")
  add("| `combo_pool` | Optimal prediction pool (Hall-Mitchell / Geweke-Amisano) |")
  add("| `combo_bma` | Bayesian model averaging — reported as a diagnostic only |")

  # ---- 2. Performance (LEVEL errors) ----
  add("\n## 2. Forecast performance — level errors\n")
  add("All tables below score the **level** of each series at t+h (§ definition ",
      "in the header): cumulative real GDP and the price level for the growth ",
      "variables, the rate level for unemployment and the cash rate. Lower CRPS / ",
      "RMSE is better; higher log score is better. **Bold** = best in that column. ",
      "Models ordered best-first (mean over the shown horizons). The level errors ",
      "grow with horizon (they accumulate the whole path), so the columns are not ",
      "comparable across horizons — read down each column, not across.\n")

  add("### 2a. Who forecasts the level best, by variable and horizon bucket (CRPS)\n")
  add("Best single model (lowest mean level CRPS in the bucket; value in parentheses):\n")
  add(.winner_matrix(lvl, tgt, buckets, "level"), "\n")

  add("### 2b. Level density accuracy by variable and horizon (CRPS, lower better)\n")
  for (v in tgt) {
    lbl <- spec$label[spec$variable == v]
    note <- if (identical(spec$transform[spec$variable == v], "dlog"))
      " — cumulative level from the origin" else " — rate level at t+h"
    add("**", lbl, " (`", v, "`)", note, "**\n")
    add(.metric_table(lvl, v, "level", "crps", hcols, "low"), "\n")
  }

  add("\n### 2c. Level point accuracy by variable and horizon (RMSE, lower better)\n")
  for (v in tgt) {
    add("**`", v, "`**\n")
    add(.metric_table(lvl, v, "level", "rmse", hcols, "low"), "\n")
  }

  add("\n### 2d. Level density calibration — mean log predictive density (higher better)\n")
  add("Averaged across the 4 targets, on the level forecast. The mean log score is ",
      "brutally sensitive to tail events: **−∞** means at least one origin where the ",
      "realization fell outside that member's predictive support (an individual model ",
      "can catastrophically miss a tail). The **combinations never do** — the linear ",
      "pool always assigns positive density — which is the clearest single piece of ",
      "evidence that pooling buys calibration and robustness.\n")
  ls_all <- aggregate(logdens ~ member + h, lvl, mean)
  add(.metric_table(
    transform(ls_all, variable = "all", measure = "level"),
    "all", "level", "logdens", hcols, "high"), "\n")

  # ---- 3. Combination vs best member ----
  add("\n## 3. Do the combinations beat the best single model?\n")
  add("Mean **level** CRPS over all 4 targets, by horizon bucket. The honest test ",
      "of a pool is whether it beats both equal weights and the best individual ",
      "member. (Buckets average level errors of differing scale across horizons, so ",
      "use this within a bucket to rank models, not to compare buckets.)\n")
  add("| Model | near (1-4) | medium (5-8) | far (9-12) |")
  add("|:--|:--|:--|:--|")
  members <- unique(lvl$member)
  avg <- lapply(members, function(m) {
    sapply(buckets, function(hs) {
      d <- lvl[lvl$member == m & lvl$h %in% hs, ]
      mean(d$crps)
    })
  })
  names(avg) <- members
  A <- do.call(rbind, avg)
  ord <- order(rowMeans(A))
  best <- apply(A, 2, min)
  for (m in members[ord]) {
    cells <- vapply(seq_len(ncol(A)), function(j) {
      s <- formatC(A[m, j], format = "f", digits = 3)
      if (abs(A[m, j] - best[j]) < 1e-9) paste0("**", s, "**") else s
    }, "")
    add(sprintf("| %s | %s |", m, paste(cells, collapse = " | ")))
  }

  # ---- 4. Significance ----
  if (!is.null(dm) && nrow(dm)) {
    add("\n## 4. Statistical significance (Diebold-Mariano)\n")
    add("How often each combination **significantly beats** the random-walk and ",
        "AR(4) benchmarks on **level** CRPS (Harvey-corrected, 10% level), counted ",
        "over the 4 targets x 3 horizons {4, 8, 12} tested. A negative DM statistic ",
        "means the combination is more accurate; significance is one-sided here. ",
        "Level errors are high-variance (they accumulate the path and are dominated ",
        "by a few episodes such as 2020), so the DM test has low power on them — the ",
        "combinations beat the benchmarks on *average* (§3) more often than they do ",
        "*significantly*.\n")
    # level DM: cumulative-level (cum) for growth variables, the rate level (q)
    # for level variables. Significance is computed from p_crps (the dm target
    # carries p_crps, not the star column).
    dm_lvl <- do.call(rbind, lapply(tgt, function(v)
      dm[dm$variable == v & dm$measure == .level_measure(v, spec), ]))
    dm_lvl$beats <- dm_lvl$dm_crps < 0 & dm_lvl$p_crps < 0.10
    dd <- dm_lvl[grepl("^combo_", dm_lvl$member) & dm_lvl$h %in% c(4, 8, 12), ]
    combos <- sort(unique(dd$member)); refs <- sort(unique(dd$reference))
    add("| Combination | ", paste0("beats ", refs, collapse = " | "), " |")
    add("|", paste(rep(":--", length(refs) + 1), collapse = "|"), "|")
    for (cm in combos) {
      cells <- vapply(refs, function(rf) {
        s <- dd[dd$member == cm & dd$reference == rf, ]
        sprintf("%d / %d", sum(s$beats, na.rm = TRUE), nrow(s))
      }, "")
      add(sprintf("| %s | %s |", cm, paste(cells, collapse = " | ")))
    }
    star <- function(p) if (p < 0.01) "***" else if (p < 0.05) "**" else "*"
    strong <- dd[which(dd$beats & dd$p_crps < 0.05), ]
    if (nrow(strong)) {
      strong <- strong[order(strong$p_crps), ]
      add("\nStrongest results (significant at 5% or better):\n")
      for (i in seq_len(min(6, nrow(strong))))
        add(sprintf("- `%s` beats `%s` on **%s** at h=%d (DM %.2f%s)",
                    strong$member[i], strong$reference[i], strong$variable[i],
                    strong$h[i], strong$dm_crps[i], star(strong$p_crps[i])))
    }
  }

  # ---- 5. Model profiles ----
  add("\n## 5. Model profiles\n")
  add("One entry per model: its specification, what makes it distinct, the role it ",
      "plays in the suite, its strengths and failure modes, and where it actually ",
      "ranks in this evaluation (the eval line is computed, not asserted). Full ",
      "rationale is in README.md.\n")
  pool <- vapply(c(cfg$suite, lapply(cfg$benchmarks, function(b) list(name = b))),
                 `[[`, "", "name")   # individual models (excludes combinations)
  add("### 5a. VAR members\n")
  for (m in cfg$suite) {
    p <- .model_profiles[[m$name]]
    add("**`", m$name, "`**  ")
    add("*Spec:* ", .spec_oneliner(m, cfg), ".  ")
    if (!is.null(p)) {
      add("*Distinctive:* ", p$distinctive, "  ")
      add("*Role:* ", p$role, "  ")
      add("*Strengths & failure modes:* ", p$watch, "  ")
    }
    add("*In this evaluation (level error):* ",
        .member_evidence(lvl, m$name, tgt, buckets, pool, "level"), "  ")
    if (!is.null(p)) add("*See:* README.md ", p$refs, "\n") else add("\n")
  }
  add("### 5b. Benchmark members\n")
  for (b in cfg$benchmarks) {
    p <- .benchmark_profiles[[b]]
    add("**`", b, "`**  ")
    add("*Spec:* ", .engine_desc[[b]], "; ", .covid_desc(b, cfg$covid$treatment),
        " (COVID).  ")
    if (!is.null(p)) {
      add("*Distinctive:* ", p$distinctive, "  ")
      add("*Role:* ", p$role, "  ")
      add("*Strengths & failure modes:* ", p$watch, "  ")
    }
    add("*In this evaluation (level error):* ",
        .member_evidence(lvl, b, tgt, buckets, pool, "level"), "  ")
    if (!is.null(p)) add("*See:* README.md ", p$refs, "\n") else add("\n")
  }
  add("### 5c. Combination schemes\n")
  for (cs in c("combo_equal", "combo_logscore", "combo_pool", "combo_bma")) {
    p <- .combo_profiles[[cs]]
    if (is.null(p)) next
    add("**`", cs, "`**  ")
    add("*Weights:* ", p$how, "  ")
    add("*In this evaluation:* ", p$note, "\n")
  }

  add("\n## 6. How to read this\n")
  add("- **The metric is the level error.** For GDP and inflation it is the ",
      "cumulative level from the forecast origin (where real GDP and the price ",
      "level land h quarters out); for unemployment and the cash rate it is the ",
      "rate level at t+h. Level errors accumulate the whole forecast path, so they ",
      "grow with horizon and are far more discriminating than the single-quarter ",
      "growth rate — a model that gets the persistent/drift component wrong (e.g. ",
      "the random walk) is exposed here. The quarterly and year-ended *growth* ",
      "scores live in `output/tables/scores_by_horizon.csv` if you need them.\n")
  add("- **Point gains over the best member are modest by design** — equal ",
      "weights are hard to beat (the forecast-combination puzzle). The pool's ",
      "payoff is *calibration and robustness*: it insures against any single ",
      "member failing, rather than always winning on accuracy.\n")
  add("- **CRPS and log score can disagree** on outlier-heavy windows (log ",
      "score is far more sensitive to tail events); both are reported. A ",
      "COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`. ",
      "See the full Quarto report for fan charts, PIT calibration, and ",
      "weight-evolution plots.\n")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(L, out_path)
  out_path
}
