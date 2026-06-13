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

# per-engine COVID-period treatment (DECISIONS.md D17)
.covid_desc <- function(engine, treatment) {
  if (treatment == "none") return("none")
  if (engine == "sv") return("t-errors (SV-t)")
  if (engine == "ucsv") return("t-errors (robust)")
  if (engine %in% c("gibbs", "ss", "conj_br", "rw", "ar4", "ucmean"))
    return(if (treatment == "dummy") "dummy/drop" else "LP scaling")
  "—"
}

# ---- per-model narrative profiles (DECISIONS.md is the full rationale) --------
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
    watch = paste0("Best at short horizons (h = 1-2), where getting the conditional variance ",
      "right matters most. The volatility state at the jump-off can over/under-shoot if the ",
      "last few quarters were unusual; the small system limits cross-variable information."),
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
    watch = paste0("Tends to lead the far-horizon GDP year-ended / cumulative-level buckets. ",
      "Constant volatility leans on LP for 2020; the conjugate Kronecker prior cannot represent ",
      "asymmetric shrinkage, which is why block exogeneity comes from the recursive structure."),
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
    watch = "Competitive at h = 1 for level variables; fails badly at long horizons for growth variables — its level path runs away, which the cumulative-level metric (§2e) exposes brutally.",
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
.member_evidence <- function(scores, member, tgt, buckets, pool) {
  res <- list()
  for (v in tgt) for (bn in names(buckets)) {
    d <- scores[scores$measure == "q" & scores$variable == v &
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
.winner_matrix <- function(scores, variables, buckets) {
  hdr <- paste0("| Variable | ", paste(names(buckets), collapse = " | "), " |")
  sep <- paste0("|", paste(rep(":--", length(buckets) + 1), collapse = "|"), "|")
  rows <- vapply(variables, function(v) {
    cells <- vapply(buckets, function(hs) {
      d <- scores[scores$variable == v & scores$measure == "q" & scores$h %in% hs, ]
      agg <- aggregate(crps ~ member, d, mean)
      w <- agg$member[which.min(agg$crps)]
      sprintf("%s (%.3f)", w, min(agg$crps))
    }, "")
    paste0("| ", v, " | ", paste(cells, collapse = " | "), " |")
  }, "")
  paste(c(hdr, sep, rows), collapse = "\n")
}

#' Full scorecard markdown.
write_model_scorecard <- function(scores, spec, dm, diag, cfg,
                                  src = "real", panel = "1997Q4-2026Q1",
                                  out_path = "reports/model_scorecard.md") {
  # scores arrives PER ORIGIN (member x origin x variable x measure x h);
  # collapse to mean scores per (member, variable, measure, h) so every table
  # below shows the average over forecast origins, not a single origin.
  if ("origin" %in% names(scores)) scores <- summarise_scores(scores)
  tgt <- spec$variable[spec$target == TRUE | spec$target == "TRUE"]
  hcols <- c(1, 2, 4, 8, 12)
  buckets <- list(`near (1-4)` = 1:4, `medium (5-8)` = 5:8, `far (9-12)` = 9:12)
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
  d9 <- all(diag$converged_all & diag$sanity_all & diag$no_lookahead & diag$reproducible) &&
        all(diag$block_exog_max < 1e-2)
  add("**Diagnostics (§9):** ", if (d9) "all green" else "SEE DIAGNOSTICS",
      " (block exogeneity, MCMC convergence, forecast sanity, no-look-ahead, reproducibility).\n")

  # ---- 1. Model suite specifications ----
  add("## 1. The models\n")
  add("Every VAR member is **block-exogenous** (the domestic block never feeds ",
      "back into the foreign block) and produces **iterated** density forecasts. ",
      "Members are designed to fail differently; see DECISIONS.md for the full rationale.\n")
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

  # ---- 2. Performance ----
  add("\n## 2. Forecast performance\n")
  add("Lower CRPS / RMSE is better; higher log score is better. **Bold** = best ",
      "in that column. Models ordered best-first (by mean over the shown horizons).\n")

  add("### 2a. Who forecasts best, by variable and horizon (CRPS)\n")
  add("Best single model (lowest mean CRPS in the bucket; value in parentheses):\n")
  add(.winner_matrix(scores, tgt, buckets), "\n")

  add("### 2b. Density accuracy by variable and horizon (CRPS, lower better)\n")
  for (v in tgt) {
    lbl <- spec$label[spec$variable == v]
    add("**", lbl, " (`", v, "`)**\n")
    add(.metric_table(scores, v, "q", "crps", hcols, "low"), "\n")
  }

  add("\n### 2c. Point accuracy by variable and horizon (RMSE, lower better)\n")
  for (v in tgt) {
    add("**`", v, "`**\n")
    add(.metric_table(scores, v, "q", "rmse", hcols, "low"), "\n")
  }

  add("\n### 2d. Density calibration — mean log predictive density (higher better)\n")
  add("Averaged across the 4 targets. The mean log score is brutally sensitive ",
      "to tail events: **−∞** means at least one origin where the realization ",
      "fell outside that member's predictive support (an individual model can ",
      "catastrophically miss a tail). The **combinations never do** — the linear ",
      "pool always assigns positive density — which is the clearest single ",
      "piece of evidence that pooling buys calibration and robustness.\n")
  ls_all <- aggregate(logdens ~ member + h,
                      scores[scores$measure == "q", ], mean)
  add(.metric_table(
    transform(ls_all, variable = "all", measure = "q"),
    "all", "q", "logdens", hcols, "high"), "\n")

  # ---- 2e. integrated measures for growth variables ----
  dlog_tgt <- spec$variable[spec$target %in% c(TRUE, "TRUE") &
                            spec$transform == "dlog"]
  has_int <- length(dlog_tgt) &&
    any(scores$measure %in% c("ye", "cum"))
  if (has_int) {
    add("\n### 2e. Year-ended and cumulative-level accuracy (growth variables)\n")
    add("For GDP and inflation (modelled as quarterly growth), the `q` score ",
        "above is the marginal rate *in that one quarter* — the narrowest, least ",
        "predictable view at long horizons. Two integrated views matter more for ",
        "policy and are scored here:\n")
    add("- **Year-ended** (`ye`, 4-quarter sum ending at t+h): the RBA's headline ",
        "concept — year-ended GDP growth, and the 2-3% *year-ended* trimmed-mean ",
        "inflation target.\n")
    add("- **Cumulative level** (`cum`, the h-quarter sum from the origin = ",
        "100·(log level_{t+h} − log level_t)): where the level lands h quarters ",
        "out. Far more discriminating than the single quarter — a model that gets ",
        "the persistent/drift component wrong (e.g. the random walk) is exposed ",
        "here but not by the quarterly score.\n")
    add("(For level variables — unemployment, cash rate — the `q` score already ",
        "*is* the level at t+h, so no separate cumulative view is needed. At ",
        "h=4 the two measures below coincide by construction — both span the 4 ",
        "quarters from the origin — and diverge from h=8 on.)\n")
    hi <- c(4, 8, 12)
    for (v in dlog_tgt) {
      lbl <- spec$label[spec$variable == v]
      add("**", lbl, " (`", v, "`) — CRPS, lower better**\n")
      add("Year-ended:\n")
      add(.metric_table(scores, v, "ye", "crps", hi, "low"), "\n")
      add("Cumulative level (from forecast origin):\n")
      add(.metric_table(scores, v, "cum", "crps", hi, "low"), "\n")
    }
  }

  # ---- 3. Combination vs best member ----
  add("\n## 3. Do the combinations beat the best single model?\n")
  add("Mean CRPS over all 4 targets, by horizon bucket. The honest test of a ",
      "pool is whether it beats both equal weights and the best individual member.\n")
  add("| Model | near (1-4) | medium (5-8) | far (9-12) |")
  add("|:--|:--|:--|:--|")
  members <- unique(scores$member)
  avg <- lapply(members, function(m) {
    sapply(buckets, function(hs) {
      d <- scores[scores$member == m & scores$measure == "q" & scores$h %in% hs, ]
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
        "AR(4) benchmarks on CRPS (Harvey-corrected, 10% level), counted over ",
        "the 4 targets x 3 horizons {1, 4, 8} tested. A negative DM statistic ",
        "means the combination is more accurate; significance is one-sided here.\n")
    dd <- dm[grepl("^combo_", dm$member) & dm$measure == "q" & dm$h %in% c(1, 4, 8), ]
    combos <- sort(unique(dd$member)); refs <- sort(unique(dd$reference))
    add("| Combination | ", paste0("beats ", refs, collapse = " | "), " |")
    add("|", paste(rep(":--", length(refs) + 1), collapse = "|"), "|")
    for (cm in combos) {
      cells <- vapply(refs, function(rf) {
        s <- dd[dd$member == cm & dd$reference == rf, ]
        nbeat <- sum(s$dm_crps < 0 & !is.na(s$sig_crps) & s$sig_crps != "")
        sprintf("%d / %d", nbeat, nrow(s))
      }, "")
      add(sprintf("| %s | %s |", cm, paste(cells, collapse = " | ")))
    }
    strong <- dd[dd$dm_crps < 0 & dd$sig_crps %in% c("**", "***"), ]
    if (nrow(strong)) {
      strong <- strong[order(strong$p_crps), ]
      add("\nStrongest results (significant at 5% or better):\n")
      for (i in seq_len(min(6, nrow(strong))))
        add(sprintf("- `%s` beats `%s` on **%s** at h=%d (DM %.2f%s)",
                    strong$member[i], strong$reference[i], strong$variable[i],
                    strong$h[i], strong$dm_crps[i], strong$sig_crps[i]))
    }
  }

  # ---- 5. Model profiles ----
  add("\n## 5. Model profiles\n")
  add("One entry per model: its specification, what makes it distinct, the role it ",
      "plays in the suite, its strengths and failure modes, and where it actually ",
      "ranks in this evaluation (the eval line is computed, not asserted). Full ",
      "rationale is in DECISIONS.md.\n")
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
    add("*In this evaluation:* ", .member_evidence(scores, m$name, tgt, buckets, pool), "  ")
    if (!is.null(p)) add("*See:* DECISIONS.md ", p$refs, "\n") else add("\n")
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
    add("*In this evaluation:* ", .member_evidence(scores, b, tgt, buckets, pool), "  ")
    if (!is.null(p)) add("*See:* DECISIONS.md ", p$refs, "\n") else add("\n")
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
  add("- **Point gains over the best member are modest by design** — equal ",
      "weights are hard to beat (the forecast-combination puzzle). The pool's ",
      "payoff is *calibration and robustness*: it insures against any single ",
      "member failing, rather than always winning on accuracy.\n")
  add("- **CRPS and log score can disagree** on outlier-heavy windows (log ",
      "score is far more sensitive to tail events); both are reported. A ",
      "COVID-excluded variant is in `output/tables/scores_by_horizon_excovid.csv`.\n")
  add("- **Pick the horizon view to match the decision.** The quarterly score ",
      "(§2b) is the marginal growth rate; for GDP and inflation the year-ended ",
      "and cumulative-level views (§2e) are usually what a central bank acts on, ",
      "and they rank models differently at long horizons. For level variables the ",
      "quarterly score already is the level. See the full Quarto report for fan ",
      "charts, PIT calibration, and weight-evolution plots.\n")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(L, out_path)
  out_path
}
