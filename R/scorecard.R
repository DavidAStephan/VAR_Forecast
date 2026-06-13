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

  add("\n## 5. How to read this\n")
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
