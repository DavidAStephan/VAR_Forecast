# transforms.R -- raw levels -> model units, plus shape/coverage self-checks ----

#' Apply the spec transform to one raw series.
apply_transform <- function(x, transform) {
  switch(transform,
    dlog     = c(NA, 100 * diff(log(x))),
    loglevel = 100 * log(x),
    level    = x,
    stop("unknown transform: ", transform)
  )
}

#' Transform the raw panel into model units. Drops the leading NA row created
#' by differencing and trims to the common non-missing sample.
transform_data <- function(raw, spec) {
  df <- raw$data
  out <- data.frame(date = df$date)
  for (i in seq_len(nrow(spec))) {
    v <- spec$variable[i]
    stopifnot(v %in% names(df))
    out[[v]] <- apply_transform(df[[v]], spec$transform[i])
  }
  # common sample: first row where everything is finite, through last
  ok <- apply(out[, -1, drop = FALSE], 1, function(r) all(is.finite(r)))
  first <- which(ok)[1]
  last  <- max(which(ok))
  out <- out[first:last, , drop = FALSE]
  # interior NAs (mixed publication lags in real data): trim trailing rows
  # until the panel is balanced rather than imputing
  while (nrow(out) > 0 &&
         !all(is.finite(as.matrix(out[nrow(out), -1])))) out <- out[-nrow(out), ]
  if (anyNA(out[, -1])) {
    bad <- names(out[, -1])[colSums(!is.finite(as.matrix(out[, -1]))) > 0]
    stop("interior missing values in: ", paste(bad, collapse = ", "))
  }
  rownames(out) <- NULL
  attr(out, "dgp") <- attr(raw, "dgp")
  out
}

#' Self-checks on the transformed panel (data-layer diagnostics).
check_data <- function(td, spec, min_quarters = 80) {
  checks <- list()
  checks$balanced  <- !anyNA(td[, -1])
  checks$coverage  <- nrow(td) >= min_quarters
  checks$ordering  <- identical(names(td)[-1], spec$variable)
  rates <- spec$variable[spec$transform == "level"]
  checks$rate_range <- all(vapply(rates, function(v)
    all(td[[v]] > -5 & td[[v]] < 30), logical(1)))
  growth <- spec$variable[spec$transform == "dlog"]
  checks$growth_range <- all(vapply(growth, function(v)
    all(abs(td[[v]]) < 50), logical(1)))
  failed <- names(checks)[!unlist(checks)]
  if (length(failed)) stop("data checks failed: ", paste(failed, collapse = ", "))
  log_info("data checks passed ({nrow(td)} quarters, {ncol(td)-1} vars)")
  invisible(checks)
}

#' Year-ended transform of a quarterly dlog series: rolling 4-quarter sum.
#' For 'level' variables the year-ended concept is the level itself.
year_ended <- function(x, transform) {
  if (transform == "dlog") {
    n <- length(x)
    if (n < 4) return(rep(NA_real_, n))
    c(rep(NA_real_, 3), x[4:n] + x[3:(n - 1)] + x[2:(n - 2)] + x[1:(n - 3)])
  } else x
}
