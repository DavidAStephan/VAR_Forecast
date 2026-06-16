# aaa_logging.R -- logging facade (makes `logger` an OPTIONAL dependency) -------
#
# The suite logs through bare log_info()/log_warn()/... calls rather than
# `logger::...` directly, so the `logger` package is no longer required to run.
#   * When `logger` is installed (the usual dev setup) every call forwards to it
#     unchanged, preserving the configured layout and namespace behaviour.
#   * When it is NOT installed (e.g. the RBA work environment) the calls fall
#     back to base-R console messages in the same "HH:MM:SS [LEVEL] msg" layout
#     setup_logging() configures for logger, with glue interpolation and the
#     same severity threshold.
# The "aaa_" prefix makes this file sort first in R/, so the facade exists
# before anything (or any parallel worker, via ensure_project_loaded()) calls
# it. It is intentionally NOT in config_hash()'s estimation-file list: logging
# is a pure side effect and must not invalidate the OOS cache.

# Honour an explicit opt-out (force the fallback even where logger exists) so
# the work-environment path can be exercised on a machine that has logger.
.logger_available <- function()
  !isTRUE(getOption("soe.no_logger", FALSE)) &&
    requireNamespace("logger", quietly = TRUE)

# fallback severity threshold (numeric, logger's scale); default INFO.
.log_state <- new.env(parent = emptyenv())
.log_state$threshold <- 400L
.LOG_SEVERITY <- c(TRACE = 100L, DEBUG = 200L, INFO = 400L, SUCCESS = 350L,
                   WARN = 500L, ERROR = 600L, FATAL = 700L)

.log_fallback <- function(level, severity, ..., .envir) {
  if (severity < .log_state$threshold) return(invisible(NULL))
  tmpl <- paste0(...)
  msg <- if (requireNamespace("glue", quietly = TRUE))
    tryCatch(as.character(glue::glue(tmpl, .envir = .envir)),
             error = function(e) tmpl)
  else tmpl
  message(sprintf("%s [%s] %s", format(Sys.time(), "%H:%M:%S"), level,
                  paste(msg, collapse = " ")))
  invisible(NULL)
}

# Build one log_<level>() function: delegate to logger if present (forwarding
# the caller's environment so glue templates resolve), else use the fallback.
.mk_log <- function(level, severity) {
  force(level); force(severity)
  lc <- tolower(level)
  function(...) {
    if (.logger_available()) {
      fn <- get(paste0("log_", lc), envir = asNamespace("logger"))
      fn(..., .topenv = parent.frame())
    } else .log_fallback(level, severity, ..., .envir = parent.frame())
  }
}

log_trace   <- .mk_log("TRACE",   100L)
log_debug   <- .mk_log("DEBUG",   200L)
log_info    <- .mk_log("INFO",    400L)
log_success <- .mk_log("SUCCESS", 350L)
log_warn    <- .mk_log("WARN",    500L)
log_error   <- .mk_log("ERROR",   600L)
log_fatal   <- .mk_log("FATAL",   700L)

# setup_logging() calls these; mirror logger's API for the fallback.
log_threshold <- function(level = "INFO", ...) {
  if (.logger_available()) return(logger::log_threshold(level, ...))
  sev <- .LOG_SEVERITY[[toupper(as.character(level))]]
  if (!is.null(sev)) .log_state$threshold <- sev
  invisible(NULL)
}
log_layout <- function(layout, ...) {
  if (.logger_available()) return(logger::log_layout(layout, ...))
  invisible(NULL)   # the fallback uses a fixed layout matching the project's
}
layout_glue_generator <- function(format = NULL, ...) {
  if (.logger_available()) return(logger::layout_glue_generator(format = format, ...))
  function(...) ""  # dummy; unused by the fallback
}
