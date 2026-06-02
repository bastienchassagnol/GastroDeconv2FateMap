#' Patch vscDebugger promise inspection for R 4.6+
#'
#' vscDebugger 0.5.6 uses removed C APIs (`PRVALUE`, `PRCODE`, `Rf_findVar`)
#' that segfault on R 4.6.0. The debugger then crashes while refreshing the
#' variables pane (often on an unrelated line such as `stats::cor()`).
#'
#' Until vscDebugger is rebuilt for R 4.6, replace `isPromise()` with a safe
#' R implementation that skips the broken `.Call(c_is_promise, ...)`.
#'
#' @return `TRUE` if the patch was applied, `FALSE` otherwise (invisible).
#' @keywords internal
patch_vscdebugger_promises_r46 <- function() {
  if (getRversion() < "4.6.0") {
    return(invisible(FALSE))
  }
  if (!requireNamespace("vscDebugger", quietly = TRUE)) {
    return(invisible(FALSE))
  }

  if (isTRUE(getOption("vsc.r46.promise_patch", FALSE))) {
    return(invisible(TRUE))
  }

  ns <- asNamespace("vscDebugger")
  if (bindingIsLocked("isPromise", ns)) {
    unlockBinding("isPromise", ns)
  }
  assign(
    "isPromise",
    function(name, env, strict = TRUE) {
      name_chr <- as.character(name)
      if (!exists(name_chr, envir = env, inherits = FALSE)) {
        return(FALSE)
      }
      # No public R-level promise probe yet; use `get()` in getVarInEnv instead.
      FALSE
    },
    envir = ns
  )
  options(vsc.r46.promise_patch = TRUE)

  invisible(TRUE)
}
