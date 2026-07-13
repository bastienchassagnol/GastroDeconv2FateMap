source("renv/activate.R")

# Ensure the VSCode-R extension temp directory exists at session start.
# The extension writes workspace.json here; the subdirectory is not created
# automatically when R recycles its tempdir between sessions.
local({
  vscode_tmp <- file.path(tempdir(), "vscode-R")
  if (!dir.exists(vscode_tmp)) dir.create(vscode_tmp, recursive = TRUE)
})

omnideconv_python <- path.expand(
  "~/.local/share/r-miniconda/envs/r-omnideconv/bin/python"
)
if (
  !nzchar(Sys.getenv("RETICULATE_PYTHON")) && file.exists(omnideconv_python)
) {
  Sys.setenv(RETICULATE_PYTHON = omnideconv_python)
}

# vscDebugger defaults for a smoother interactive experience
options(
  vsc.defaultDebugMode = "workspace",
  vsc.defaultAllowGlobalDebugging = TRUE,
  vsc.defaultIncludePackageScopes = TRUE,
  vsc.defaultOverwriteSource = TRUE,
  vsc.defaultOverwritePrint = TRUE,
  vsc.defaultOverwriteMessage = TRUE,
  vsc.setBreakpointsInStack = TRUE,
  vsc.showInternalFrames = FALSE,
  vsc.trySilent = TRUE,
  vsc.previewPromises = FALSE,
  vsc.showPromiseDetails = FALSE,
  vsc.evaluateActiveBindings = FALSE
)

# vscDebugger 0.5.6 C promise helpers segfault on R 4.6.0 (see GitHub issue #202).
local({
  project_root <- if (requireNamespace("renv", quietly = TRUE)) {
    renv::project()
  } else {
    getwd()
  }
  if (length(project_root) == 0L || !nzchar(project_root)) {
    project_root <- getwd()
  }
  patch_path <- file.path(project_root, "R", "vscdebugger_r46_patch.R")
  if (!file.exists(patch_path)) {
    return(invisible(NULL))
  }
  source(patch_path, local = TRUE)
  if (requireNamespace("vscDebugger", quietly = TRUE)) {
    patch_vscdebugger_promises_r46()
  }
  setHook(
    packageEvent("vscDebugger", "onLoad"),
    function(...) {
      patch_vscdebugger_promises_r46()
    }
  )
})

# VS Code / Cursor R session watcher (R 4.6+).
# - Do not source the Python .venv here; that is handled by uv for Python only.
# - Source the vscode-R init script so variables appear in the R workspace pane.
# - R 4.6 may not call `.First`; call `init_last()` explicitly after sourcing.
local({
  is_ide_terminal <- Sys.getenv("TERM_PROGRAM") %in% c("vscode", "cursor")
  if (!interactive() || !is_ide_terminal) {
    return(invisible(NULL))
  }

  # Skip attach for one-shot `R -e` / `Rscript` invocations (would hang).
  args <- commandArgs(trailingOnly = FALSE)
  if (any(grepl("^--file=", args)) || any(grepl("^-e$", args))) {
    return(invisible(NULL))
  }

  # vscode-R init.R only recognises TERM_PROGRAM == "vscode".
  if (identical(Sys.getenv("TERM_PROGRAM"), "cursor")) {
    Sys.setenv(TERM_PROGRAM = "vscode")
  }

  init_r <- path.expand("~/.vscode-R/init.R")
  if (!file.exists(init_r)) {
    return(invisible(NULL))
  }

  tryCatch(
    {
      source(init_r, local = FALSE)
      if (exists("init_last", mode = "function", inherits = TRUE)) {
        init_last()
      } else if (exists(".vsc.attach", mode = "function", inherits = TRUE)) {
        .vsc.attach()
      }
    },
    error = function(err) {
      message(
        "vscode-R session watcher could not attach: ",
        conditionMessage(err)
      )
    }
  )
  invisible(NULL)
})
