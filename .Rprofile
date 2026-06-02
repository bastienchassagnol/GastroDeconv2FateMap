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

# R 4.6 no longer calls `.First.sys()`. Ensure VS Code session watcher can
# still attach when an interactive terminal session starts.
# local({
#   if (identical(Sys.getenv("TERM_PROGRAM"), "vscode")) {
#     old_first <- get0(".First", envir = .GlobalEnv, mode = "function")
#     assign(
#       ".First",
#       function() {
#         if (is.function(old_first)) {
#           old_first()
#         }
#         if (exists(".vsc.attach", mode = "function", inherits = TRUE)) {
#           try(.vsc.attach(), silent = TRUE)
#         }
#       },
#       envir = .GlobalEnv
#     )
#   }
# })
