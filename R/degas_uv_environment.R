# ============================================================================
# DEGAS: lightweight uv-managed Python environment ----
# ============================================================================

#' Create or reuse a `uv`-managed Python environment for DEGAS
#'
#' @description
#' DEGAS (Johnson et al. 2022) requires an old TensorFlow build
#' (`tensorflow==2.4.1`, Python <= 3.8) that is incompatible with the
#' repository's main `uv` project (Python >= 3.13, see `pyproject.toml`).
#' A disposable virtual environment is created with
#' `uv venv --no-project` so the root `pyproject.toml` is not applied, then
#' packages are installed with `uv pip install --link-mode=copy`.
#'
#' @param project_dir Character scalar. Directory holding the dedicated `.venv`
#'   (e.g. `.uv/degas`); created if missing.
#' @param python_version Character scalar. Python version passed to
#'   `uv venv --python`. Default `"3.8.18"`.
#' @param packages Character vector of pinned `pip`-style requirements.
#' @param force Logical. Recreate the environment even if the interpreter
#'   already exists. Default `FALSE`.
#' @param verbose Logical. Print `uv` command output. Default `TRUE`.
#'
#' @return Character scalar: path to the environment's Python executable,
#'   suitable for `DEGAS.pyloc` (preferred over `reticulate::use_python()` when
#'   `RETICULATE_PYTHON` is set in the shell).
#'
#' @seealso `SigBridgeR::DoDEGAS()`, `SigBridgeR::Screen()`
#'
#' @importFrom glue glue
#' @export
create_degas_uv_env <- function(
    project_dir = file.path(".uv", "degas"),
    python_version = "3.8.18",
    packages = c(
      "tensorflow==2.4.1",
      "protobuf==3.20.3",
      "numpy==1.19.5"
    ),
    force = FALSE,
    verbose = TRUE
) {
  uv_bin <- Sys.which("uv")
  if (!nzchar(uv_bin)) {
    stop(
      "`uv` is not installed or is not on PATH; install uv rather than ",
      "falling back to Conda for the DEGAS Python dependencies.",
      call. = FALSE
    )
  }

  venv_dir <- file.path(project_dir, ".venv")
  python_bin <- file.path(
    venv_dir,
    if (.Platform$OS.type == "windows") "Scripts" else "bin",
    if (.Platform$OS.type == "windows") "python.exe" else "python"
  )

  if (force && dir.exists(venv_dir)) {
    unlink(venv_dir, recursive = TRUE, force = TRUE)
  }

  if (!file.exists(python_bin)) {
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)

    # `--no-project` avoids inheriting the repo's Python >= 3.13 constraint.
    .run_uv_checked(
      uv_bin,
      c("venv", "--no-project", "--python", python_version, venv_dir),
      verbose = verbose
    )
    # `--link-mode=copy` avoids cross-filesystem hardlink warnings on `/mnt`.
    .run_uv_checked(
      uv_bin,
      c(
        "pip", "install",
        "--no-project",
        "--link-mode=copy",
        "--python", python_bin,
        packages
      ),
      verbose = verbose
    )
  } else if (verbose) {
    message(glue::glue("Reusing existing DEGAS uv environment at {venv_dir}"))
  }

  stopifnot(file.exists(python_bin))
  python_bin
}

#' Run a `uv` subcommand and stop on non-zero exit status
#'
#' @param uv_bin Character scalar. Path to the `uv` executable.
#' @param args Character vector of arguments passed to `uv`.
#' @param verbose Logical. Print captured output on success.
#'
#' @return The captured command output, invisibly.
#'
#' @keywords internal
.run_uv_checked <- function(uv_bin, args, verbose = TRUE) {
  output <- system2(uv_bin, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  status <- if (is.null(status)) 0L else as.integer(status)

  if (!identical(status, 0L)) {
    stop(
      "Command failed: ", uv_bin, " ", paste(args, collapse = " "), "\n",
      paste(output, collapse = "\n"),
      call. = FALSE
    )
  }
  if (verbose) {
    message(paste(output, collapse = "\n"))
  }
  invisible(output)
}
