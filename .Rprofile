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
