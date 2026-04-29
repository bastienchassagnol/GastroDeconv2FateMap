source("renv/activate.R")

omnideconv_python <- path.expand(
  "~/.local/share/r-miniconda/envs/r-omnideconv/bin/python"
)
if (
  !nzchar(Sys.getenv("RETICULATE_PYTHON")) && file.exists(omnideconv_python)
) {
  Sys.setenv(RETICULATE_PYTHON = omnideconv_python)
}
