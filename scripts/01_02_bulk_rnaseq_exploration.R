# ==========================================================================
# 0. Output naming convention ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/bulk-omics"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
