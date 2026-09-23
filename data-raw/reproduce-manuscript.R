#!/usr/bin/env Rscript

# Reproducibility driver for the JSS manuscript.
#
# This script is deliberately outside the package's runtime API. It renders the
# manuscript, imports any completed simulation RDS files, and writes public
# reader-facing tables plus internal execution metadata. It never launches the
# long-running simulation arrays; submit those separately with the documented
# sbatch files. The script is intended for a supplement, an attachment, or the
# companion GitHub repository rather than for routine package use.

repo_root <- normalizePath(Sys.getenv("RECURSURVTESTS_ROOT", getwd()), mustWork = TRUE)
manuscript <- file.path(repo_root, "vignettes", "jss-manuscript.qmd")
if (!file.exists(manuscript)) {
  stop("Run this script from the repository root or set RECURSURVTESTS_ROOT")
}

required <- c("quarto", "Rscript")
if (!nzchar(Sys.which("quarto"))) {
  stop("Quarto is required to render the manuscript")
}

if (!requireNamespace("here", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE) ||
    !requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install the manuscript dependencies: here, knitr, and pkgload")
}

output_dir <- Sys.getenv("MANUSCRIPT_OUTPUT_DIR", file.path(repo_root, "docs", "jss-manuscript"))
if (!grepl("^(/|[A-Za-z]:)", output_dir)) output_dir <- file.path(repo_root, output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Rendering ", manuscript)
status <- system2("quarto", c("render", manuscript), stdout = "", stderr = "")
if (!identical(status, 0L)) {
  stop("Quarto rendering failed with status ", status)
}

message("Manuscript artifacts written to ", output_dir)
message("Internal status files are under ", file.path(output_dir, "_internal"))
