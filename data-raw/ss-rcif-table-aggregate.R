# Aggregate completed Sivadasan--Sankaran simulation cells.
# Run from the package root after one or more tasks finish:
# Rscript data-raw/ss-rcif-table-aggregate.R [results_dir] [output_csv]
#   [published_values_csv]
# Optional published CSV columns: table,n,time,cause,abs_bias,mse.

args <- commandArgs(trailingOnly = TRUE)
indir <- if (length(args)) args[1] else "data-raw/ss-rcif-tables/results"
outfile <- if (length(args) > 1L) args[2] else
  "data-raw/ss-rcif-tables/summary.csv"
files <- sort(list.files(indir, pattern = "^table-[0-9]+-n[0-9]+-reps[0-9]+[.]rds$",
                         full.names = TRUE))
if (!length(files)) stop("No table simulation results in ", indir)

runs <- lapply(files, readRDS)
keys <- vapply(runs, function(x) paste(x$metadata$table, x$metadata$n, sep = ":"),
               character(1))
if (anyDuplicated(keys)) stop("Multiple files for the same table/sample-size cell")
raw <- do.call(rbind, lapply(runs, `[[`, "raw"))
if (!nrow(raw)) stop("No successful estimates")

cell <- split(raw, interaction(raw$table, raw$n, raw$method, raw$cause,
                               raw$time, drop = TRUE))
summary <- do.call(rbind, lapply(cell, function(x) {
  e <- x$estimate - x$truth
  sq <- e^2
  bias <- mean(e)
  data.frame(
    table = x$table[1], n = x$n[1], censoring = x$censoring[1],
    censor_limit = x$censor_limit[1], subject_weight = x$subject_weight[1],
    method = x$method[1], cause = x$cause[1], time = x$time[1],
    n_success = length(e), truth = x$truth[1], mean_estimate = mean(x$estimate),
    bias = bias, abs_bias = abs(bias), mean_absolute_error = mean(abs(e)),
    mse = mean(sq), bias_mc_se = stats::sd(e) / sqrt(length(e)),
    mse_mc_se = stats::sd(sq) / sqrt(length(e)),
    mean_singleton_subjects = mean(x$singleton_subjects),
    mean_probability_conservation_error =
      mean(x$probability_conservation_max_abs)
  )
}))
summary <- summary[order(summary$table, summary$n, summary$method,
                         summary$cause, summary$time), , drop = FALSE]
rownames(summary) <- NULL
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, outfile, row.names = FALSE)

inventory <- data.frame(
  table = vapply(runs, function(x) x$metadata$table, integer(1)),
  n = vapply(runs, function(x) x$metadata$n, integer(1)),
  requested_reps = vapply(runs, function(x) x$metadata$requested_reps, integer(1)),
  failed_fits = vapply(runs, function(x) length(x$failures), integer(1)),
  output_file = files
)
utils::write.csv(inventory,
                 file.path(dirname(outfile), "run-inventory.csv"), row.names = FALSE)
if (length(args) > 2L) {
  published <- utils::read.csv(args[3], stringsAsFactors = FALSE)
  needed <- c("table", "n", "time", "cause", "abs_bias", "mse")
  if (!all(needed %in% names(published))) {
    stop("Published CSV needs columns: ", paste(needed, collapse = ","))
  }
  if (anyDuplicated(published[c("table", "n", "time", "cause")])) {
    stop("Duplicate published benchmark cells")
  }
  comparison <- merge(summary, published, by = c("table", "n", "time", "cause"),
                      suffixes = c("_local", "_published"))
  comparison$abs_bias_difference <- comparison$abs_bias_local -
    comparison$abs_bias_published
  comparison$mse_difference <- comparison$mse_local -
    comparison$mse_published
  utils::write.csv(comparison,
    file.path(dirname(outfile), "published-comparison.csv"), row.names = FALSE)
  message("Compared ", nrow(comparison), " local rows with published benchmarks")
}
message("Wrote ", nrow(summary), " summary rows from ", length(files),
        " of 24 planned cells to ", outfile)
