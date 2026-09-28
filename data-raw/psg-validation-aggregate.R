# Combine task-level raw p-values from psg-validation-simulations.R.
# Run after every Slurm array task has completed:
# Rscript data-raw/psg-validation-aggregate.R
# For Equation 6 results, use the same switch used by the simulation array:
# PSG_ZHAO_VARIANCE_METHOD=zhao_eq6 Rscript data-raw/psg-validation-aggregate.R
library(tidyverse)

variance_method <- match.arg(
  Sys.getenv("PSG_ZHAO_VARIANCE_METHOD", "pooled_risk"),
  c("pooled_risk", "zhao_eq6")
)
results_dir <- if (variance_method == "pooled_risk") {
  file.path("data-raw", "psg-validation-results")
} else {
  file.path("data-raw", "psg-validation-results-zhao_eq6")
}
summary_dir <- if (variance_method == "pooled_risk") "data-raw" else results_dir
files <- sort(
  list.files(
    results_dir,
    pattern = "^psg-validation-results-task[0-9]+\\.rds$",
    full.names = TRUE
  )
)
if (!length(files))
  stop("No task-level PSG result files found in data-raw/.")
task_result <- map(files, readRDS, .progress = TRUE)
task_methods <- vapply(task_result, function(x) {
  if (is.null(x$variance_method)) "pooled_risk" else x$variance_method
}, character(1))
if (any(task_methods != variance_method))
  stop("Task-level results do not match PSG_ZHAO_VARIANCE_METHOD.")
raw_result <- map_df(task_result, `[[`, "raw")
key <- raw_result[c("design", "alternative", "replicate", "method")]
if (anyDuplicated(key))
  stop("Duplicate scenario-replicate-method rows found.")
result <- stats::aggregate(p.value ~ design + alternative + method, raw_result, function(x)
  mean(x < .05))
names(result)[names(result) == "p.value"] <- "rejection_rate"
result$n_sim <- as.integer(stats::aggregate(p.value ~ design + alternative + method, raw_result, length)$p.value)
result$monte_carlo_se <- sqrt(result$rejection_rate * (1 - result$rejection_rate) / result$n_sim)
B_value <- unique(vapply(task_result, function(x)
  unique(x$summary$B), integer(1)))
if (length(B_value) != 1L)
  stop("Task-level results use different values of B.")
result$B <- B_value
saveRDS(list(summary = result, raw = raw_result,
             variance_method = variance_method),
        file.path(summary_dir, "psg-validation-results.rds"))
utils::write.csv(result, file.path(summary_dir, "psg-validation-results.csv"), row.names = FALSE)
grDevices::png(file.path(summary_dir, "psg-validation-rejection-rates.png"), 1400, 800)
graphics::barplot(
  result$rejection_rate,
  names.arg = paste(result$design, result$alternative, result$method, sep = "\n"),
  las = 2,
  ylim = c(0, 1),
  ylab = "Rejection rate",
  main = "PSG-shaped recurrent gap-time simulation"
)
graphics::abline(h = .05,
                 lty = 2,
                 col = "firebrick")
grDevices::dev.off()
print(result)
