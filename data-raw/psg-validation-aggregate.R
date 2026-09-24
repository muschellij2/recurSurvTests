# Combine task-level raw p-values from psg-validation-simulations.R.
# Run after every Slurm array task has completed:
# Rscript data-raw/psg-validation-aggregate.R

files <- sort(list.files("data-raw", "^psg-validation-results-task[0-9]+\\.rds$",
  full.names = TRUE))
if (!length(files)) stop("No task-level PSG result files found in data-raw/.")
task_result <- lapply(files, readRDS)
raw_result <- do.call(rbind, lapply(task_result, `[[`, "raw"))
key <- raw_result[c("design", "alternative", "replicate", "method")]
if (anyDuplicated(key)) stop("Duplicate scenario-replicate-method rows found.")
result <- stats::aggregate(p.value ~ design + alternative + method, raw_result,
  function(x) mean(x < .05))
names(result)[names(result) == "p.value"] <- "rejection_rate"
result$n_sim <- as.integer(stats::aggregate(p.value ~ design + alternative + method,
  raw_result, length)$p.value)
result$monte_carlo_se <- sqrt(result$rejection_rate * (1 - result$rejection_rate) / result$n_sim)
B_value <- unique(vapply(task_result, function(x) unique(x$summary$B), integer(1)))
if (length(B_value) != 1L) stop("Task-level results use different values of B.")
result$B <- B_value
saveRDS(list(summary = result, raw = raw_result), "data-raw/psg-validation-results.rds")
utils::write.csv(result, "data-raw/psg-validation-results.csv", row.names = FALSE)
grDevices::png("data-raw/psg-validation-rejection-rates.png", 1400, 800)
graphics::barplot(result$rejection_rate, names.arg = paste(result$design, result$alternative, result$method, sep = "\n"),
  las = 2, ylim = c(0, 1), ylab = "Rejection rate", main = "PSG-shaped recurrent gap-time simulation")
graphics::abline(h = .05, lty = 2, col = "firebrick")
grDevices::dev.off()
print(result)
