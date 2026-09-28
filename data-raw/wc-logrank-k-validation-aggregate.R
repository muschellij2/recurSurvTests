# Combine completed K-group validation array tasks.
# Run from the package root:
#   Rscript data-raw/wc-logrank-k-validation-aggregate.R

result_dir <- file.path("data-raw", "wc-logrank-k-validation-results")
files <- sort(list.files(
  result_dir,
  pattern = "^wc-logrank-k-validation-results-task[0-9]+[.]rds$",
  full.names = TRUE
))
if (!length(files)) stop("No task-level K-group results found in ", result_dir)

task_id <- as.integer(sub(
  "^.*-task([0-9]+)[.]rds$", "\\1", basename(files)
))
if (anyNA(task_id) || anyDuplicated(task_id)) stop("Invalid or duplicate task IDs")
order_id <- order(task_id)
files <- files[order_id]
task_id <- task_id[order_id]
if (!identical(task_id, seq_len(max(task_id)))) {
  stop("Task files are incomplete; found IDs: ", paste(task_id, collapse = ", "))
}

task_results <- lapply(files, readRDS)
for (i in seq_along(task_results)) {
  x <- task_results[[i]]
  if (!is.data.frame(x) ||
      !all(c("scenario", "rejection_rate", "mc_se", "n_sim") %in% names(x)) ||
      any(!is.finite(x$rejection_rate)) || any(x$rejection_rate < 0 | x$rejection_rate > 1)) {
    stop("Malformed task result: ", basename(files[i]))
  }
}
task_table <- do.call(rbind, lapply(seq_along(task_results), function(i) {
  x <- task_results[[i]]
  x$task_id <- task_id[i]
  x
}))
if (anyDuplicated(task_table[c("task_id", "scenario")])) {
  stop("Duplicate task-scenario rows")
}

grouped <- split(task_table, task_table$scenario)
summary <- do.call(rbind, lapply(names(grouped), function(scenario) {
  x <- grouped[[scenario]]
  n <- sum(x$n_sim)
  rejected <- sum(x$rejection_rate * x$n_sim)
  if (abs(rejected - round(rejected)) > 1e-6) {
    stop("Task rejection rates do not imply integer rejection counts")
  }
  rejected <- as.integer(round(rejected))
  rate <- rejected / n
  ci <- stats::binom.test(rejected, n)$conf.int
  data.frame(
    scenario = scenario,
    rejection_rate = rate,
    mc_se = sqrt(rate * (1 - rate) / n),
    n_sim = n,
    rejected = rejected,
    ci_low = unname(ci[1]),
    ci_high = unname(ci[2]),
    null_p_value = if (scenario == "null") {
      stats::binom.test(rejected, n, p = 0.05)$p.value
    } else {
      NA_real_
    }
  )
}))
rownames(summary) <- NULL

saveRDS(list(summary = summary, task_results = task_table),
        file.path(result_dir, "wc-logrank-k-validation-results.rds"))
utils::write.csv(summary,
                  file.path(result_dir, "wc-logrank-k-validation-results.csv"),
                  row.names = FALSE)
utils::write.csv(task_table,
                  file.path(result_dir, "wc-logrank-k-validation-task-manifest.csv"),
                  row.names = FALSE)
print(summary)
