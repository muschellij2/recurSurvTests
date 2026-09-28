# Aggregate the two-group WRS power simulation.
# Run from the package root after all tasks finish:
#   Rscript data-raw/wc-logrank-power/aggregate.R

result_dir <- file.path("data-raw", "wc-logrank-power", "results")
files <- sort(list.files(
  result_dir,
  pattern = "^wc-logrank-power-task-[0-9]+-of-[0-9]+[.]rds$",
  full.names = TRUE
))
if (!length(files)) stop("No task results found in ", result_dir)
objects <- lapply(files, readRDS)
settings <- objects[[1]]$settings
task_ids <- vapply(objects, function(x) x$settings$task_id, integer(1))
if (anyDuplicated(task_ids) || !identical(sort(task_ids), seq_len(settings$n_tasks))) {
  stop("Expected one output for each task 1:", settings$n_tasks)
}
for (x in objects) {
  if (!identical(x$settings$n_sim_per_cell, settings$n_sim_per_cell) ||
      !identical(x$settings$n_tasks, settings$n_tasks) ||
      !identical(x$settings$seed, settings$seed) ||
      !identical(x$settings$alpha, settings$alpha)) {
    stop("Task outputs have mixed simulation settings")
  }
}
raw <- do.call(rbind, lapply(objects, `[[`, "raw"))
key <- raw[c("cell_id", "replicate")]
if (anyDuplicated(key)) stop("Duplicate cell-replicate rows found")
expected_cells <- 80L
if (!setequal(unique(raw$cell_id), seq_len(expected_cells))) {
  stop("Expected all 80 design cells")
}
counts <- aggregate(
  cbind(reject_rho0 = as.integer(raw$p_rho0 < settings$alpha),
        reject_rho1 = as.integer(raw$p_rho1 < settings$alpha)) ~
    cell_id + distribution + heterogeneity + n_per_group + log_time_shift + multiplier,
  data = raw,
  FUN = sum
)
cell_n <- aggregate(p_rho0 ~ cell_id, raw, length)
counts$n_sim <- as.integer(cell_n$p_rho0[match(counts$cell_id, cell_n$cell_id)])
if (any(counts$n_sim != settings$n_sim_per_cell)) {
  stop("At least one cell has an unexpected replicate count")
}
for (method in c("rho0", "rho1")) {
  reject_col <- paste0("reject_", method)
  counts[[paste0("power_", method)]] <- counts[[reject_col]] / counts$n_sim
  counts[[paste0("mc_se_", method)]] <- sqrt(
    counts[[paste0("power_", method)]] * (1 - counts[[paste0("power_", method)]]) / counts$n_sim
  )
  ci <- t(vapply(seq_len(nrow(counts)), function(i) {
    stats::binom.test(counts[[reject_col]][i], counts$n_sim[i])$conf.int
  }, numeric(2)))
  counts[[paste0("ci_low_", method)]] <- ci[, 1]
  counts[[paste0("ci_high_", method)]] <- ci[, 2]
}
counts$alpha <- settings$alpha
counts <- counts[order(counts$distribution, counts$heterogeneity,
                       counts$n_per_group, counts$log_time_shift), ]
rownames(counts) <- NULL

saveRDS(list(settings = settings, summary = counts, raw = raw),
        file.path("data-raw", "wc-logrank-power", "results.rds"))
utils::write.csv(counts,
  file.path("data-raw", "wc-logrank-power", "results.csv"), row.names = FALSE)

grDevices::png(file.path("data-raw", "wc-logrank-power", "power-curves.png"),
               width = 1800, height = 1400, res = 160)
old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
cols <- c("25" = "#0072B2", "50" = "#009E73", "100" = "#D55E00", "200" = "#CC79A7")
for (dist in settings$distributions) for (het in names(settings$heterogeneity_ranges)) {
  d <- counts[counts$distribution == dist & counts$heterogeneity == het, ]
  graphics::plot(NA, xlim = range(d$multiplier), ylim = c(0, 1),
    xlab = "Treatment gap-time multiplier", ylab = "Rejection probability",
    main = paste(dist, "baseline;", het, "frailty"))
  graphics::abline(h = settings$alpha, col = "grey50", lty = 3)
  for (n in sort(unique(d$n_per_group))) for (method in c("rho0", "rho1")) {
    q <- d[d$n_per_group == n, ]
    q <- q[order(q$multiplier), ]
    power_col <- paste0("power_", method)
    graphics::lines(q$multiplier, q[[power_col]], col = cols[as.character(n)],
                    lty = if (method == "rho0") 1 else 2, lwd = 1.4)
    graphics::points(q$multiplier, q[[power_col]], col = cols[as.character(n)],
                     pch = if (method == "rho0") 16 else 1, cex = .8)
  }
}
graphics::par(old_par)
grDevices::dev.off()
print(counts)
