# Two-group power and size study for Luo-Huang WRS G*_rho.
# Run a task from the package root with WC_POWER_* variables set.

if (Sys.getenv("RUN_WC_POWER_SIMULATION") != "true") {
  stop("Set RUN_WC_POWER_SIMULATION=true to run this simulation.")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools first.")
}
devtools::load_all(".", quiet = TRUE)

n_sim <- as.integer(Sys.getenv("WC_POWER_N_SIM", "5000"))
n_tasks <- as.integer(Sys.getenv("WC_POWER_N_TASKS", "20"))
task_id <- as.integer(Sys.getenv("WC_POWER_TASK_ID", "1"))
seed <- as.integer(Sys.getenv("WC_POWER_SEED", "20260928"))
alpha <- as.numeric(Sys.getenv("WC_POWER_ALPHA", "0.05"))
output_dir <- Sys.getenv("WC_POWER_RESULTS_DIR",
                         file.path("data-raw", "wc-logrank-power", "results"))
if (anyNA(c(n_sim, n_tasks, task_id, seed, alpha)) || n_sim < 1L ||
    n_tasks < 1L || task_id < 1L || task_id > n_tasks ||
    alpha <= 0 || alpha >= 1) {
  stop("Invalid WC_POWER_* simulation setting.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Zhao-style recurrent histories provide a transparent, censored gap-time
# generator: subject-specific multiplicative heterogeneity, staggered entry,
# and administrative study end. The generator is used here as a data generator;
# the tests evaluated are the separate Luo-Huang G*_rho implementations.
simulate_gaps <- getFromNamespace(".zhao_simulated_gap_data", "recurSurvTests")
heterogeneity_ranges <- list(low = c(0.5, 1.5), high = c(0.01, 1.99))
grid <- expand.grid(
  distribution = c("exponential", "weibull"),
  heterogeneity = names(heterogeneity_ranges),
  n_per_group = c(25L, 50L, 100L, 200L),
  log_time_shift = c(-0.35, -0.20, 0, 0.20, 0.35),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$cell_id <- seq_len(nrow(grid))
cells <- grid[seq(task_id, nrow(grid), by = n_tasks), , drop = FALSE]

started <- proc.time()[["elapsed"]]
total_replicates <- nrow(cells) * n_sim
raw <- vector("list", total_replicates)
counter <- 0L
progress <- utils::txtProgressBar(min = 0, max = total_replicates, style = 3)
report_every <- max(1L, ceiling(total_replicates / 20L))
format_duration <- function(seconds) {
  seconds <- max(0, round(seconds))
  sprintf("%02d:%02d:%02d", seconds %/% 3600L,
          (seconds %% 3600L) %/% 60L, seconds %% 60L)
}
for (cell_index in seq_len(nrow(cells))) {
  cell <- cells[cell_index, ]
  for (replicate_id in seq_len(n_sim)) {
    counter <- counter + 1L
    set.seed(seed + cell$cell_id * 100000L + replicate_id)
    dat <- simulate_gaps(
      n_per_group = cell$n_per_group,
      distribution = cell$distribution,
      heterogeneity = heterogeneity_ranges[[cell$heterogeneity]],
      followup = 180,
      staggered_entry = TRUE,
      treatment_time_multiplier = exp(cell$log_time_shift)
    )
    p0 <- wc_logrank(dat, group = "group", episode = "episode", rho = 0)$p.value
    p1 <- wc_logrank(dat, group = "group", episode = "episode", rho = 1)$p.value
    raw[[counter]] <- data.frame(
      cell_id = cell$cell_id,
      distribution = cell$distribution,
      heterogeneity = cell$heterogeneity,
      n_per_group = cell$n_per_group,
      log_time_shift = cell$log_time_shift,
      multiplier = exp(cell$log_time_shift),
      replicate = replicate_id,
      p_rho0 = p0,
      p_rho1 = p1
    )
    utils::setTxtProgressBar(progress, counter)
    if (counter == 1L || counter %% report_every == 0L ||
        counter == total_replicates) {
      elapsed <- proc.time()[["elapsed"]] - started
      remaining <- elapsed / counter * (total_replicates - counter)
      message(sprintf(
        "\nTask %d/%d: %d/%d replicates; current cell %d/%d (%s, %s heterogeneity, n/arm=%d, shift=%+.2f); elapsed %s; ETA %s",
        task_id, n_tasks, counter, total_replicates, cell_index,
        nrow(cells), cell$distribution, cell$heterogeneity,
        cell$n_per_group, cell$log_time_shift,
        format_duration(elapsed), format_duration(remaining)
      ))
    }
  }
}
close(progress)
raw <- do.call(rbind, raw)
outfile <- file.path(output_dir, sprintf("wc-logrank-power-task-%02d-of-%02d.rds", task_id, n_tasks))
saveRDS(list(
  settings = list(n_sim_per_cell = n_sim, n_tasks = n_tasks, task_id = task_id,
    seed = seed, alpha = alpha, followup = 180, staggered_entry = TRUE,
    distributions = c("exponential", "weibull"),
    heterogeneity_ranges = heterogeneity_ranges,
    sample_sizes_per_group = c(25L, 50L, 100L, 200L),
    log_time_shifts = c(-0.35, -0.20, 0, 0.20, 0.35)),
  raw = raw
), outfile)
message("Saved ", nrow(raw), " simulation rows to ", outfile)
