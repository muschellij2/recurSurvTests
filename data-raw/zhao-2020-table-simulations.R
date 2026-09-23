# Zhao et al. (2020) table-reproduction simulations
#
# This is intentionally NOT run during package tests or builds. It is a
# long-running diagnostic: Table 2 alone has 12 cells x 3 tests x 100,000
# simulated datasets. Run explicitly, from the package root, with
#
#   RUN_ZHAO_100K=true Rscript data-raw/zhao-2020-table-simulations.R
#
# Runtime controls:
#
# * ZHAO_VARIANCE_METHOD: pooled_risk (default) or zhao_eq6.
# * ZHAO_PARALLEL: serial (default) or multicore.
# * ZHAO_N_CORES: local workers for multicore mode (default: 1).
# * ZHAO_TASK_ID and ZHAO_N_TASKS: one-based external array-task partition.
#   Each task writes its own RDS file, so this works with Slurm/PBS arrays or
#   separate machines. It can be combined with ZHAO_PARALLEL=multicore.
# * ZHAO_CHUNK_SIZE: replicates per progress-bar unit (default: 100).
# * ZHAO_N_NULL and ZHAO_N_POWER: optional smaller dry-run counts.

if (Sys.getenv("RUN_ZHAO_100K") != "true") {
  stop("Refusing to start 100,000-replicate simulations. Set RUN_ZHAO_100K=true to run.")
}


n_null <- as.integer(Sys.getenv("ZHAO_N_NULL", "100000"))
n_power <- as.integer(Sys.getenv("ZHAO_N_POWER", "10000"))
n_per_group <- 100L
followup <- 180
variance_method <- Sys.getenv("ZHAO_VARIANCE_METHOD", "pooled_risk")
parallel_mode <- match.arg(Sys.getenv("ZHAO_PARALLEL", "serial"), c("serial", "multicore"))
n_cores <- as.integer(Sys.getenv("ZHAO_N_CORES", "1"))
task_id <- as.integer(Sys.getenv("ZHAO_TASK_ID", "1"))
n_tasks <- as.integer(Sys.getenv("ZHAO_N_TASKS", "1"))
chunk_size <- as.integer(Sys.getenv("ZHAO_CHUNK_SIZE", "100"))

suffix <- if (n_tasks == 1L) "" else sprintf("-task-%03d-of-%03d", task_id, n_tasks)
outfile = paste0("data-raw/zhao-2020-table-simulations/zhao-2020-table-simulations", suffix, ".rds")
if (!file.exists(outfile)) {

  if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools to run this script.")
  devtools::load_all(".", quiet = TRUE)

  if (anyNA(c(n_null, n_power, n_cores, task_id, n_tasks, chunk_size)) ||
      n_null < 1L || n_power < 1L || n_cores < 1L || n_tasks < 1L ||
      task_id < 1L || task_id > n_tasks || chunk_size < 1L) {
    stop("Invalid ZHAO_* runtime control.")
  }
  if (parallel_mode == "multicore" && .Platform$OS.type == "windows") {
    stop("multicore mode requires a Unix-like platform; use external array tasks on Windows.")
  }
  distributions <- c("exponential", "weibull", "lognormal", "loglogistic")
  heterogeneity <- list(low = c(0.5, 1.5), medium = c(0.1, 1.9), high = c(0.01, 1.99))
  tests <- c("logrank", "peto_prentice", "gehan_breslow")

  run_chunk <- function(job, n_sim, seed) {
    z_range <- heterogeneity[[job$heterogeneity]]
    ans <- zhao_rank_simulation(
      n_sim = n_sim, n_per_group = n_per_group, scenario = "paper",
      distribution = job$distribution, heterogeneity = z_range, followup = followup,
      alpha = 0.025, test = job$test, variance_method = variance_method,
      alternative = "less", treatment_time_multiplier = job$multiplier, seed = seed
    )
    list(n = n_sim, rejection = sum(ans$p.value < 0.025, na.rm = TRUE),
         z = ans$z, mean_completed_events = ans$mean_completed_events)
  }

  combine_chunks <- function(job, chunks) {
    z <- unlist(lapply(chunks, `[[`, "z"), use.names = FALSE)
    n <- length(z)
    rejection <- sum(vapply(chunks, `[[`, numeric(1), "rejection"))
    event_mean <- weighted.mean(vapply(chunks, `[[`, numeric(1), "mean_completed_events"),
                                vapply(chunks, `[[`, numeric(1), "n"))
    data.frame(
      table = job$table, distribution = job$distribution,
      heterogeneity = job$heterogeneity, z_low = heterogeneity[[job$heterogeneity]][1],
      z_high = heterogeneity[[job$heterogeneity]][2], test = job$test,
      n_sim = n, treatment_time_multiplier = job$multiplier,
      empirical_rejection = rejection / n,
      monte_carlo_se = sqrt((rejection / n) * (1 - rejection / n) / n),
      mean_z = mean(z), sd_z = stats::sd(z), mean_completed_events = event_mean
    )
  }

  run_jobs <- function(jobs) {
    jobs <- jobs[seq(task_id, nrow(jobs), by = n_tasks), , drop = FALSE]
    chunks <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) {
      n <- jobs$n_sim[i]
      data.frame(job = i, n_sim = pmin(chunk_size, n - seq(0, n - 1, by = chunk_size)),
                 chunk = seq_len(ceiling(n / chunk_size)))
    }))
    progress <- utils::txtProgressBar(min = 0, max = nrow(chunks), style = 3)
    on.exit(close(progress), add = TRUE)
    output <- vector("list", nrow(jobs))
    for (i in seq_len(nrow(jobs))) output[[i]] <- list()
    # Run in waves so the progress bar advances after each completed multicore
    # batch, rather than only when an entire table cell has completed.
    for (start in seq(1L, nrow(chunks), by = n_cores)) {
      ii <- start:min(start + n_cores - 1L, nrow(chunks))
      work <- lapply(ii, function(k) {
        row <- chunks[k, ]
        list(job = row$job, chunk = row$chunk, n_sim = row$n_sim,
             seed = 1000000L * task_id + 10000L * row$job + row$chunk)
      })
      do_work <- function(x) {
        x$result <- run_chunk(jobs[x$job, ], x$n_sim, x$seed)
        x
      }
      if (parallel_mode == "multicore" && length(work) > 1L) {
        work <- parallel::mclapply(work, do_work, mc.cores = min(n_cores, length(work)))
      } else {
        work <- lapply(work, do_work)
      }
      for (x in work) output[[x$job]][[x$chunk]] <- x$result
      utils::setTxtProgressBar(progress, max(ii))
    }
    do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) combine_chunks(jobs[i, ], output[[i]])))
  }

  # Zhao Table 1: distribution of standardized null statistics at medium
  # heterogeneity. The paper reports 100,000 replicates and the empirical SD.
  table1_jobs <- expand.grid(distribution = distributions, test = tests,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  table1_jobs$table <- "table1"; table1_jobs$heterogeneity <- "medium"
  table1_jobs$multiplier <- 1; table1_jobs$n_sim <- n_null

  # Zhao Table 2: one-sided 0.025 null rejection rates, 100,000 replicates.
  table2_jobs <- expand.grid(distribution = distributions, heterogeneity = names(heterogeneity),
                             test = tests, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  table2_jobs$table <- "table2"; table2_jobs$multiplier <- 1; table2_jobs$n_sim <- n_null

  # Zhao Table 3: 10,000 one-sided power replicates. Each quoted alternative
  # changes a distribution location/scale by 0.25 or 0.5 on the log-time scale.
  table3_jobs <- expand.grid(distribution = distributions, heterogeneity = names(heterogeneity),
                             log_time_shift = c(0.25, 0.5), test = tests,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  table3_jobs$table <- "table3"; table3_jobs$multiplier <- exp(table3_jobs$log_time_shift)
  table3_jobs$n_sim <- n_power

  all_jobs <- rbind(table1_jobs[, c("table", "distribution", "heterogeneity", "test", "multiplier", "n_sim")],
                    table2_jobs[, c("table", "distribution", "heterogeneity", "test", "multiplier", "n_sim")],
                    table3_jobs[, c("table", "distribution", "heterogeneity", "test", "multiplier", "n_sim")])
  all_results <- run_jobs(all_jobs)

  result <- list(
    source = "Zhao et al. (2020), Tables 1--3",
    settings = list(n_null = n_null, n_power = n_power, n_per_group = n_per_group,
                    followup = followup, variance_method = variance_method,
                    parallel_mode = parallel_mode, n_cores = n_cores, task_id = task_id, n_tasks = n_tasks,
                    chunk_size = chunk_size,
                    note = "The paper used survsim::rec.ev.surv(); this package uses its own transparent gap-time generator, so exact equality is not expected."),
    table1 = all_results[all_results$table == "table1", , drop = FALSE],
    table2 = all_results[all_results$table == "table2", , drop = FALSE],
    table3 = all_results[all_results$table == "table3", , drop = FALSE]
  )

  suffix <- if (n_tasks == 1L) "" else sprintf("-task-%03d-of-%03d", task_id, n_tasks)
  saveRDS(result, paste0("data-raw/zhao-2020-table-simulations/zhao-2020-table-simulations", suffix, ".rds"))
  print(result$table1)
  print(result$table2)
  print(result$table3)
}
