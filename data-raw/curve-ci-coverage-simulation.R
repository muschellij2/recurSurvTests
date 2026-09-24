# Coverage simulation for subject-bootstrap recurrent survival intervals.
#
# This script is intentionally opt-in because the default configuration uses a
# large Monte-Carlo reference population plus 1,000 data sets x 999 cluster
# bootstrap fits.  Run, for example:
#
#   RUN_CURVE_CI_COVERAGE=true Rscript data-raw/curve-ci-coverage-simulation.R
#
# The truth is the large-sample limit of each *implemented estimator* under the
# same recurrent-record design.  This is deliberate.  With terminal-gap data,
# treating all rows as IID exponential observations does not in general define
# the PSH target, so exp(-lambda t) would be an invalid coverage benchmark.

if (!identical(tolower(Sys.getenv("RUN_CURVE_CI_COVERAGE")), "true")) {
  message("Set RUN_CURVE_CI_COVERAGE=true to run this coverage simulation.")
} else {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(recurSurvTests)
  }

  simulate_histories <- function(n_subjects, completed_gaps = 5L,
                                 followup = 8, rate = 1) {
    out <- lapply(seq_len(n_subjects), function(i) {
      elapsed <- 0
      rows <- list()
      episode <- 0L
      while (elapsed < followup) {
        # PSG-like scored epochs make the event-time grid finite. This also
        # makes a large (for example, 100,000-subject) Monte-Carlo reference
        # computationally
        # feasible for the WC weighted-risk-set calculation.
        y <- ceiling(stats::rexp(1, rate = rate) / 0.25) * 0.25
        if (elapsed + y >= followup) {
          episode <- episode + 1L
          rows[[length(rows) + 1L]] <- data.frame(
            id = i, episode = episode, gap = followup - elapsed, status = 0
          )
          break
        }
        episode <- episode + 1L
        rows[[length(rows) + 1L]] <- data.frame(
          id = i, episode = episode, gap = y, status = 1
        )
        elapsed <- elapsed + y
        if (episode >= completed_gaps) {
          # A terminal administrative gap preserves the required final-censor
          # convention while yielding unequal recurrent-event counts.
          episode <- episode + 1L
          rows[[length(rows) + 1L]] <- data.frame(
            id = i, episode = episode, gap = followup - elapsed, status = 0
          )
          break
        }
      }
      do.call(rbind, rows)
    })
    do.call(rbind, out)
  }

  truth_n <- as.integer(Sys.getenv("CURVE_CI_TRUTH_N", "25000"))
  task_id <- as.integer(Sys.getenv("CURVE_CI_TASK_ID", "1"))
  n_tasks <- as.integer(Sys.getenv("CURVE_CI_N_TASKS", "1"))
  set.seed(20260918)
  started_at <- proc.time()[["elapsed"]]
  times <- c(0.5, 1, 2)
  truth_data <- simulate_histories(truth_n)
  truth <- lapply(c("wc", "psh"), function(estimator) {
    fit <- if (estimator == "wc") {
      wc_surv(truth_data, episode = "episode")
    } else {
      psh_surv(truth_data, episode = "episode")
    }
    recurSurvTests:::.recurrent_surv_step(fit, times)
  })
  names(truth) <- c("wc", "psh")
  truth_elapsed <- proc.time()[["elapsed"]] - started_at
  format_duration <- function(seconds) {
    seconds <- max(0, round(seconds))
    sprintf("%02d:%02d:%02d", seconds %/% 3600L,
            (seconds %% 3600L) %/% 60L, seconds %% 60L)
  }
  message(sprintf("Reference curves complete in %s.", format_duration(truth_elapsed)))

  one_replicate <- function(rep_id, n_subjects = 100L, B = 999L) {
    d <- simulate_histories(n_subjects)
    do.call(rbind, lapply(c("wc", "psh"), function(estimator) {
      do.call(rbind, lapply(c("pointwise", "simultaneous"), function(interval) {
        ci <- recurrent_surv_ci(
          d, estimator = estimator, episode = "episode", times = times,
          interval = interval, B = B
        )
        data.frame(
          replicate = rep_id, estimator = estimator, interval = interval,
          time = times, truth = truth[[estimator]],
          covered = truth[[estimator]] >= ci$curve$lower &
            truth[[estimator]] <= ci$curve$upper
        )
      }))
    }))
  }

  n_total <- as.integer(Sys.getenv("CURVE_CI_N_SIM", "1000"))
  n_per_task <- ceiling(n_total / n_tasks)
  n_sim <- max(0L, min(n_per_task, n_total - (task_id - 1L) * n_per_task))
  B <- as.integer(Sys.getenv("CURVE_CI_B", "999"))
  first_replicate <- (task_id - 1L) * n_per_task + 1L
  progress_bar <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  completed <- 0L
  results <- do.call(rbind, lapply(seq_len(n_sim), function(i) {
    value <- one_replicate(first_replicate + i - 1L, B = B)
    completed <<- completed + 1L
    utils::setTxtProgressBar(progress_bar, completed)
    report_every <- max(1L, ceiling(n_sim / 20L))
    if (completed == 1L || completed %% report_every == 0L || completed == n_sim) {
      simulation_elapsed <- proc.time()[["elapsed"]] - started_at - truth_elapsed
      remaining <- simulation_elapsed / completed * (n_sim - completed)
      message(sprintf("\nCurve-CI simulations: %d/%d; elapsed %s; ETA %s; projected task total %s",
        completed, n_sim, format_duration(simulation_elapsed), format_duration(remaining),
        format_duration(truth_elapsed + simulation_elapsed + remaining)))
    }
    value
  }))
  close(progress_bar)
  pointwise <- aggregate(covered ~ estimator + interval + time, results, mean)
  split_key <- interaction(results$replicate, results$estimator,
                           results$interval, drop = TRUE)
  joint <- do.call(rbind, lapply(split(results, split_key), function(x) {
    data.frame(replicate = x$replicate[1], estimator = x$estimator[1],
               interval = x$interval[1], covered = all(x$covered))
  }))
  simultaneous <- aggregate(covered ~ estimator + interval, joint, mean)
  names(pointwise)[names(pointwise) == "covered"] <- "coverage"
  names(simultaneous)[names(simultaneous) == "covered"] <- "joint_coverage"
  pointwise$mc_se <- sqrt(pointwise$coverage * (1 - pointwise$coverage) / n_sim)
  simultaneous$mc_se <- sqrt(
    simultaneous$joint_coverage * (1 - simultaneous$joint_coverage) / n_sim
  )

  suffix <- if (n_tasks > 1L) paste0("-task", task_id) else ""
  saveRDS(list(truth = truth, results = results, pointwise = pointwise,
               simultaneous = simultaneous),
          paste0("data-raw/curve-ci-coverage-results", suffix, ".rds"))
  utils::write.csv(pointwise, paste0("data-raw/curve-ci-pointwise-coverage", suffix, ".csv"),
                   row.names = FALSE)
  utils::write.csv(simultaneous, paste0("data-raw/curve-ci-simultaneous-coverage", suffix, ".csv"),
                   row.names = FALSE)
}
