# PSG-shaped validation simulations. Run explicitly from the package root:
# RUN_PSG_VALIDATION=true Rscript data-raw/psg-validation-simulations.R
if (Sys.getenv("RUN_PSG_VALIDATION") != "true") stop("Set RUN_PSG_VALIDATION=true to run.")
n_total <- as.integer(Sys.getenv("PSG_N_SIM", "1000"))
task_id <- as.integer(Sys.getenv("PSG_TASK_ID", "1")); n_tasks <- as.integer(Sys.getenv("PSG_N_TASKS", "1"))
B <- as.integer(Sys.getenv("PSG_B", "999"))
seed <- as.integer(Sys.getenv("PSG_SEED", "20260922"))
suffix <- if (n_tasks > 1L) paste0("-task", task_id) else ""
output_file <- paste0("data-raw/psg-validation-results", suffix, ".rds")

# Task-level RDS files are the canonical outputs.  This makes interrupted
# arrays safely resumable: a scheduler retry does not replace completed work.
if (file.exists(output_file)) {
  message("PSG validation output already exists: ", output_file, "; skipping.")
  quit(status = 0L)
}

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools first.")
devtools::load_all(".", quiet = TRUE)

format_duration <- function(seconds) {
  seconds <- max(0, round(seconds))
  sprintf("%02d:%02d:%02d", seconds %/% 3600L,
          (seconds %% 3600L) %/% 60L, seconds %% 60L)
}

# A visit has shared frailty, a visit effect, persistent burst/recovery states,
# an empirically calibrated overnight recording duration, and a terminal
# censored gap.  Newcastle PSG records had 420--801 minutes (median 565) and
# 30--199 contiguous scored-stage bouts (median 84).  The parameters below
# target that scale without requiring the external dataset at run time.
simulate_visit <- function(id, pair, arm, frailty, multiplier = 1) {
  recording <- pmin(800, pmax(420, stats::rlnorm(1, log(565), .23)))
  visit_effect <- stats::rlnorm(1, 0, .18)
  state <- stats::rbinom(1, 1, .5); total <- 0; gap <- numeric()
  repeat {
    state_scale <- if (state == 1) .55 else 1.8
    next_gap <- stats::rlnorm(1, log(3.6 * frailty * visit_effect * multiplier * state_scale), 1.1)
    if (total + next_gap >= recording) { gap <- c(gap, recording - total); break }
    gap <- c(gap, next_gap); total <- total + next_gap
    if (stats::runif(1) > .85) state <- 1 - state
  }
  data.frame(id = id, pair = pair, arm = arm, episode = seq_along(gap), gap = gap,
    status = c(rep(1L, length(gap) - 1L), 0L))
}

simulate_psg <- function(paired = FALSE, multiplier = 1, n_control = 7,
                         n_treatment = 21, n_pairs = 14) {
  if (!paired) {
    arm <- c(rep("control", n_control), rep("treatment", n_treatment))
    frailty <- stats::rlnorm(length(arm), 0, .55)
    return(do.call(rbind, lapply(seq_along(arm), function(i) simulate_visit(
      i, NA_integer_, arm[i], frailty[i], if (arm[i] == "treatment") multiplier else 1))))
  }
  # Newcastle has one recording per participant, so this paired design remains
  # a hypothetical within-person sensitivity analysis.  Fourteen pairs retain
  # the observed 28-visit scale.
  frailty <- stats::rlnorm(n_pairs, 0, .55)
  do.call(rbind, unlist(lapply(seq_len(n_pairs), function(i) list(
    simulate_visit(2L * i - 1L, i, "control", frailty[i]),
    simulate_visit(2L * i, i, "treatment", frailty[i], multiplier)
  )), recursive = FALSE))
}

one_run <- function(paired, multiplier) {
  d <- simulate_psg(paired, multiplier)
  if (!paired) return(c(
    wc_rho0 = wc_logrank(d, group = "arm", episode = "episode", rho = 0)$p.value,
    wc_rho1 = wc_logrank(d, group = "arm", episode = "episode", rho = 1)$p.value,
    zhao = zhao_rank_test(d, group = "arm", episode = "episode")$p.value
  ))
  c(
    paired_perm_rho0 = wc_paired_permutation(d, "pair", "arm", episode = "episode", rho = 0, B = B)$p.value,
    paired_perm_rho1 = wc_paired_permutation(d, "pair", "arm", episode = "episode", rho = 1, B = B)$p.value,
    zhao_unpaired_sensitivity = zhao_rank_test(d, group = "arm", episode = "episode")$p.value
  )
}

scenario <- expand.grid(design = c("independent_unequal_arms", "paired_visits"),
  alternative = c("null", "longer_bouts", "shorter_bouts"), stringsAsFactors = FALSE)
scenario$multiplier <- c(1, 1.25, .8)[match(scenario$alternative, c("null", "longer_bouts", "shorter_bouts"))]
total_work <- n_total * nrow(scenario)
work_per_task <- ceiling(total_work / n_tasks)
first_work <- (task_id - 1L) * work_per_task + 1L
work_id <- if (first_work > total_work) integer() else {
  seq.int(first_work, min(total_work, first_work + work_per_task - 1L))
}
n_sim <- length(work_id)
progress_bar <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
started_at <- proc.time()[["elapsed"]]
report_progress <- function(done) {
  utils::setTxtProgressBar(progress_bar, done)
  report_every <- max(1L, ceiling(n_sim / 20L))
  if (done == 1L || done %% report_every == 0L || done == n_sim) {
    elapsed <- proc.time()[["elapsed"]] - started_at
    remaining <- elapsed / done * (n_sim - done)
    message(sprintf("\nPSG scenario-replicates: %d/%d; elapsed %s; ETA %s; projected total %s",
      done, n_sim, format_duration(elapsed), format_duration(remaining),
      format_duration(elapsed + remaining)))
  }
}
raw_result <- do.call(rbind, lapply(seq_along(work_id), function(i) {
  cell <- work_id[i]
  # Interleave designs and alternatives so every array chunk has a comparable
  # mix of the inexpensive independent and costly paired-permutation work.
  scenario_id <- (cell - 1L) %% nrow(scenario) + 1L
  replicate_id <- (cell - 1L) %/% nrow(scenario) + 1L
  s <- scenario[scenario_id, ]
  set.seed(seed + cell - 1L)
  value <- one_run(s$design == "paired_visits", s$multiplier)
  report_progress(i)
  data.frame(design = s$design, alternative = s$alternative,
    replicate = replicate_id, method = names(value), p.value = unname(value))
}))
close(progress_bar)
result <- stats::aggregate(p.value ~ design + alternative + method, raw_result,
  function(x) mean(x < .05))
names(result)[names(result) == "p.value"] <- "rejection_rate"
result$n_sim <- as.integer(stats::aggregate(p.value ~ design + alternative + method,
  raw_result, length)$p.value)
result$monte_carlo_se <- sqrt(result$rejection_rate * (1 - result$rejection_rate) / result$n_sim)
result$B <- B
saveRDS(list(summary = result, raw = raw_result), output_file)
if (n_tasks == 1L) {
  utils::write.csv(result, "data-raw/psg-validation-results.csv", row.names = FALSE)
  grDevices::png("data-raw/psg-validation-rejection-rates.png", 1400, 800)
  graphics::barplot(result$rejection_rate, names.arg = paste(result$design, result$alternative, result$method, sep = "\n"),
    las = 2, ylim = c(0, 1), ylab = "Rejection rate", main = "PSG-shaped recurrent gap-time simulation")
  graphics::abline(h = .05, lty = 2, col = "firebrick")
  grDevices::dev.off()
}
print(result)
