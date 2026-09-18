# PSG-shaped validation simulations. Run explicitly from the package root:
# RUN_PSG_VALIDATION=true Rscript data-raw/psg-validation-simulations.R
if (Sys.getenv("RUN_PSG_VALIDATION") != "true") stop("Set RUN_PSG_VALIDATION=true to run.")
if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools first.")
devtools::load_all(".", quiet = TRUE)
n_total <- as.integer(Sys.getenv("PSG_N_SIM", "1000"))
task_id <- as.integer(Sys.getenv("PSG_TASK_ID", "1")); n_tasks <- as.integer(Sys.getenv("PSG_N_TASKS", "1"))
n_sim <- ceiling(n_total / n_tasks)
B <- as.integer(Sys.getenv("PSG_B", "999"))
set.seed(as.integer(Sys.getenv("PSG_SEED", "20260922")) + task_id - 1L)

# A visit has shared frailty, a visit effect, persistent burst/recovery states,
# unequal 360--480 minute recording duration, and a terminal censored gap.
simulate_visit <- function(id, pair, arm, frailty, multiplier = 1) {
  recording <- stats::runif(1, 360, 480); visit_effect <- stats::rlnorm(1, 0, .18)
  state <- stats::rbinom(1, 1, .5); total <- 0; gap <- numeric()
  repeat {
    state_scale <- if (state == 1) .55 else 1.8
    next_gap <- stats::rlnorm(1, log(2.4 * frailty * visit_effect * multiplier * state_scale), .4)
    if (total + next_gap >= recording) { gap <- c(gap, recording - total); break }
    gap <- c(gap, next_gap); total <- total + next_gap
    if (stats::runif(1) > .85) state <- 1 - state
  }
  data.frame(id = id, pair = pair, arm = arm, episode = seq_along(gap), gap = gap,
    status = c(rep(1L, length(gap) - 1L), 0L))
}

simulate_psg <- function(paired = FALSE, multiplier = 1, n_control = 30, n_treatment = 45) {
  if (!paired) {
    arm <- c(rep("control", n_control), rep("treatment", n_treatment))
    frailty <- stats::rlnorm(length(arm), 0, .55)
    return(do.call(rbind, lapply(seq_along(arm), function(i) simulate_visit(
      i, NA_integer_, arm[i], frailty[i], if (arm[i] == "treatment") multiplier else 1))))
  }
  frailty <- stats::rlnorm(n_control, 0, .55)
  do.call(rbind, unlist(lapply(seq_len(n_control), function(i) list(
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
result <- do.call(rbind, lapply(seq_len(nrow(scenario)), function(i) {
  s <- scenario[i, ]; p <- replicate(n_sim, one_run(s$design == "paired_visits", s$multiplier))
  rate <- rowMeans(p < .05)
  data.frame(design = s$design, alternative = s$alternative, method = names(rate),
    rejection_rate = unname(rate), monte_carlo_se = sqrt(rate * (1 - rate) / n_sim), n_sim = n_sim, B = B)
}))
suffix <- if (n_tasks > 1L) paste0("-task", task_id) else ""
saveRDS(result, paste0("data-raw/psg-validation-results", suffix, ".rds"))
utils::write.csv(result, paste0("data-raw/psg-validation-results", suffix, ".csv"), row.names = FALSE)
grDevices::png(paste0("data-raw/psg-validation-rejection-rates", suffix, ".png"), 1400, 800)
graphics::barplot(result$rejection_rate, names.arg = paste(result$design, result$alternative, result$method, sep = "\n"),
  las = 2, ylim = c(0, 1), ylab = "Rejection rate", main = "PSG-shaped recurrent gap-time simulation")
graphics::abline(h = .05, lty = 2, col = "firebrick")
grDevices::dev.off()
print(result)
