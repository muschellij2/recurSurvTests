# Zhao et al. (2020) table-reproduction simulations
#
# This is intentionally NOT run during package tests or builds. It is a
# long-running diagnostic: Table 2 alone has 12 cells x 3 tests x 100,000
# simulated datasets. Run explicitly, from the package root, with
#
#   RUN_ZHAO_100K=true Rscript data-raw/zhao-2020-table-simulations.R
#
# Set ZHAO_VARIANCE_METHOD to zhao_eq6 to evaluate the literal Equation 6 path,
# or leave its calibrated pooled_risk default.

if (Sys.getenv("RUN_ZHAO_100K") != "true") {
  stop("Refusing to start 100,000-replicate simulations. Set RUN_ZHAO_100K=true to run.")
}

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools to run this script.")
devtools::load_all(".", quiet = TRUE)

n_null <- 100000L
n_power <- 10000L
n_per_group <- 100L
followup <- 180
variance_method <- Sys.getenv("ZHAO_VARIANCE_METHOD", "pooled_risk")
distributions <- c("exponential", "weibull", "lognormal", "loglogistic")
heterogeneity <- list(low = c(0.5, 1.5), medium = c(0.1, 1.9), high = c(0.01, 1.99))
tests <- c("logrank", "peto_prentice", "gehan_breslow")

run_cell <- function(distribution, z_range, test, n_sim, multiplier = 1,
                     alpha = 0.025, seed) {
  ans <- zhao_rank_simulation(
    n_sim = n_sim, n_per_group = n_per_group, scenario = "paper",
    distribution = distribution, heterogeneity = z_range, followup = followup,
    alpha = alpha, test = test, variance_method = variance_method,
    alternative = "less", treatment_time_multiplier = multiplier, seed = seed
  )
  data.frame(
    distribution = distribution, z_low = z_range[1], z_high = z_range[2],
    test = test, n_sim = n_sim, treatment_time_multiplier = multiplier,
    empirical_rejection = ans$empirical_type1_error,
    monte_carlo_se = ans$monte_carlo_se, mean_z = ans$mean_z, sd_z = ans$sd_z,
    mean_completed_events = ans$mean_completed_events
  )
}

# Zhao Table 1: distribution of standardized null statistics at medium
# heterogeneity. The paper reports 100,000 replicates and the empirical SD.
table1_jobs <- expand.grid(distribution = distributions, test = tests,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
table1 <- lapply(seq_len(nrow(table1_jobs)), function(i) {
  j <- table1_jobs[i, ]
  run_cell(j$distribution, heterogeneity$medium, j$test, n_null, seed = 1000 + i)
})

# Zhao Table 2: one-sided 0.025 null rejection rates, 100,000 replicates.
table2_jobs <- expand.grid(distribution = distributions, heterogeneity = names(heterogeneity),
  test = tests, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
table2 <- lapply(seq_len(nrow(table2_jobs)), function(i) {
  j <- table2_jobs[i, ]
  run_cell(j$distribution, heterogeneity[[j$heterogeneity]], j$test, n_null, seed = 2000 + i)
})

# Zhao Table 3: 10,000 one-sided power replicates. Each quoted alternative
# changes a distribution location/scale by 0.25 or 0.5 on the log-time scale.
table3_jobs <- expand.grid(distribution = distributions, heterogeneity = names(heterogeneity),
  log_time_shift = c(0.25, 0.5), test = tests,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
table3 <- lapply(seq_len(nrow(table3_jobs)), function(i) {
  j <- table3_jobs[i, ]
  run_cell(j$distribution, heterogeneity[[j$heterogeneity]], j$test, n_power,
    multiplier = exp(j$log_time_shift), seed = 3000 + i)
})

result <- list(
  source = "Zhao et al. (2020), Tables 1--3",
  settings = list(n_null = n_null, n_power = n_power, n_per_group = n_per_group,
    followup = followup, variance_method = variance_method,
    note = "The paper used survsim::rec.ev.surv(); this package uses its own transparent gap-time generator, so exact equality is not expected."),
  table1 = do.call(rbind, table1),
  table2 = do.call(rbind, table2),
  table3 = do.call(rbind, table3)
)

saveRDS(result, "data-raw/zhao-2020-table-simulations.rds")
print(result$table1)
print(result$table2)
print(result$table3)
