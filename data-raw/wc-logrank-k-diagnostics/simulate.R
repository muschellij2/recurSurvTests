# Diagnose the K-group null calibration as subject count and allocation change.
# Run as a four-task array from the package root.
if (Sys.getenv("RUN_WC_K_DIAGNOSTIC") != "true") {
  stop("Set RUN_WC_K_DIAGNOSTIC=true to run this diagnostic simulation.")
}
if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools first.")
devtools::load_all(".", quiet = TRUE)

n_sim <- as.integer(Sys.getenv("WC_K_DIAGNOSTIC_N_SIM", "5000"))
task_id <- as.integer(Sys.getenv("WC_K_DIAGNOSTIC_TASK_ID", "1"))
seed <- as.integer(Sys.getenv("WC_K_DIAGNOSTIC_SEED", "20260928"))
if (anyNA(c(n_sim, task_id, seed)) || n_sim < 1L || !task_id %in% 1:4) {
  stop("Invalid diagnostic settings; task ID must be 1, 2, 3, or 4.")
}

# These allocations isolate the original small-sample unbalanced design,
# allocation balance at the same N, and increasing independent-cluster counts.
designs <- list(
  original_unbalanced = c(A = 30L, B = 45L, C = 60L),
  balanced_same_N = c(A = 45L, B = 45L, C = 45L),
  unbalanced_2xN = c(A = 60L, B = 90L, C = 120L),
  unbalanced_4xN = c(A = 120L, B = 180L, C = 240L)
)
arm_n <- designs[[task_id]]

one_subject <- function(id, arm) {
  total <- 0
  state <- stats::rbinom(1, 1, 0.5)
  frailty <- stats::rlnorm(1, 0, 0.55)
  end <- floor(stats::runif(1, 360, 480) / 0.5) * 0.5
  gap <- numeric()
  repeat {
    scale <- if (state == 1) 0.55 else 1.8
    next_gap <- ceiling(stats::rlnorm(1, log(2.4 * frailty * scale), 0.4) / 0.5) * 0.5
    if (total + next_gap >= end) {
      gap <- c(gap, end - total)
      break
    }
    gap <- c(gap, next_gap)
    total <- total + next_gap
    if (stats::runif(1) > 0.85) state <- 1 - state
  }
  data.frame(id = id, arm = arm, episode = seq_along(gap), gap = gap,
             status = c(rep(1L, length(gap) - 1L), 0L))
}
simulate_null <- function() {
  arms <- rep(names(arm_n), arm_n)
  do.call(rbind, lapply(seq_along(arms), function(i) one_subject(i, arms[i])))
}

raw <- vector("list", n_sim)
for (replicate_id in seq_len(n_sim)) {
  set.seed(seed + task_id * 1000000L + replicate_id)
  fit <- wc_logrank_k(simulate_null(), group = "arm", episode = "episode")
  q <- fit$df
  n <- sum(arm_n)
  # This F reference is included only as a finite-sample diagnostic. The test
  # under evaluation remains the documented chi-square Wald reference.
  p_f <- if (q > 0L && n > q) {
    stats::pf(((n - q) / (q * n)) * fit$chisq, q, n - q, lower.tail = FALSE)
  } else {
    NA_real_
  }
  raw[[replicate_id]] <- data.frame(
    design = names(designs)[task_id], n_total = n, n_A = arm_n[1],
    n_B = arm_n[2], n_C = arm_n[3], replicate = replicate_id,
    chisq = fit$chisq, df = q, p_chisq = fit$p.value, p_f = p_f
  )
}
raw <- do.call(rbind, raw)
outdir <- file.path("data-raw", "wc-logrank-k-diagnostics", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(settings = list(n_sim = n_sim, task_id = task_id, seed = seed,
  alpha = 0.05, arm_sizes = arm_n), raw = raw),
  file.path(outdir, sprintf("task-%02d-of-04.rds", task_id)))
message("Saved ", nrow(raw), " null replicates for ", names(designs)[task_id])
