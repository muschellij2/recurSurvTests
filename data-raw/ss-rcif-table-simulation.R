# Reproduce the six Sivadasan--Sankaran (2023) simulation designs.
# From the package root, for one of 24 table/sample-size cells:
# RUN_SS_RCIF_TABLE_SIMULATION=true SS_RCIF_TASK_ID=1 SS_RCIF_REPS=1000 \
#   Rscript data-raw/ss-rcif-table-simulation.R
# Task 1-4 = Table 1 n=50,100,200,500; 5-8 = Table 2; etc.

if (Sys.getenv("RUN_SS_RCIF_TABLE_SIMULATION") != "true") {
  stop("Set RUN_SS_RCIF_TABLE_SIMULATION=true to run")
}
task_id <- as.integer(Sys.getenv("SS_RCIF_TASK_ID", "1"))
reps <- as.integer(Sys.getenv("SS_RCIF_REPS", "1000"))
seed <- as.integer(Sys.getenv("SS_RCIF_SEED", "20231007"))
outdir <- Sys.getenv("SS_RCIF_RESULTS_DIR", "data-raw/ss-rcif-tables/results")
if (anyNA(c(task_id, reps, seed)) || !task_id %in% 1:24 || reps < 1L) {
  stop("Invalid SS_RCIF_TASK_ID, SS_RCIF_REPS, or SS_RCIF_SEED")
}

source("R/utils.R")
source("R/sivadasan_sankaran.R")

design <- data.frame(
  table = 1:6,
  censoring = c("uniform", "uniform", "uniform", "uniform", "fixed", "fixed"),
  censor_limit = c(3, 5, 3, 5, 3, 5),
  subject_weight = c("followup", "followup", "one", "one", "followup", "followup"),
  stringsAsFactors = FALSE
)
scenario <- design[(task_id - 1L) %/% 4L + 1L, , drop = FALSE]
n <- c(50L, 100L, 200L, 500L)[(task_id - 1L) %% 4L + 1L]
times <- seq(0.2, 1.6, by = 0.2)
lambda <- c(A = 0.1, B = 0.12)
lambda_total <- sum(lambda)
truth <- outer(times, lambda / lambda_total,
  function(t, p) p * (1 - 1 / (1 + lambda_total * t^2)))
colnames(truth) <- names(lambda)

methods <- list(
  literal_right = list(cause_weighting = "literal_eq15", cif_method = "paper",
                       survival_side = "right"),
  literal_left = list(cause_weighting = "literal_eq15", cif_method = "paper",
                      survival_side = "left"),
  shared_right = list(cause_weighting = "shared", cif_method = "paper",
                      survival_side = "right"),
  shared_decrement = list(cause_weighting = "shared",
                          cif_method = "exponential_decrement")
)

simulate_data <- function() {
  subjects <- vector("list", n)
  for (i in seq_len(n)) {
    z <- stats::rexp(1, rate = 1)
    censor_time <- if (scenario$censoring == "fixed") scenario$censor_limit else
      stats::runif(1, min = 0, max = scenario$censor_limit)
    elapsed <- 0
    rows <- list()
    repeat {
      gap <- sqrt(stats::rexp(1) / (z * lambda_total))
      if (elapsed + gap >= censor_time) {
        rows[[length(rows) + 1L]] <- data.frame(
          id = i, gap = censor_time - elapsed, status = 0L,
          cause = NA_character_
        )
        break
      }
      rows[[length(rows) + 1L]] <- data.frame(
        id = i, gap = gap, status = 1L,
        cause = sample(names(lambda), 1L, prob = lambda)
      )
      elapsed <- elapsed + gap
    }
    subjects[[i]] <- do.call(rbind, rows)
  }
  do.call(rbind, subjects)
}

extract_cif <- function(fit) {
  vapply(names(lambda), function(cause) {
    j <- findInterval(times, fit$curve$time)
    c(0, fit$cause_curves[[cause]]$F)[j + 1L]
  }, numeric(length(times)))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, sprintf("table-%02d-n%03d-reps%04d.rds",
  scenario$table, n, reps))
if (file.exists(outfile)) stop("Output already exists: ", outfile)

raw <- vector("list", reps * length(methods) * length(lambda) * length(times))
failures <- character()
cursor <- 0L
started <- proc.time()[["elapsed"]]
for (b in seq_len(reps)) {
  set.seed(seed + task_id * 100000L + b)
  dat <- simulate_data()
  events <- dat[dat$status == 1L, , drop = FALSE]
  singleton_count <- vapply(names(lambda), function(cause) {
    counts <- table(factor(events$id[events$cause == cause], levels = seq_len(n)))
    sum(counts == 1L)
  }, integer(1))
  if (length(unique(events$cause)) < 2L) {
    failures <- c(failures, sprintf("replicate %d: fewer than two observed causes", b))
    next
  }
  for (method in names(methods)) {
    fit <- try(do.call(ss_rcif, c(
      list(data = dat, subject_weight = scenario$subject_weight), methods[[method]]
    )), silent = TRUE)
    if (inherits(fit, "try-error")) {
      failures <- c(failures, sprintf("replicate %d, %s: %s", b, method,
                                      as.character(fit)[1]))
      next
    }
    estimate <- extract_cif(fit)
    for (cause in names(lambda)) for (j in seq_along(times)) {
      cursor <- cursor + 1L
      raw[[cursor]] <- data.frame(
        task_id = task_id, table = scenario$table, n = n, replicate = b,
        censoring = scenario$censoring, censor_limit = scenario$censor_limit,
        subject_weight = scenario$subject_weight, method = method,
        cause = cause, time = times[j], truth = truth[j, cause],
        estimate = estimate[j, cause],
        singleton_subjects = singleton_count[cause],
        completed_events = nrow(events),
        probability_conservation_max_abs = fit$probability_conservation_max_abs
      )
    }
  }
  if (b %% max(1L, ceiling(reps / 20L)) == 0L || b == reps) {
    message("Task ", task_id, ": ", b, "/", reps, " replicates")
  }
}
raw <- if (cursor) do.call(rbind, raw[seq_len(cursor)]) else data.frame()
metadata <- list(
  task_id = task_id, table = scenario$table, n = n, requested_reps = reps,
  seed = seed, censoring = scenario$censoring,
  censor_limit = scenario$censor_limit, subject_weight = scenario$subject_weight,
  methods = methods, times = times, lambda = lambda,
  source_md5 = tools::md5sum(c("R/utils.R", "R/sivadasan_sankaran.R",
                               "data-raw/ss-rcif-table-simulation.R")),
  elapsed_seconds = proc.time()[["elapsed"]] - started,
  run_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string
)
saveRDS(list(metadata = metadata, raw = raw, failures = failures), outfile)
message("Saved ", cursor, " estimates; ", length(failures),
        " failures; ", round(metadata$elapsed_seconds, 1), " seconds: ", outfile)
