# Quick diagnostic for the Sivadasan--Sankaran two-cause Weibull frailty model.
# Run from the package root:
# Rscript data-raw/ss-rcif-quick-simulation.R [reps] [n] [censor_max]
# Defaults are deliberately small; e.g. use 1000 200 5 for a fuller comparison.

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 100L
n <- if (length(args) > 1L) as.integer(args[2]) else 100L
censor_max <- if (length(args) > 2L) as.numeric(args[3]) else 3
if (is.na(reps) || reps < 1L || is.na(n) || n < 2L ||
    is.na(censor_max) || !censor_max %in% c(3, 5)) {
  stop("reps must be >= 1, n must be >= 2, and censor_max must be 3 or 5")
}

source("R/utils.R")
source("R/sivadasan_sankaran.R")
set.seed(20231007)

lambda <- c(A = 0.1, B = 0.12)
total <- sum(lambda)
times <- c(1, 2)
truth <- outer(times, lambda / total,
               function(t, p) p * (1 - 1 / (1 + total * t^2)))
colnames(truth) <- names(lambda)

simulate_subjects <- function(n, censoring) {
  out <- vector("list", n)
  for (i in seq_len(n)) {
    z <- stats::rexp(1)
    limit <- if (censoring == "fixed") censor_max else
      stats::runif(1, 0, censor_max)
    elapsed <- 0
    rows <- list()
    repeat {
      gap <- sqrt(stats::rexp(1) / (z * total))
      if (elapsed + gap >= limit) {
        rows[[length(rows) + 1L]] <- data.frame(
          id = i, gap = limit - elapsed, status = 0L, cause = NA_character_
        )
        break
      }
      rows[[length(rows) + 1L]] <- data.frame(
        id = i, gap = gap, status = 1L,
        cause = sample(names(lambda), 1L, prob = lambda)
      )
      elapsed <- elapsed + gap
    }
    out[[i]] <- do.call(rbind, rows)
  }
  do.call(rbind, out)
}

settings <- list(
  literal_paper = list(cause_weighting = "literal_eq15", cif_method = "paper"),
  shared_paper = list(cause_weighting = "shared", cif_method = "paper"),
  shared_decrement = list(cause_weighting = "shared",
                          cif_method = "exponential_decrement")
)

estimate <- function(dat, setting) {
  fit <- do.call(ss_rcif, c(list(data = dat), setting))
  vapply(names(lambda), function(cause) {
    vapply(times, function(t) {
      j <- findInterval(t, fit$curve$time)
      if (j == 0L) 0 else fit$cause_curves[[cause]]$F[j]
    }, numeric(1))
  }, numeric(length(times)))
}

results <- list()
for (censoring in c("fixed", "uniform")) {
  estimates <- array(NA_real_, c(reps, length(times), length(lambda),
                                 length(settings)),
                     dimnames = list(NULL, times, names(lambda), names(settings)))
  for (b in seq_len(reps)) {
    dat <- simulate_subjects(n, censoring)
    if (length(unique(dat$cause[dat$status == 1L])) < 2L) next
    for (method in names(settings)) {
      estimates[b, , , method] <- estimate(dat, settings[[method]])
    }
  }
  for (method in names(settings)) {
    for (cause in names(lambda)) {
      for (t in seq_along(times)) {
        x <- estimates[, t, cause, method]
        x <- x[is.finite(x)]
        target <- truth[t, cause]
        results[[length(results) + 1L]] <- data.frame(
          censoring = paste0(censoring, censor_max), n = n, reps = length(x),
          method = method, cause = cause, time = times[t], truth = target,
          mean = mean(x), bias = mean(x - target),
          mse = mean((x - target)^2)
        )
      }
    }
  }
}
print(do.call(rbind, results), row.names = FALSE, digits = 4)
