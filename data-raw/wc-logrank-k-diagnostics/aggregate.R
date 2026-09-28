# Aggregate the K-group sample-size null diagnostic.
result_dir <- file.path("data-raw", "wc-logrank-k-diagnostics", "results")
files <- sort(list.files(result_dir, pattern = "^task-[0-9]+-of-04[.]rds$",
                         full.names = TRUE))
if (length(files) != 4L) stop("Expected all four diagnostic task files")
objects <- lapply(files, readRDS)
ids <- vapply(objects, function(x) x$settings$task_id, integer(1))
if (!identical(sort(ids), 1:4)) stop("Invalid or duplicate task IDs")
raw <- do.call(rbind, lapply(objects, `[[`, "raw"))
alpha <- objects[[1]]$settings$alpha
summary <- do.call(rbind, lapply(split(raw, raw$design), function(x) {
  n <- nrow(x)
  k_chi <- sum(x$p_chisq < alpha)
  k_f <- sum(x$p_f < alpha)
  ci_chi <- stats::binom.test(k_chi, n)$conf.int
  ci_f <- stats::binom.test(k_f, n)$conf.int
  data.frame(
    design = x$design[1], n_total = x$n_total[1], n_sim = n,
    size_chisq = k_chi / n,
    mc_se_chisq = sqrt((k_chi / n) * (1 - k_chi / n) / n),
    ci_low_chisq = ci_chi[1], ci_high_chisq = ci_chi[2],
    p_nominal_chisq = stats::binom.test(k_chi, n, p = alpha)$p.value,
    size_f = k_f / n,
    mc_se_f = sqrt((k_f / n) * (1 - k_f / n) / n),
    ci_low_f = ci_f[1], ci_high_f = ci_f[2],
    p_nominal_f = stats::binom.test(k_f, n, p = alpha)$p.value
  )
}))
rownames(summary) <- NULL
saveRDS(list(summary = summary, raw = raw),
        file.path("data-raw", "wc-logrank-k-diagnostics", "results.rds"))
utils::write.csv(summary,
  file.path("data-raw", "wc-logrank-k-diagnostics", "results.csv"),
  row.names = FALSE)
print(summary)
