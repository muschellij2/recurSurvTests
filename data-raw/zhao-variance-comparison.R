# Compare the published-table summaries for pooled-risk and printed Eq. 6
# variance conventions. Generate each report first, from the matching task RDS
# directory, then run this script from the package root.

pooled_file <- file.path("docs", "zhao-validation", "comparison.csv")
eq6_file <- file.path("docs", "zhao-validation-zhao-eq6", "comparison.csv")
if (!file.exists(pooled_file) || !file.exists(eq6_file)) {
  stop("Generate both pooled-risk and zhao_eq6 validation reports first.")
}
pooled <- read.csv(pooled_file, stringsAsFactors = FALSE)
eq6 <- read.csv(eq6_file, stringsAsFactors = FALSE)
keys <- c("table", "distribution", "heterogeneity", "test", "log_time_shift", "task_id")
if (anyDuplicated(pooled[keys]) || anyDuplicated(eq6[keys])) {
  stop("Duplicate table cells in a variance report")
}
names(pooled)[names(pooled) == "estimate"] <- "estimate_pooled"
names(eq6)[names(eq6) == "estimate"] <- "estimate_eq6"
names(pooled)[names(pooled) == "comparison"] <- "comparison_pooled"
names(eq6)[names(eq6) == "comparison"] <- "comparison_eq6"
names(pooled)[names(pooled) == "variance_method"] <- "variance_method_pooled"
names(eq6)[names(eq6) == "variance_method"] <- "variance_method_eq6"
common <- merge(
  pooled[c(keys, "paper_value", "available", "design_eligible", "estimate_pooled",
           "local_mc_se", "comparison_pooled", "variance_method_pooled")],
  eq6[c(keys, "estimate_eq6", "local_mc_se", "comparison_eq6", "variance_method_eq6")],
  by = keys, all = TRUE, sort = FALSE, suffixes = c("_pooled", "_eq6")
)
if (anyNA(common$estimate_pooled) || anyNA(common$estimate_eq6)) {
  stop("The two reports do not contain the same complete set of estimates")
}
if (any(common$variance_method_pooled != "pooled_risk") ||
    any(common$variance_method_eq6 != "zhao_eq6")) {
  stop("Unexpected variance-method metadata")
}
common$eq6_minus_pooled <- common$estimate_eq6 - common$estimate_pooled
common$abs_error_pooled <- abs(common$estimate_pooled - common$paper_value)
common$abs_error_eq6 <- abs(common$estimate_eq6 - common$paper_value)
common$eq6_closer_to_paper <- common$abs_error_eq6 < common$abs_error_pooled

out_dir <- file.path("docs", "zhao-validation-variance-comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(common, file.path(out_dir, "variance-comparison.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(common, common$table), function(x) {
  data.frame(
    table = x$table[1], cells = nrow(x),
    eq6_closer_cells = sum(x$eq6_closer_to_paper),
    pooled_mean_absolute_error = mean(x$abs_error_pooled),
    eq6_mean_absolute_error = mean(x$abs_error_eq6),
    mean_eq6_minus_pooled = mean(x$eq6_minus_pooled)
  )
}))
rownames(summary) <- NULL
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summary)
