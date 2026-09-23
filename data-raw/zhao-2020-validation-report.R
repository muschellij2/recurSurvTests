# Run from the package root; does not run or modify simulations.
# Rscript data-raw/zhao-2020-validation-report.R [input_directory] [output_directory]
# Requires only base R for CSV/Markdown/PNG; rmarkdown + Pandoc also produce HTML.

zhao_validation <- function(input = "data-raw/zhao-2020-table-simulations",
                            output = "docs/zhao-validation") {
  reference <- "data-raw/zhao-2020-published-tables.csv"
  paper <- read.csv(reference, stringsAsFactors = FALSE)
  key <- function(d) paste(d$table, d$distribution, d$heterogeneity, d$test,
                           sprintf("%.8f", d$log_time_shift), sep = "/")
  stopifnot(nrow(paper) == 120L, !anyDuplicated(key(paper)))
  # Preserve the exact job ordering in the simulation array driver.
  dist <- c("exponential", "weibull", "lognormal", "loglogistic")
  tests <- c("logrank", "peto_prentice", "gehan_breslow")
  t1 <- expand.grid(distribution = dist, test = tests, stringsAsFactors = FALSE)
  t1$table <- "table1"; t1$heterogeneity <- "medium"; t1$log_time_shift <- 0
  t2 <- expand.grid(distribution = dist, heterogeneity = c("low", "medium", "high"),
                    test = tests, stringsAsFactors = FALSE)
  t2$table <- "table2"; t2$log_time_shift <- 0
  t3 <- expand.grid(distribution = dist, heterogeneity = c("low", "medium", "high"),
                    log_time_shift = c(.25, .5), test = tests, stringsAsFactors = FALSE)
  t3$table <- "table3"
  cols <- c("table", "distribution", "heterogeneity", "test", "log_time_shift")
  jobs <- rbind(t1[cols], t2[cols], t3[cols])
  paper$task_id <- match(key(paper), key(jobs))
  stopifnot(!anyNA(paper$task_id))
  files <- sort(list.files(input, pattern = "^zhao-2020-table-simulations.*[.]rds$", full.names = TRUE))
  if (!length(files)) stop("No simulation RDS files found in ", input)
  objects <- lapply(files, readRDS)
  common <- c("n_null", "n_power", "n_per_group", "followup", "variance_method", "n_tasks")
  settings <- objects[[1]]$settings
  if (any(vapply(common, function(nm) is.null(settings[[nm]]), logical(1)))) stop("Missing settings")
  for (x in objects) for (nm in common) {
    if (!isTRUE(all.equal(x$settings[[nm]], settings[[nm]]))) stop("Mixed settings: ", nm)
  }
  if (settings$n_per_group != 100 || settings$followup != 180 || settings$n_tasks != 120) {
    stop("Expected the published 100-per-group, 180-day design and 120-task array")
  }
  ids <- vapply(objects, function(x) x$settings$task_id, numeric(1))
  if (anyNA(ids) || any(!ids %in% 1:120) || anyDuplicated(ids)) stop("Invalid or duplicate task IDs")
  rows <- lapply(seq_along(objects), function(i) {
    x <- objects[[i]]
    for (nm in c("table1", "table2", "table3")) {
      if (!is.data.frame(x[[nm]]) || any(x[[nm]]$table != nm)) stop("Malformed table in ", files[i])
    }
    d <- do.call(rbind, x[c("table1", "table2", "table3")])
    if (nrow(d) != 1L) stop("Expected one cell per array file: ", files[i])
    d$log_time_shift <- round(log(d$treatment_time_multiplier), 8)
    if (!identical(key(d), key(jobs[ids[i], , drop = FALSE]))) stop("Task/cell mismatch: ", files[i])
    expected_range <- switch(d$heterogeneity, low = c(.5, 1.5), medium = c(.1, 1.9), high = c(.01, 1.99))
    if (!isTRUE(all.equal(c(d$z_low, d$z_high), expected_range))) stop("Heterogeneity mismatch")
    n <- if (d$table == "table3") settings$n_power else settings$n_null
    numeric_cols <- c("n_sim", "empirical_rejection", "monte_carlo_se", "mean_z", "sd_z", "mean_completed_events")
    if (!all(numeric_cols %in% names(d)) || any(!is.finite(as.matrix(d[numeric_cols])))) stop("Nonfinite or absent summaries")
    if (d$n_sim != n || n < 2 || n != as.integer(n) || d$empirical_rejection < 0 ||
        d$empirical_rejection > 1 || d$sd_z <= 0 || d$mean_completed_events < 0) stop("Invalid summary")
    r <- d$empirical_rejection * n
    if (abs(r - round(r)) > 1e-6) stop("Rejection proportion does not imply an integer count")
    expected_se <- sqrt(d$empirical_rejection * (1 - d$empirical_rejection) / n)
    if (abs(d$monte_carlo_se - expected_se) > 1e-10) stop("Monte Carlo SE mismatch")
    d$generator_version <- if (is.null(x$settings$generator_version)) "legacy_unversioned" else x$settings$generator_version
    d$task_id <- ids[i]; d$input_file <- basename(files[i]); d
  })
  observed <- do.call(rbind, rows)
  if (anyDuplicated(key(observed))) stop("Duplicate simulation cells")
  d <- merge(paper, observed, by = c(cols, "task_id"), all.x = TRUE, sort = FALSE)
  stopifnot(nrow(d) == 120L, sum(!is.na(d$n_sim)) == length(files))
  d <- d[order(d$task_id), ]; rownames(d) <- NULL
  d$available <- !is.na(d$n_sim)
  d$variance_method <- settings$variance_method
  d$design_eligible <- d$distribution != "loglogistic" |
    (!is.na(d$generator_version) & d$generator_version == "zhao_table1_loglogistic_v2")
  d$estimate <- ifelse(d$table == "table1", d$sd_z, d$empirical_rejection)
  d$local_mc_se <- ifelse(d$table == "table1", d$sd_z / sqrt(2 * (d$n_sim - 1)), d$monte_carlo_se)
  d$paper_mc_se <- d$paper_value / sqrt(2 * (d$paper_n_sim - 1))
  rates <- d$table != "table1"
  d$paper_mc_se[rates] <- sqrt(d$paper_value[rates] * (1 - d$paper_value[rates]) / d$paper_n_sim[rates])
  d$difference <- d$estimate - d$paper_value
  d$combined_mc_se <- sqrt(d$local_mc_se^2 + d$paper_mc_se^2)
  d$comparison_margin <- qnorm(.975) * d$combined_mc_se + .0005
  d$comparison <- ifelse(!d$available, "missing", ifelse(!d$design_eligible, "design mismatch",
    ifelse(abs(d$difference) > d$comparison_margin, "discrepant", "compatible")))
  # Wilson intervals for local rejection probabilities (also reported for Table 1).
  z <- qnorm(.975); p <- d$empirical_rejection; n <- d$n_sim
  center <- (p + z^2 / (2*n)) / (1 + z^2/n)
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
  d$rejection_ci_low <- center - half; d$rejection_ci_high <- center + half
  d$nominal_p <- NA_real_
  ix <- which(d$table == "table2" & d$available & d$design_eligible)
  d$nominal_p[ix] <- vapply(ix, function(i) binom.test(round(p[i]*n[i]), n[i], p=.025)$p.value, numeric(1))
  d$nominal_p_holm <- NA_real_
  # Use all 36 planned Table 2 comparisons in multiplicity correction.
  d$nominal_p_holm[ix] <- p.adjust(d$nominal_p[ix], "holm", n = 36L)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(output, "comparison.csv"), row.names = FALSE, na = "")
  manifest <- data.frame(task_id = ids, file = basename(files), md5 = unname(tools::md5sum(files)))
  write.csv(manifest, file.path(output, "input-manifest.csv"), row.names = FALSE)
  provenance <- c(reference, "data-raw/zhao-2020-published-tables.md", "data-raw/zhao-2020-validation-report.R",
                  "data-raw/zhao-2020-table-simulations.R", "R/zhao_rank_tests.R")
  write.csv(data.frame(file = provenance, md5 = unname(tools::md5sum(provenance))),
            file.path(output, "source-manifest.csv"), row.names = FALSE)
  capture.output(str(settings), sessionInfo(), file = file.path(output, "session-info.txt"))
  eligible <- d$available & d$design_eligible
  counts <- do.call(rbind, lapply(c("table1", "table2", "table3"), function(t) {
    a <- d[d$table == t, ]
    data.frame(Table = t, Planned = nrow(a), Available = sum(a$available),
               Comparable = sum(a$available & a$design_eligible),
               Compatible = sum(a$comparison == "compatible"), Discrepant = sum(a$comparison == "discrepant"))
  }))
  # Export a standalone figure; missing and design-mismatched cells are excluded.
  png(file.path(output, "comparison.png"), width = 1800, height = 650, res = 150)
  par(mfrow = c(1,3), mar = c(4,4,3,1))
  colors <- c(logrank = "#0072B2", peto_prentice = "#D55E00", gehan_breslow = "#009E73")
  for (t in c("table1", "table2", "table3")) {
    a <- d[d$table == t & eligible, ]
    lim <- if (nrow(a)) range(c(a$paper_value, a$estimate - 1.96*a$local_mc_se,
                              a$estimate + 1.96*a$local_mc_se)) else c(0,1)
    plot(a$paper_value, a$estimate, xlim = lim, ylim = lim, pch = 19, col = colors[a$test],
         xlab = "Zhao published value", ylab = "Local estimate", main = switch(t,
           table1 = "Table 1: null statistic SD", table2 = "Table 2: type I error", table3 = "Table 3: power"))
    abline(0, 1, col = "grey50", lty = 2)
    segments(a$paper_value, a$estimate - 1.96*a$local_mc_se, a$paper_value,
             a$estimate + 1.96*a$local_mc_se, col = colors[a$test])
    legend("topleft", c("LR", "PP", "GB"), col = colors, pch = 19, bty = "n", cex = .7)
  }
  dev.off()
  fmt <- function(x, digits = 4) ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
  md_table <- function(x) {
    x[] <- lapply(x, as.character)
    c(paste0("| ", paste(names(x), collapse = " | "), " |"),
      paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
      apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")), "")
  }
  missing <- d$task_id[!d$available]
  worst <- d[which(eligible & d$table == "table3"), ]
  worst <- worst[order(abs(worst$difference), decreasing = TRUE), ]
  lines <- c("# Validation against Zhao et al. (2020) simulations", "",
    paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE), "."), "",
    "## Assessment", "",
    sprintf("**Published-table reproduction is not established.** %d of 120 planned cells are available; %d are comparable under the design audit below. Of these, %d are compatible and %d are discrepant using the exploratory Monte Carlo criterion defined below.",
      sum(d$available), sum(eligible), sum(d$comparison == "compatible"), sum(d$comparison == "discrepant")), "",
    sprintf("These results assess `zhao_rank_test(variance_method = \"%s\")` against the paper's **New method** columns. They do not establish that the implementation reproduces the printed Equation 6 variance, and they do not validate `wc_logrank()` (the separate Luo–Huang test).", settings$variance_method), "",
    md_table(counts),
    sprintf("Among available, design-comparable Table 2 cells, %d of %d differ from nominal 0.025 after Holm adjustment over all 36 planned cells (exact two-sided binomial tests, familywise 0.05). Agreement with a paper estimate and calibration to nominal alpha are separate questions.",
      sum(d$nominal_p_holm < .05, na.rm = TRUE), length(ix)), "",
    if (nrow(worst)) sprintf("Across available comparable power cells, local minus published power ranges from %+.2f to %+.2f percentage points. These differences require investigation of both the data-generating design and the implementation; this report does not identify a unique cause.",
      100 * min(worst$difference), 100 * max(worst$difference)) else "", "",
    if (nrow(worst)) sprintf("The largest available power difference is %s / %s heterogeneity / %s / log-time shift %.2f: local %.4f versus Zhao %.3f (difference %+.4f).",
      worst$distribution[1], worst$heterogeneity[1], worst$test[1], worst$log_time_shift[1],
      worst$estimate[1], worst$paper_value[1], worst$difference[1]) else "No comparable power cells are available.", "",
    "![Published versus local estimates](comparison.png)", "",
    "Vertical bars show approximate local 95% Monte Carlo intervals; the diagonal denotes equality. The classification uses uncertainty from both simulations plus paper rounding. Missing cells and the log-logistic design mismatch are excluded from this figure.", "",
    "## Source and benchmark interpretation", "",
    "Zhao Q, Zhang B, LaValley MP, Massaro JM, Lunetta KL, Chang M (2020). *Extended Rank Tests for Analyzing Recurrent Event Data*. Statistics in Biopharmaceutical Research 12(1):90–98. [DOI](https://doi.org/10.1080/19466315.2019.1601596). Benchmarks were transcribed from the locally supplied publisher PDF: Table 1, page 93; Table 2, page 94; Table 3, page 95. Only the New method LR, PP, and GB columns are used.", "",
    "**Table 2 scale ambiguity:** the caption says `(*100)`, but the cells visibly print 27, 26, 25, etc. Section 3.2 specifies one-sided alpha 0.025 and describes controlled type I error. This report interprets 27 as 0.027, dividing the printed values by 1000. This inference is not author-confirmed. All Table 2 benchmark comparisons are conditional on it; the raw printed values are preserved in the reference CSV.", "",
    "LR = log-rank, PP = Peto–Prentice, GB = Gehan–Breslow. These tests target marginal recurrent **gap-time survival**, not cumulative event burden over calendar time or recurrent cause-specific cumulative incidence.", "",
    "## Design audit", "",
    "The paper uses 100 subjects per group, study end 180, entry U(0,180), and subject gap multipliers Z drawn from U(0.5,1.5), U(0.1,1.9), or U(0.01,1.99). Tables 1–2 use 100,000 datasets per cell; Table 3 uses 10,000. Local metadata is checked against the subject count and follow-up; actual replicate counts appear below. Local seeds and the generator are different from the paper's `survsim::rec.ev.surv()` implementation. The paper also describes noninformative censoring without supplying its full generating parameters in Section 3; the local generator uses administrative censoring only. Exact design equivalence cannot be assumed.", "",
    "For exponential, Weibull, and log-normal baseline gaps, the current source matches the distributions printed in Table 1: rate exp(-4); Weibull shape 2 and scale exp(4); and log-normal meanlog 4, sdlog 0.5. Table 3 multiplies treatment gaps by exp(0.25) or exp(0.5). For Weibull shape 2, these correctly correspond to changes of 0.5 or 1 in -log(lambda).", "",
    "**Log-logistic mismatch:** Table 1 specifies S(t) = 1 / [1 + (lambda*t)^(1/gamma)], lambda = exp(-4), gamma = 0.5. An inverse draw is `exp(4) * (U^(-1) - 1)^0.5`, with mean exp(4)*pi/2, approximately 85.76. The legacy generator used `exp(-4) * (U^(-0.5) - 1)^0.5`, with mean exp(-4)*pi/4, approximately 0.01438. Both the scale and inverse transform were incorrect. The production generator is now corrected. New outputs tagged `zhao_table1_loglogistic_v2` are eligible for comparison; unversioned log-logistic outputs remain excluded as design mismatch. Other distributions are unaffected by this correction.", "",
    "All supplied outputs in this report must share the same variance method. `pooled_risk` uses a pooled-risk residual; `zhao_eq6` is the separate group-risk residual convention printed in Equation 6. Calibration of one does not settle the validity or intended interpretation of the other.", "",
    "### Why every fourth array task can exhaust memory", "",
    "Task IDs 4, 8, 12, and so on select log-logistic gaps. Under the legacy incorrect draw, a renewal approximation gives roughly 1.37, 2.05, and 3.35 million completed events per null dataset at low, medium, and high heterogeneity: 200*90*E(1/Z)/E(base gap), with E(1/Z)=log(b/a)/(b-a). Almost all continuous event times are distinct. Each dense 200-by-event-count double matrix therefore occupies about 2.05, 3.05, or 4.98 GiB. The rank test retains roughly six such matrices, with additional temporary arrays. Across four workers, those six matrices alone account for approximately 49, 73, or 120 GiB. These are approximate typical null-design allocations, not guaranteed peak bounds. A 25 GiB allocation is consequently insufficient for the legacy generator. The nested risk-set calculations also become extremely slow.", "",
    "A separate corrected-generator sizing pilot (`data-raw/zhao-2020-memory-pilot.R`) draws five null datasets at each heterogeneity level in an isolated R environment. The September 23, 2026 local macOS run produced 222–676 completed events per dataset and a maximum resident set size of 196,788,224 bytes (about 188 MiB) for one R process. This is a small pilot, not a cluster or full-run peak measurement. With the corrected production generator, retaining 16–25 GiB total per four-worker task is a conservative starting allocation; measure cluster MaxRSS on a short pilot before scaling up. Increasing memory for the incorrect generator would not make its results a validation of Zhao's simulations.", "",
    "## Monte Carlo comparison method", "",
    "For rejection rates, the local Monte Carlo SE is sqrt(p*(1-p)/n); local 95% intervals use Wilson's formula. The paper SE uses its reported probability and replicate count. A cell is flagged discrepant when |local - paper| exceeds 1.96*sqrt(SE_local^2 + SE_paper^2) + 0.0005. The final term allows half a unit of paper rounding. Independent simulation runs are assumed. This is an exploratory cellwise screen, not an equivalence test or a multiplicity-adjusted global pass/fail rule. Compatible means no discrepancy detected by this screen, not proven equality.", "",
    "Table 1 compares empirical SDs, using the normal-theory approximation SE(SD) = SD/sqrt(2*(n-1)) for both runs. This is approximate because replicate-level fourth moments are unavailable. It does not test tail behavior or the full null distribution. Raw replicate z scores and p values were not retained in the array RDS files. Although nonfinite summary moments are rejected by this aggregator, the original driver's rejection count used `na.rm = TRUE`; failed-p-value counts cannot be independently audited from these summaries.", "",
    "## Completeness and provenance", "",
    paste0("Missing task IDs: ", if (length(missing)) paste(missing, collapse = ", ") else "none", "."), "",
    "The aggregator rejects duplicate task IDs, duplicate cells, inconsistent settings, unexpected task-to-cell mappings, invalid counts, nonfinite moments, and inconsistent Monte Carlo SEs. All 120 planned cells are retained in the tables and CSV; missing results are not treated as zero rejections. `input-manifest.csv` records input MD5 checksums; `source-manifest.csv` records the current reference, driver, implementation, and reporting sources. These are current-source checksums: the original RDS files do not record the execution commit or source hashes, so the historical code cannot be proven from these files. The design audit distinguishes the legacy generator from the corrected version recorded in new output metadata. Runtime settings and R session information are in `session-info.txt`.", "",
    "## Cell-by-cell comparison", "",
    "Differences are local minus published, on the SD scale for Table 1 and probability scale for Tables 2–3 (multiply by 100 for percentage points). A dash means the task is unavailable. Heterogeneity is low U(0.5,1.5), medium U(0.1,1.9), or high U(0.01,1.99).")
  for (t in c("table1", "table2", "table3")) {
    a <- d[d$table == t, ]
    tab <- data.frame(Task = a$task_id, Distribution = a$distribution,
      Heterogeneity = a$heterogeneity, Test = c(logrank="LR", peto_prentice="PP", gehan_breslow="GB")[a$test],
      Shift = fmt(a$log_time_shift, 2), N = ifelse(a$available, a$n_sim, "—"),
      Zhao = fmt(a$paper_value, 3), Local = fmt(a$estimate), Difference = fmt(a$difference),
      MCSE = fmt(a$local_mc_se), Result = a$comparison, check.names = FALSE)
    if (t == "table1") tab$Mean_z <- fmt(a$mean_z)
    if (t == "table2") {
      tab$`Local 95% CI` <- ifelse(a$available, paste0(fmt(a$rejection_ci_low), "–", fmt(a$rejection_ci_high)), "—")
      tab$`Nominal Holm p` <- fmt(a$nominal_p_holm)
    }
    lines <- c(lines, "", paste0("### ", switch(t, table1="Table 1: null statistic SD", table2="Table 2: type I error", table3="Table 3: power")), "", md_table(tab))
  }
  lines <- c(lines, "## Next steps implied by this validation", "",
    "1. Recover the remaining non-log-logistic task outputs and regenerate this report.",
    "2. Rerun the corrected log-logistic cells with `sbatch --array=4-120:4%20 data-raw/zhao-2020-table-simulations.sbatch`. The script preserves the original 120-task partition and new outputs identify the corrected generator. Archive any legacy log-logistic outputs before replacement. Other distributions need not be rerun for this generator fix.",
    "3. Investigate reproducible discrepancies in the comparable cells, checking the generator/censoring design and variance convention against the paper and, if obtainable, author code. Resolve the Table 2 scaling ambiguity with the authors.",
    "4. Run targeted paired-seed comparisons of the pooled-risk and Equation 6 conventions before making an exact-reproduction claim. Matching nominal size alone does not validate the reported power results.", "",
    "## Regenerate", "", "From the package root:", "", "```sh",
    "Rscript data-raw/zhao-2020-validation-report.R", "```", "",
    "Optional arguments select input and output directories. This command reads existing results and does not launch simulations. It writes Markdown, CSVs, a PNG figure, provenance files, and (when rmarkdown/Pandoc are available) a standalone HTML report.")
  md <- file.path(output, "report.md")
  writeLines(lines, md, useBytes = TRUE)
  if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
    rmarkdown::pandoc_convert(normalizePath(md), from = "markdown", to = "html5",
      output = file.path(normalizePath(output), "report.html"),
      options = c("--standalone", "--self-contained", "--metadata", "pagetitle=Zhao simulation validation"))
  }
  print(counts)
  message("Report written to ", md)
  invisible(d)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 2L) stop("Usage: Rscript data-raw/zhao-2020-validation-report.R [input] [output]")
  do.call(zhao_validation, as.list(args))
}
