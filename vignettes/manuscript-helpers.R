# Reporting helpers; no simulation arrays are launched by manuscript rendering.
manuscript_zhao <- function(root, input = Sys.getenv("ZHAO_RESULTS_DIR", "data-raw/zhao-2020-table-simulations")) {
  if (!grepl("^(/|[A-Za-z]:)", input)) input <- file.path(root, input)
  reference <- file.path(root, "data-raw/zhao-2020-published-tables.csv")
  files <- sort(list.files(input, pattern = "^zhao-2020-table-simulations.*[.]rds$", full.names = TRUE))
  if (length(files)) {
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "data-raw/zhao-2020-validation-report.R"), envir = env)
    d <- env$zhao_validation(input = input, reference = reference, write_outputs = FALSE)
  } else {
    # Preserve the complete published design when no simulation files are supplied.
    d <- read.csv(reference, stringsAsFactors = FALSE)
    for (nm in c("n_sim", "estimate", "local_mc_se", "difference", "comparison_margin",
                 "rejection_ci_low", "rejection_ci_high", "mean_z", "nominal_p_holm")) d[[nm]] <- NA_real_
    d$available <- d$design_eligible <- FALSE
    d$comparison <- "missing"
    d$variance_method <- "not available"
  }
  attr(d, "input_manifest") <- data.frame(file = basename(files), md5 = unname(tools::md5sum(files)))
  d
}

manuscript_zhao_status <- function(d) {
  do.call(rbind, lapply(c("table1", "table2", "table3"), function(nm) {
    a <- d[d$table == nm, ]
    data.frame(Table = sub("table", "", nm), Planned = nrow(a), Available = sum(a$available),
      Comparable = sum(a$available & a$design_eligible),
      Missing = sum(!a$available), Excluded = sum(a$available & !a$design_eligible),
      Replicates_per_available_cell = if (any(a$available))
        paste(unique(a$n_sim[a$available]), collapse = ", ") else "Pending",
      Discrepant = sum(a$comparison == "discrepant"))
  }))
}

# Wide tables retain the row layout and LR/PP/GB column order of Zhao Tables 1-3.
manuscript_zhao_table <- function(d, table) {
  a <- d[d$table == table, ]
  row_cols <- if (table == "table1") "distribution" else c("heterogeneity", "distribution")
  if (table == "table3") row_cols <- c(row_cols, "log_time_shift")
  base <- unique(a[row_cols])
  ord <- lapply(row_cols, function(nm) switch(nm,
    distribution = match(base[[nm]], c("exponential", "weibull", "lognormal", "loglogistic")),
    heterogeneity = match(base[[nm]], c("low", "medium", "high")), base[[nm]]))
  base <- base[do.call(order, ord), , drop = FALSE]
  key <- function(x) do.call(paste, c(x[row_cols], sep = "/"))
  for (test in c("logrank", "peto_prentice", "gehan_breslow")) {
    b <- a[a$test == test, ]; b <- b[match(key(base), key(b)), ]
    abbr <- c(logrank = "LR", peto_prentice = "PP", gehan_breslow = "GB")[[test]]
    base[[paste(abbr, "Zhao")]] <- sprintf("%.3f", b$paper_value)
    base[[paste(abbr, "local (MCSE)")]] <- ifelse(!b$available, "Pending",
      ifelse(!b$design_eligible, "Excluded", sprintf("%.4f (%.4f)%s", b$estimate, b$local_mc_se,
        ifelse(b$comparison == "discrepant", " *", ""))))
  }
  names(base)[match(row_cols, names(base))] <- c(distribution="Distribution",
    heterogeneity="Heterogeneity", log_time_shift="Log-time shift")[row_cols]
  rownames(base) <- NULL
  base
}

manuscript_zhao_plot <- function(d) {
  old <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  on.exit(par(old))
  colors <- c(logrank = "#0072B2", peto_prentice = "#D55E00", gehan_breslow = "#009E73")
  for (nm in c("table1", "table2", "table3")) {
    a <- d[d$table == nm & d$available & d$design_eligible, ]
    if (!nrow(a)) {
      plot.new(); title(nm); text(.5, .5, "Results pending"); next
    }
    lower <- a$estimate - 1.96 * a$local_mc_se
    upper <- a$estimate + 1.96 * a$local_mc_se
    if (nm != "table1") { lower <- a$rejection_ci_low; upper <- a$rejection_ci_high }
    lim <- range(c(a$paper_value, lower, upper))
    plot(a$paper_value, a$estimate, xlim = lim, ylim = lim,
      pch = 19, col = colors[a$test], xlab = "Zhao published value", ylab = "Local estimate",
      main = c(table1="Null statistic SD", table2="Type I error", table3="Power")[[nm]])
    abline(0, 1, lty = 2, col = "grey50")
    segments(a$paper_value, lower, a$paper_value, upper, col = colors[a$test])
    legend("topleft", c("LR", "PP", "GB"), col = colors, pch = 19, bty = "n", cex = .8)
  }
}

manuscript_reference_data <- function(seed, n, event_range, rate) {
  set.seed(seed)
  counts <- sample(event_range, n, replace = TRUE)
  do.call(rbind, lapply(seq_len(n), function(i) {
    k <- counts[i]
    data.frame(id = i, episode = seq_len(k + 1L), gap = rexp(k + 1L, rate),
               status = c(rep(1L, k), 0L))
  }))
}

manuscript_reference_checks <- function() {
  scenarios <- list(
    `100 subjects, 1-5 events` = manuscript_reference_data(20260918, 100, 1:5, 1/12),
    `80 subjects, 0-6 events` = manuscript_reference_data(20260919, 80, 0:6, 1/9)
  )
  rows <- list()
  for (scenario in names(scenarios)) for (estimator in c("WC", "PSH")) {
    dat <- scenarios[[scenario]]
    ours <- if (estimator == "WC") recurSurvTests::wc_surv(dat, episode="episode") else
      recurSurvTests::psh_surv(dat, episode="episode")
    for (pkg in c("survrec", "newTestSurvRec")) {
      available <- requireNamespace(pkg, quietly = TRUE)
      difference <- NA_real_
      if (available) {
        if (pkg == "survrec") {
          fun <- getExportedValue(pkg, if (estimator == "WC") "wc.fit" else "psh.fit")
          ref <- fun(survrec::Survr(dat$id, dat$gap, dat$status))
        } else {
          fun <- getExportedValue(pkg, if (estimator == "WC") "WC.fit" else "PSH.fit")
          ref <- fun(as.matrix(dat[c("id", "gap", "status")]))
        }
        # Fail on a mismatched grid rather than silently recycling curve values.
        stopifnot(length(ours$time) == length(ref$time),
                  isTRUE(all.equal(ours$time, as.numeric(ref$time), tolerance=1e-7)))
        difference <- max(abs(ours$surv - ref$survfunc))
      }
      rows[[length(rows)+1L]] <- data.frame(Scenario=scenario, Estimator=estimator,
        Reference=pkg, Version=if(available) as.character(utils::packageVersion(pkg)) else "not installed",
        Max_absolute_difference=difference)
    }
  }
  do.call(rbind, rows)
}
