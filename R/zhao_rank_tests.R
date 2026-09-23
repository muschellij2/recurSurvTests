#' Zhao et al. extended recurrent gap-time rank test
#'
#' Tests equality of two Wang--Chang marginal recurrent gap-time survival
#' distributions using the extended rank-test construction of Zhao et al.
#' (2020).  Unlike [wc_logrank()], this function uses Zhao et al.'s adjusted
#' subject-level residuals (their Equation 6) to estimate the robust variance.
#'
#' @param data A data frame in long recurrent gap-time format.
#' @param group Character scalar naming a subject-level two-group variable.
#' @param id,gap,status,episode See [wc_logrank()].
#' @param test Rank weight: `"logrank"` (LR), `"gehan_breslow"` (GB), or
#'   `"peto_prentice"` (PP).
#' @param variance_method Robust variance residual. The default `"pooled_risk"`
#'   uses a pooled risk-set denominator and has residuals that sum exactly to
#'   the score. `"zhao_eq6"` follows the group-specific denominator printed in
#'   Zhao et al.'s Equation 6 as a sensitivity/reproduction implementation.
#' @param alternative Alternative for the normal-score p-value.
#' @param tau Optional upper truncation time.
#'
#' @return A list with the unstandardized score `W`, robust variance `variance`,
#' standard error, normal score `z`, p-value, subject residuals, and event-time
#' details. `group1` is the treatment group, determined by second appearance
#' among subjects; a positive score corresponds to more events in `group1`.
#'
#' @details
#' The implementation follows Equations 3, 6, and 7 of Zhao et al. (2020).
#' The paper presents its calculation without tied observed times. Here, tied
#' events are handled by aggregating their Wang--Chang weighted increments at a
#' common time; this is a natural extension, but is not separately derived in
#' that paper. Peto--Prentice weights use the pooled Wang--Chang survival
#' estimate `exp(-cumulative hazard)`, as defined in the paper.
#'
#' @references
#' Zhao Q, Zhang B, LaValley MP, Massaro JM, Lunetta KL, Chang M (2020).
#' Extended rank tests for analyzing recurrent event data. *Statistics in
#' Biopharmaceutical Research*, 12(1), 90-98.
#' \doi{10.1080/19466315.2019.1601596}.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:4, each = 3), episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4),
#'   arm = rep(c("A", "A", "B", "B"), each = 3)
#' )
#' zhao_rank_test(dat, group = "arm", episode = "episode")
#'
#' @export
zhao_rank_test <- function(data, group, id = "id", gap = "gap",
                           status = "status", episode = NULL,
                           test = c("logrank", "gehan_breslow", "peto_prentice"),
                           variance_method = c("pooled_risk", "zhao_eq6"),
                           alternative = c("two.sided", "greater", "less"),
                           tau = Inf) {
  test <- match.arg(test)
  variance_method <- match.arg(variance_method)
  alternative <- match.arg(alternative)
  d <- .validate_gap_data(data, id, gap, status, episode)
  if (!group %in% names(d)) stop("Missing group column: ", group)

  recs <- .subject_records(d, id, gap, status, group = group)
  lev <- unique(vapply(recs, function(r) as.character(r$group), character(1)))
  if (length(lev) != 2L) stop("zhao_rank_test currently requires exactly two groups")
  group0 <- lev[1]
  group1 <- lev[2]
  Z <- as.numeric(vapply(recs, function(r) as.character(r$group), character(1)) == group1)
  times <- .event_grid(recs, tau)
  if (!length(times)) stop("No completed recurrent gaps within tau")

  cc <- .wrs_components(recs, times)
  risk <- colSums(cc$R_i)
  dN <- colSums(cc$dN_i)
  risk1 <- colSums(cc$R_i * Z)
  risk_by_group <- rbind(colSums(cc$R_i * (1 - Z)), risk1)
  frac1 <- ifelse(risk > 0, risk1 / risk, 0)
  d_lambda <- ifelse(risk > 0, dN / risk, 0)
  cumhaz <- cumsum(d_lambda)
  survival <- exp(-cumhaz)
  weight <- switch(
    test,
    logrank = rep(1, length(times)),
    gehan_breslow = risk / length(recs),
    peto_prentice = survival
  )

  # Equation 3 score, written as the subject-weighted event increments.
  centered_group <- outer(Z, frac1, "-")
  score_contribution <- rowSums(sweep(cc$dN_i, 2, weight, "*") * centered_group)
  W <- sum(score_contribution)

  # Equation 6: first term is the score contribution; the second term is
  # the adjustment required for the robust subject-level residual.
  subject_group_risk <- (1 - Z) %o% risk_by_group[1, ] + Z %o% risk_by_group[2, ]
  denominator <- if (variance_method == "zhao_eq6") subject_group_risk else {
    matrix(risk, nrow = length(recs), ncol = length(times), byrow = TRUE)
  }
  adjustment_fraction <- cc$R_i / denominator
  adjustment_fraction[!is.finite(adjustment_fraction)] <- 0
  adjustment <- rowSums(
    sweep(adjustment_fraction, 2, dN * weight, "*") * centered_group
  )
  residual <- score_contribution - adjustment
  variance <- sum(residual^2)
  se <- sqrt(variance)
  z <- if (se > 0) W / se else NA_real_
  p <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
    greater = stats::pnorm(z, lower.tail = FALSE),
    less = stats::pnorm(z)
  )

  list(
    method = paste0("Zhao et al. extended ", test, " rank test"),
    test = test,
    variance_method = variance_method,
    alternative = alternative,
    group0 = group0,
    group1 = group1,
    n_subjects = length(recs),
    W = W,
    variance = variance,
    se = se,
    z = z,
    p.value = p,
    subject_residuals = data.frame(
      id = vapply(recs, function(r) as.character(r$id), character(1)),
      group = vapply(recs, function(r) as.character(r$group), character(1)),
      score_contribution = score_contribution,
      adjustment = adjustment,
      residual = residual
    ),
    detail = data.frame(
      time = times, risk = risk, risk_group1 = risk1, dN = dN,
      dLambda = d_lambda, survival = survival, weight = weight
    )
  )
}

#' Simulate Zhao et al. null-size experiments
#'
#' Generates recurrent gap-time data under the simulation design in Section 3
#' of Zhao et al. (2020) and repeatedly applies [zhao_rank_test()]. This is a
#' compact, reproducible diagnostic for nominal type-I-error calibration; it is
#' not intended to reproduce the paper's 100,000-replicate tables during a
#' routine package test.
#'
#' @param n_sim Number of simulated datasets.
#' @param n_per_group Subjects per group.
#' @param scenario `"paper"` reproduces the broad design choices in Zhao et
#'   al.'s Section 3. `"psg_200"` represents a fixed 8-hour PSG recording in
#'   minutes, calibrated for approximately `target_events` completed events per
#'   participant.
#' @param distribution Baseline gap-time distribution used in Zhao et al.'s
#'   Table 1. The paper's parameterization is used: exponential rate
#'   `exp(-4)`; Weibull shape 2 and scale `exp(4)` (equivalent to
#'   `lambda = exp(-8)` in `exp(-lambda * t^2)`); log-normal `meanlog = 4`,
#'   `sdlog = 0.5`; and log-logistic `lambda = exp(-4)`, `gamma = 0.5`.
#' @param heterogeneity Two endpoints of the subject multiplier `Z ~ U(a, b)`.
#' @param followup Administrative study end time.
#' @param target_events Target mean number of completed events per participant
#'   for `scenario = "psg_200"`.
#' @param time_resolution Optional gap-time discretization. The PSG default is
#'   0.5 minutes (30-second scoring epochs); leave `NULL` for continuous times.
#' @param alpha Nominal test level.
#' @param test Passed to [zhao_rank_test()].
#' @param variance_method Passed to [zhao_rank_test()].
#' @param alternative Passed to [zhao_rank_test()]. Zhao et al.'s Tables 2--3
#'   use the one-sided `"less"` alternative when `group = "treatment"` has
#'   longer gaps (a lower recurrent-event hazard).
#' @param treatment_time_multiplier Multiplicative factor applied to treatment
#'   gap times. The default `1` generates the null; values above one generate
#'   longer treatment gaps for power simulations.
#' @param seed Optional random seed.
#'
#' @return A list containing the empirical rejection rate, its Monte-Carlo
#' standard error, normal-score moments, p-values, and simulation settings.
#'
#' @references Zhao et al. (2020), Section 3 and Tables 1--3.
#'
#' @examples
#' \donttest{
#' zhao_rank_simulation(n_sim = 100, n_per_group = 50, seed = 1)
#' }
#'
#' @export
zhao_rank_simulation <- function(
    n_sim = 1000, n_per_group = 100,
    scenario = c("paper", "psg_200"),
    distribution = c("exponential", "weibull", "lognormal", "loglogistic"),
    heterogeneity = NULL, followup = NULL, target_events = 200,
    time_resolution = NULL,
    alpha = 0.05, test = c("logrank", "gehan_breslow", "peto_prentice"),
    variance_method = c("pooled_risk", "zhao_eq6"),
    alternative = c("two.sided", "greater", "less"),
    treatment_time_multiplier = 1, seed = NULL) {
  scenario <- match.arg(scenario)
  distribution <- match.arg(distribution)
  test <- match.arg(test)
  variance_method <- match.arg(variance_method)
  alternative <- match.arg(alternative)
  if (is.null(heterogeneity)) heterogeneity <- if (scenario == "paper") c(0.1, 1.9) else c(0.75, 1.25)
  if (is.null(followup)) followup <- if (scenario == "paper") 180 else 8 * 60
  if (is.null(time_resolution) && scenario == "psg_200") time_resolution <- 0.5
  if (!is.null(seed)) set.seed(seed)
  if (length(n_sim) != 1 || n_sim < 1 || n_sim != as.integer(n_sim)) {
    stop("n_sim must be a positive integer")
  }
  if (length(n_per_group) != 1 || n_per_group < 2 || n_per_group != as.integer(n_per_group)) {
    stop("n_per_group must be an integer >= 2")
  }
  if (length(heterogeneity) != 2 || heterogeneity[1] <= 0 || heterogeneity[2] < heterogeneity[1]) {
    stop("heterogeneity must contain positive endpoints in increasing order")
  }
  if (!is.numeric(followup) || length(followup) != 1 || followup <= 0) stop("followup must be > 0")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1)")
  if (!is.numeric(treatment_time_multiplier) || length(treatment_time_multiplier) != 1 || treatment_time_multiplier <= 0) {
    stop("treatment_time_multiplier must be > 0")
  }
  if (!is.null(time_resolution) && (!is.numeric(time_resolution) || length(time_resolution) != 1 || time_resolution <= 0)) {
    stop("time_resolution must be NULL or > 0")
  }

  if (scenario == "psg_200" && (!is.numeric(target_events) || length(target_events) != 1 || target_events <= 0)) {
    stop("target_events must be > 0 for scenario = 'psg_200'")
  }
  mean_gap <- NULL
  # Use the analytic expectation rather than consuming random-number state.
  if (scenario == "psg_200") {
    a <- heterogeneity[1]; b <- heterogeneity[2]
    inv_z_mean <- if (a == b) 1 / a else log(b / a) / (b - a)
    mean_gap <- followup * inv_z_mean / target_events
  }
  p_value <- z_value <- mean_events <- rep(NA_real_, n_sim)
  for (b in seq_len(n_sim)) {
    dat <- .zhao_simulated_gap_data(n_per_group, distribution, heterogeneity, followup,
      staggered_entry = scenario == "paper", mean_gap = mean_gap,
      time_resolution = time_resolution,
      treatment_time_multiplier = treatment_time_multiplier)
    fit <- zhao_rank_test(dat, group = "group", episode = "episode", test = test,
      variance_method = variance_method, alternative = alternative)
    p_value[b] <- fit$p.value
    z_value[b] <- fit$z
    mean_events[b] <- mean(tapply(dat$status, dat$id, sum))
  }
  reject <- p_value < alpha
  list(
    method = "Zhao et al. recurrent gap-time null simulation",
    empirical_type1_error = mean(reject, na.rm = TRUE),
    monte_carlo_se = sqrt(mean(reject, na.rm = TRUE) * (1 - mean(reject, na.rm = TRUE)) / sum(!is.na(reject))),
    mean_z = mean(z_value, na.rm = TRUE),
    sd_z = stats::sd(z_value, na.rm = TRUE),
    mean_completed_events = mean(mean_events),
    completed_events = mean_events,
    p.value = p_value,
    z = z_value,
    settings = list(n_sim = n_sim, n_per_group = n_per_group,
      scenario = scenario, distribution = distribution, heterogeneity = heterogeneity,
      followup = followup, target_events = if (scenario == "psg_200") target_events else NULL,
      mean_gap = mean_gap, time_resolution = time_resolution, alpha = alpha,
      test = test, variance_method = variance_method, alternative = alternative,
      treatment_time_multiplier = treatment_time_multiplier)
  )
}

.zhao_simulated_gap_data <- function(n_per_group, distribution, heterogeneity,
                                     followup, staggered_entry = TRUE, mean_gap = NULL,
                                     time_resolution = NULL,
                                     treatment_time_multiplier = 1) {
  n <- 2L * n_per_group
  group <- rep(c("control", "treatment"), each = n_per_group)
  multiplier <- stats::runif(n, heterogeneity[1], heterogeneity[2])
  available <- if (staggered_entry) followup - stats::runif(n, 0, followup) else rep(followup, n)
  pieces <- lapply(seq_len(n), function(i) {
    total <- 0
    gaps <- numeric()
    repeat {
      base_gap <- .zhao_draw_gap(distribution, mean_gap, time_resolution)
      next_gap <- base_gap * multiplier[i] * if (group[i] == "treatment") treatment_time_multiplier else 1
      if (total + next_gap >= available[i]) {
        gaps <- c(gaps, available[i] - total)
        status <- c(rep(1L, length(gaps) - 1L), 0L)
        break
      }
      gaps <- c(gaps, next_gap)
      total <- total + next_gap
    }
    data.frame(id = i, episode = seq_along(gaps), gap = gaps,
      status = status, group = group[i])
  })
  do.call(rbind, pieces)
}

.zhao_draw_gap <- function(distribution, mean_gap = NULL, time_resolution = NULL) {
  if (is.null(mean_gap)) {
    out <- switch(distribution,
      exponential = stats::rexp(1, rate = exp(-4)),
      weibull = stats::rweibull(1, shape = 2, scale = exp(4)),
      lognormal = stats::rlnorm(1, meanlog = 4, sdlog = 0.5),
      # Zhao Table 1: S(t) = 1 / (1 + (lambda*t)^(1/gamma)),
      # lambda = exp(-4), gamma = 0.5; invert the survival function.
      loglogistic = exp(4) * (stats::runif(1)^(-1) - 1)^0.5
    )
    return(if (is.null(time_resolution)) out else pmax(time_resolution, round(out / time_resolution) * time_resolution))
  }
  out <- switch(distribution,
    exponential = stats::rexp(1, rate = 1 / mean_gap),
    weibull = stats::rweibull(1, shape = 2, scale = mean_gap / gamma(1.5)),
    lognormal = stats::rlnorm(1, meanlog = log(mean_gap) - 0.5 * 0.5^2, sdlog = 0.5),
    loglogistic = {
      shape <- 0.5
      scale <- mean_gap * sin(pi * shape) / (shape * pi)
      scale * (stats::runif(1)^(-1) - 1)^shape
    }
  )
  if (is.null(time_resolution)) out else pmax(time_resolution, round(out / time_resolution) * time_resolution)
}
