#' Confidence intervals for recurrent-survival quantiles
#'
#' Estimates a quantile of a [wc_surv()] Wang--Chang marginal gap-time survival
#' curve or a [psh_surv()] Pena--Strawderman--Hollander curve. For example,
#' `survival_prob = 0.5` estimates the median gap time,
#' `inf {t: S(t) <= 0.5}`.
#'
#' @param data A data frame in the long gap-time format used by [wc_surv()] and
#'   [psh_surv()].
#' @param estimator Either `"wc"` for Wang--Chang or `"psh"` for PSH/pooled KM.
#' @param id,gap,status,episode,tau Passed to the selected curve estimator.
#' @param survival_prob Survival probability defining the quantile. `0.5` is
#'   the median; `0.25` is the 75th percentile of the gap-time distribution.
#' @param conf.level Confidence level in `(0, 1)`.
#' @param B Number of subject-bootstrap resamples. The default is 999; use
#'   2,000 or more for stable two-sided 95 percent tail limits when computation
#'   permits.
#' @param seed Optional random-number seed.
#' @param keep_bootstrap Logical; retain bootstrap quantiles when available.
#'
#' @return A list containing the quantile `estimate`, its `lower` and `upper`
#' confidence limits, the fitted `curve`, method details, and bootstrap
#' diagnostics. An infinite estimate or upper limit means that the fitted curve
#' did not reach `survival_prob` before the available support.
#'
#' @details
#' The bootstrap resamples complete participant histories and refits the curve.
#' It therefore preserves recurrent-gap dependence and the event-count weights
#' of the WC estimator. Independence-style normal curve limits are deliberately
#' not inverted for a quantile interval because they do not account for
#' within-subject clustering.
#'
#' A quantile interval, a pointwise curve interval, and a simultaneous curve
#' band have different coverage targets. Use [recurrent_surv_ci()] for the
#' latter two curve-level quantities.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   id = rep(1:4, each = 3), episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4)
#' )
#' recurrent_surv_quantile_ci(dat, estimator = "wc", episode = "episode",
#'                           survival_prob = 0.5, B = 999, seed = 1)
#' }
#'
#' @export
recurrent_surv_quantile_ci <- function(
    data,
    estimator = c("wc", "psh"),
    id = "id",
    gap = "gap",
    status = "status",
    episode = NULL,
    tau = Inf,
    survival_prob = 0.5,
    conf.level = 0.95,
    B = 999,
    seed = NULL,
    keep_bootstrap = FALSE) {
  estimator <- match.arg(estimator)
  if (!is.numeric(survival_prob) || length(survival_prob) != 1L ||
      !is.finite(survival_prob) || survival_prob <= 0 || survival_prob >= 1) {
    stop("survival_prob must be one finite number strictly between 0 and 1")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be one finite number strictly between 0 and 1")
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 100 ||
      B != as.integer(B)) {
    stop("Use an integer B >= 100; default B = 999, and B >= 2000 is recommended for stable 95% tails")
  }

  d <- .validate_gap_data(data, id, gap, status, episode)
  curve_fun <- if (estimator == "wc") wc_surv else psh_surv
  fit <- curve_fun(d, id = id, gap = gap, status = status,
                   episode = episode, tau = tau)
  estimate <- .recurrent_surv_quantile(fit, survival_prob)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  boot_q <- rep(NA_real_, B)
  failures <- 0L
  for (b in seq_len(B)) {
    db <- .cluster_bootstrap_df(d, id)
    ans <- try({
      fb <- curve_fun(db, id = id, gap = gap, status = status,
                      episode = episode, tau = tau)
      .recurrent_surv_quantile(fb, survival_prob)
    }, silent = TRUE)
    if (inherits(ans, "try-error") || length(ans) != 1L || is.na(ans)) {
      failures <- failures + 1L
    } else {
      boot_q[b] <- ans
    }
  }
  boot_q <- boot_q[!is.na(boot_q)]
  if (length(boot_q) < 100L) {
    stop("Too few successful bootstrap resamples")
  }
  alpha <- 1 - conf.level
  out <- list(
    method = paste(
      if (estimator == "wc") "Wang-Chang" else "PSH / pooled Kaplan-Meier",
      "subject-bootstrap percentile quantile interval"
    ),
    estimator = estimator,
    survival_prob = survival_prob,
    conf.level = conf.level,
    estimate = estimate,
    lower = as.numeric(stats::quantile(boot_q, alpha / 2, names = FALSE,
                                        type = 1)),
    upper = as.numeric(stats::quantile(boot_q, 1 - alpha / 2, names = FALSE,
                                        type = 1)),
    fit = fit,
    bootstrap_success = length(boot_q),
    bootstrap_failures = failures,
    n_infinite = sum(is.infinite(boot_q)),
    prop_infinite = mean(is.infinite(boot_q)),
    warning = if (any(is.infinite(boot_q))) {
      "Some bootstrap curves did not reach survival_prob; consider a restricted mean gap-time analysis."
    } else NULL
  )
  if (keep_bootstrap) {
    out$bootstrap_quantile <- boot_q
  }
  out
}

#' Extract a recurrent-survival quantile from a product-limit curve
#' @keywords internal
#' @noRd
.recurrent_surv_quantile <- function(fit, survival_prob) {
  .recurrent_surv_quantile_from_values(fit$time, fit$surv, survival_prob)
}

#' Extract a quantile from right-continuous survival step values
#' @keywords internal
#' @noRd
.recurrent_surv_quantile_from_values <- function(time, surv, survival_prob) {
  if (!length(time)) {
    return(Inf)
  }
  hit <- which(surv <= survival_prob)
  if (!length(hit)) Inf else time[hit[1]]
}
