#' Subject-bootstrap confidence intervals for recurrent survival curves
#'
#' Forms confidence intervals for a [wc_surv()] Wang--Chang marginal
#' gap-time survival curve or a [psh_surv()] Pena--Strawderman--Hollander
#' (pooled Kaplan--Meier) curve. Entire subject histories, rather than
#' individual gap rows, are resampled. This preserves the within-subject
#' dependence and the event-count structure that row-level resampling would
#' destroy.
#'
#' @param data A data frame in the long gap-time format used by [wc_surv()] and
#'   [psh_surv()].
#' @param estimator Either `"wc"` for Wang--Chang or `"psh"` for PSH/pooled KM.
#' @param id,gap,status,episode,tau Passed to the selected curve estimator.
#' @param times Times at which intervals are evaluated. The observed curve's
#'   event times are used by default. The fitted curve is evaluated as a
#'   right-continuous step function at supplied times.
#' @param conf.level Confidence level in `(0, 1)`.
#' @param interval Either `"pointwise"` (the standard choice for a survival
#'   curve) or `"simultaneous"`. Pointwise intervals are cluster-bootstrap
#'   percentile intervals and have the stated coverage separately at each fixed
#'   time. Simultaneous intervals are max-t bootstrap bands intended to cover
#'   all supplied times jointly.
#' @param ci_method Either `"bootstrap"` (default) or `"normal"`. The normal
#'   option forms pointwise intervals from the independence-style
#'   Nelson--Aalen/Greenwood standard-error approximation. It is included to
#'   reproduce the type of pointwise interval obtainable from legacy fitters;
#'   it is not generally appropriate for dependent recurrent gaps. Simultaneous
#'   bands require `ci_method = "bootstrap"`.
#' @param B Number of subject-level bootstrap resamples. At least 100 are
#'   required for `ci_method = "bootstrap"`; 999 or more are recommended for
#'   final analysis. The default of 999 follows common bootstrap and permutation
#'   practice. For stable two-sided 95 percent percentile limits or max-t band
#'   critical values, use 2,000 or more resamples when computation permits. It
#'   is ignored for `ci_method = "normal"`.
#' @param seed Optional random-number seed.
#' @param keep_bootstrap Logical; retain the bootstrap survival matrix in the
#'   returned object. This can be large when many time points are requested.
#'
#' @return A list with the observed `fit`, a `curve` data frame containing
#'   `time`, `surv`, `std.error`, `lower`, and `upper`, bootstrap diagnostics,
#'   and (for simultaneous bands) the max-t critical value. `bootstrap_surv` is
#'   included only when `keep_bootstrap = TRUE`.
#'
#' @details
#' Pointwise and simultaneous coverage answer different questions. A 95 percent
#' pointwise interval has approximately 95 percent coverage at one pre-specified
#' gap time; it is not a 95 percent band for the whole curve. Use
#' `interval = "simultaneous"` when the inferential statement concerns all
#' supplied times at once. The max-t band uses the bootstrap standard deviation
#' at each time and is clipped to `[0, 1]`.
#'
#' The 2.5th and 97.5th percentile limits of a two-sided 95 percent interval
#' are estimated from roughly 25 tail resamples at `B = 999` and 50 at
#' `B = 2000`. Thus 999 is a useful default, while 2,000 or more reduces
#' Monte-Carlo variability in tail limits. A confidence interval for a survival
#' quantile or a between-group curve difference is a different inferential
#' target and is not returned by this function.
#'
#' The resampling unit is the subject. This is important for both estimators:
#' PSH's pooled product limit has correlated gap rows, and Wang--Chang's
#' subject-level weights depend on the complete recurrence history. Bootstrap
#' intervals are not a substitute for checking independent censoring and the
#' estimand-specific assumptions of the selected curve.
#'
#' @references
#' Wang MC, Chang SH (1999). Nonparametric estimation of a recurrent survival
#' function. *Journal of the American Statistical Association*, 94, 146-153.
#'
#' Pena EA, Strawderman RL, Hollander M (2001). Nonparametric estimation with
#' recurrent event data. *Journal of the American Statistical Association*,
#' 96, 1299-1315. \doi{10.1198/016214501753381887}.
#'
#' @examples
#' \donttest{
#' dat <- data.frame(
#'   id = rep(1:4, each = 3), episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4)
#' )
#' ci <- recurrent_surv_ci(dat, estimator = "wc", episode = "episode",
#'                         B = 999, seed = 1)
#' ci$curve
#' }
#'
#' @export
recurrent_surv_ci <- function(
    data,
    estimator = c("wc", "psh"),
    id = "id",
    gap = "gap",
    status = "status",
    episode = NULL,
    tau = Inf,
    times = NULL,
    conf.level = 0.95,
    interval = c("pointwise", "simultaneous"),
    ci_method = c("bootstrap", "normal"),
    B = 999,
    seed = NULL,
    keep_bootstrap = FALSE) {
  estimator <- match.arg(estimator)
  interval <- match.arg(interval)
  ci_method <- match.arg(ci_method)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be one finite number strictly between 0 and 1")
  }
  if (ci_method == "bootstrap" &&
      (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 100 ||
       B != as.integer(B))) {
    stop("Use an integer B >= 100; default B = 999, and B >= 2000 is recommended for stable 95% tails")
  }
  if (ci_method == "normal" && interval == "simultaneous") {
    stop("simultaneous intervals require ci_method = 'bootstrap'")
  }

  d <- .validate_gap_data(data, id, gap, status, episode)
  curve_fun <- if (estimator == "wc") wc_surv else psh_surv
  fit <- curve_fun(d, id = id, gap = gap, status = status,
                   episode = episode, tau = tau)
  if (is.null(times)) {
    times <- fit$time
  }
  times <- sort(unique(as.numeric(times)))
  if (!length(times) || any(!is.finite(times)) || any(times < 0) ||
      any(times > tau)) {
    stop("times must be nonempty, finite, nonnegative, and no greater than tau")
  }

  observed <- .recurrent_surv_step(fit, times)
  alpha <- 1 - conf.level
  if (ci_method == "normal") {
    # This is the same independence-style Nelson--Aalen/Greenwood form used
    # by legacy recurrent-curve fitters, evaluated at the requested step times.
    var_fit <- fit$surv^2 * cumsum(fit$dN / fit$risk^2)
    jj <- findInterval(times, fit$time)
    se <- numeric(length(times))
    se[jj > 0L] <- sqrt(var_fit[jj[jj > 0L]])
    critical <- stats::qnorm(1 - alpha / 2)
    lower <- pmax(0, observed - critical * se)
    upper <- pmin(1, observed + critical * se)
    boot <- NULL
    failures <- NA_integer_
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    boot <- matrix(NA_real_, nrow = B, ncol = length(times))
    failures <- 0L
    for (b in seq_len(B)) {
      db <- .cluster_bootstrap_df(d, id)
      ans <- try({
        fb <- curve_fun(db, id = id, gap = gap, status = status,
                        episode = episode, tau = tau)
        .recurrent_surv_step(fb, times)
      }, silent = TRUE)
      if (inherits(ans, "try-error") || length(ans) != length(times) ||
          any(!is.finite(ans))) {
        failures <- failures + 1L
      } else {
        boot[b, ] <- ans
      }
    }
    boot <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(boot) < 100L) {
      stop("Too few successful bootstrap resamples")
    }

    se <- apply(boot, 2, stats::sd)
    if (interval == "pointwise") {
      lower <- apply(boot, 2, stats::quantile, probs = alpha / 2,
                     names = FALSE, type = 8)
      upper <- apply(boot, 2, stats::quantile, probs = 1 - alpha / 2,
                     names = FALSE, type = 8)
      critical <- NA_real_
    } else {
      usable <- is.finite(se) & se > sqrt(.Machine$double.eps)
      if (!any(usable)) {
        stop("Bootstrap standard errors are zero at all requested times")
      }
      centered <- sweep(boot[, usable, drop = FALSE], 2,
                        colMeans(boot[, usable, drop = FALSE]), "-")
      zmax <- apply(abs(sweep(centered, 2, se[usable], "/")), 1, max)
      critical <- as.numeric(stats::quantile(zmax, probs = conf.level,
                                              names = FALSE, type = 8))
      lower <- observed
      upper <- observed
      lower[usable] <- pmax(0, observed[usable] - critical * se[usable])
      upper[usable] <- pmin(1, observed[usable] + critical * se[usable])
    }
  }

  out <- list(
    method = paste(
      if (estimator == "wc") "Wang-Chang" else "PSH / pooled Kaplan-Meier",
      if (ci_method == "bootstrap") "subject-bootstrap" else "normal approximation",
      interval,
      if (ci_method == "normal") "intervals" else if (interval == "pointwise") "percentile intervals" else "max-t band"
    ),
    estimator = estimator,
    interval = interval,
    ci_method = ci_method,
    conf.level = conf.level,
    fit = fit,
    curve = data.frame(time = times, surv = observed, std.error = se,
                       lower = lower, upper = upper),
    critical.value = critical,
    bootstrap_success = if (is.null(boot)) NA_integer_ else nrow(boot),
    bootstrap_failures = failures
  )
  if (keep_bootstrap && !is.null(boot)) {
    out$bootstrap_surv <- boot
  }
  out
}

#' Evaluate a recurrent survival fit as a right-continuous step function
#' @keywords internal
#' @noRd
.recurrent_surv_step <- function(fit, times) {
  if (!nrow(fit)) {
    return(rep(1, length(times)))
  }
  jj <- findInterval(times, fit$time)
  ans <- rep(1, length(times))
  has_event <- jj > 0L
  ans[has_event] <- fit$surv[jj[has_event]]
  ans
}
