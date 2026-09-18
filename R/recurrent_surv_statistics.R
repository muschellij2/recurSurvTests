#' Restricted mean recurrent gap time with a subject-bootstrap confidence interval
#'
#' Estimates the area under a Wang--Chang or PSH recurrent gap-time survival
#' curve through `tau`, and obtains a percentile confidence interval by
#' resampling complete subject histories.
#'
#' @param data,id,gap,status,episode Passed to [wc_surv()] or [psh_surv()].
#' @param estimator Either `"wc"` or `"psh"`.
#' @param tau Positive restriction time.
#' @param conf.level Confidence level.
#' @param B Number of subject bootstrap resamples; 999 is the default and 2000
#'   or more is recommended for stable two-sided 95 percent tails.
#' @param seed Optional random seed.
#' @return A list with the restricted mean `estimate`, percentile limits, fitted
#'   curve, and bootstrap diagnostics.
#' @details The restricted mean gap time is `integral_0^tau S(u) du`. It remains
#' estimable when a survival median is not reached. Resampling is by subject,
#' preserving recurrent-gap clustering and WC weights.
#' @export
recurrent_surv_rmst_ci <- function(data, estimator = c("wc", "psh"),
                                   id = "id", gap = "gap", status = "status",
                                   episode = NULL, tau, conf.level = .95,
                                   B = 999, seed = NULL) {
  estimator <- match.arg(estimator)
  if (missing(tau) || !is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) stop("tau must be one positive finite number")
  if (!is.numeric(B) || length(B) != 1 || B < 100 || B != as.integer(B)) stop("Use an integer B >= 100")
  d <- .validate_gap_data(data, id, gap, status, episode)
  fun <- if (estimator == "wc") wc_surv else psh_surv
  fit <- fun(d, id, gap, status, episode, tau)
  est <- .recurrent_surv_rmst(fit, tau)
  if (!is.null(seed)) set.seed(seed)
  bb <- vapply(seq_len(B), function(b) {
    fb <- fun(.cluster_bootstrap_df(d, id), id, gap, status, episode, tau)
    .recurrent_surv_rmst(fb, tau)
  }, numeric(1))
  a <- (1 - conf.level) / 2
  list(method = paste(if (estimator == "wc") "Wang-Chang" else "PSH", "subject-bootstrap restricted mean gap time"),
       estimator = estimator, tau = tau, conf.level = conf.level, estimate = est,
       lower = unname(stats::quantile(bb, a, type = 8)),
       upper = unname(stats::quantile(bb, 1 - a, type = 8)), fit = fit,
       bootstrap_success = length(bb), bootstrap_failures = 0L)
}

#' @keywords internal
#' @noRd
.recurrent_surv_rmst <- function(fit, tau) {
  if (!nrow(fit)) return(tau)
  tt <- c(0, fit$time[fit$time < tau], tau)
  ss <- .recurrent_surv_step(fit, tt[-length(tt)])
  sum(diff(tt) * ss)
}

#' Plot recurrent-survival confidence intervals or bands
#'
#' @param x Object returned by [recurrent_surv_ci()].
#' @param ... Passed to [plot()].
#' @return Invisibly returns `x`.
#' @export
plot_recurrent_surv_ci <- function(x, ...) {
  stopifnot(is.list(x), is.data.frame(x$curve))
  z <- x$curve
  label <- if (identical(x$interval, "simultaneous")) "simultaneous confidence band" else "pointwise confidence interval"
  graphics::plot(z$time, z$surv, type = "s", ylim = c(0, 1), xlab = "Gap time", ylab = "Survival", ...)
  graphics::lines(z$time, z$lower, type = "s", lty = 2)
  graphics::lines(z$time, z$upper, type = "s", lty = 2)
  graphics::legend("topright", legend = paste0(round(100 * x$conf.level), "% ", label), bty = "n")
  invisible(x)
}
