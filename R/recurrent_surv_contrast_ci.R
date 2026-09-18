#' Subject-bootstrap confidence intervals for two recurrent survival curves
#'
#' Compares two WC or PSH curves at requested gap times using a cluster
#' bootstrap. Independent arms are resampled within arm; when `pair` is given,
#' complete paired clusters are resampled, preserving the paired visits.
#'
#' @param data Long recurrent gap-time data.
#' @param group Subject-level two-group column.
#' @param pair Optional paired-cluster column.
#' @param estimator,id,gap,status,episode,tau Passed to curve estimators.
#' @param times Evaluation times; pooled observed event times by default.
#' @param contrast `"difference"`, `"ratio"`, or `"log_ratio"`, defined as
#'   group 1 minus/divided by group 0.
#' @param interval Pointwise percentile intervals or simultaneous max-t bands.
#' @param conf.level,B,seed Usual bootstrap controls.
#' @return Observed group curves, a contrast curve with limits, and bootstrap diagnostics.
#' @export
recurrent_surv_contrast_ci <- function(data, group, pair = NULL,
    estimator = c("wc", "psh"), id = "id", gap = "gap", status = "status",
    episode = NULL, tau = Inf, times = NULL,
    contrast = c("difference", "ratio", "log_ratio"),
    interval = c("pointwise", "simultaneous"), conf.level = .95,
    B = 999, seed = NULL) {
  estimator <- match.arg(estimator); contrast <- match.arg(contrast); interval <- match.arg(interval)
  if (!group %in% names(data)) stop("Missing group column: ", group)
  if (!is.null(pair) && !pair %in% names(data)) stop("Missing pair column: ", pair)
  if (!is.numeric(B) || length(B) != 1 || B < 100 || B != as.integer(B)) stop("Use an integer B >= 100")
  d <- .validate_gap_data(data, id, gap, status, episode)
  subj <- unique(d[, c(id, group, if (!is.null(pair)) pair), drop = FALSE])
  lev <- unique(as.character(subj[[group]])); if (length(lev) != 2) stop("group must have exactly two levels")
  fun <- if (estimator == "wc") wc_surv else psh_surv
  get_curves <- function(dd) {
    lapply(lev, function(g) fun(dd[as.character(dd[[group]]) == g, , drop = FALSE], id, gap, status, episode, tau))
  }
  fits <- get_curves(d)
  if (is.null(times)) times <- sort(unique(unlist(lapply(fits, `[[`, "time"))))
  if (!length(times)) stop("No completed recurrent gaps within tau")
  calc <- function(ff) {
    a <- .recurrent_surv_step(ff[[1]], times); b <- .recurrent_surv_step(ff[[2]], times)
    if (contrast == "difference") b - a else if (contrast == "ratio") ifelse(a > 0, b / a, NA_real_) else ifelse(a > 0 && b > 0, log(b / a), NA_real_)
  }
  obs <- calc(fits); if (!is.null(seed)) set.seed(seed)
  boot <- matrix(NA_real_, B, length(times)); failures <- 0L
  for (bb in seq_len(B)) {
    if (is.null(pair)) {
      pieces <- lapply(lev, function(g) .cluster_bootstrap_df(d[as.character(d[[group]]) == g, , drop = FALSE], id))
      db <- do.call(rbind, pieces)
    } else {
      pp <- unique(as.character(d[[pair]])); draw <- sample(pp, length(pp), TRUE)
      db <- do.call(rbind, lapply(seq_along(draw), function(j) { z <- d[as.character(d[[pair]]) == draw[j], , drop = FALSE]; z[[id]] <- paste0(j, "_", z[[id]]); z }))
    }
    ans <- try(calc(get_curves(db)), silent = TRUE)
    if (inherits(ans, "try-error") || any(!is.finite(ans))) failures <- failures + 1L else boot[bb, ] <- ans
  }
  boot <- boot[stats::complete.cases(boot), , drop = FALSE]; if (nrow(boot) < 100) stop("Too few successful bootstrap resamples")
  a <- (1-conf.level)/2; se <- apply(boot, 2, stats::sd)
  if (interval == "pointwise") { lo <- apply(boot,2,stats::quantile,probs=a,type=8); hi <- apply(boot,2,stats::quantile,probs=1-a,type=8); crit <- NA_real_ } else { ok <- se > sqrt(.Machine$double.eps); z <- sweep(sweep(boot[,ok,drop=FALSE],2,colMeans(boot[,ok,drop=FALSE]),"-"),2,se[ok],"/"); crit <- unname(stats::quantile(apply(abs(z),1,max),conf.level,type=8)); lo <- obs; hi <- obs; lo[ok] <- obs[ok]-crit*se[ok]; hi[ok] <- obs[ok]+crit*se[ok] }
  list(method=paste("subject-bootstrap", interval, contrast, "of recurrent survival curves"), groups=lev, estimator=estimator, contrast=contrast, interval=interval, conf.level=conf.level, fits=stats::setNames(fits,lev), curve=data.frame(time=times, estimate=obs,std.error=se,lower=lo,upper=hi), critical.value=crit, bootstrap_success=nrow(boot), bootstrap_failures=failures)
}
