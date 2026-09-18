#' Bootstrap equal-cause test for recurrent cumulative-incidence functions
#'
#' Constructs integrated cause-specific RCIF contrasts and estimates their
#' covariance by resampling subjects with replacement. The resulting quadratic
#' form is a transparent Wald reconstruction of the equal-cause test described
#' by Sivadasan and Sankaran.
#'
#' @param data A data frame in long recurrent gap-time format.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column.
#' @param cause Character scalar naming the recurrent-event cause column.
#' @param episode Optional character scalar naming the within-subject episode
#'   order column.
#' @param subject_weight Passed to [ss_rcif()].
#' @param cause_weighting Passed to [ss_rcif()].
#' @param survival_side Passed to [ss_rcif()].
#' @param tau Optional upper truncation time.
#' @param weight Either `NULL` for unit weights, a function of the fitted RCIF
#'   curve returning one weight per event time, or a numeric vector with one
#'   value per event time. A function is recommended for bootstrap analyses
#'   because bootstrap event-time grids can vary.
#' @param B Number of subject-level bootstrap replicates. At least 20 are
#'   required; 999 or more are recommended for analysis.
#' @param seed Optional random-number seed.
#'
#' @return A list containing the observed integrated contrasts, Helmert contrast
#'   vector, Wald chi-square statistic, effective degrees of freedom, asymptotic
#'   and bootstrap-calibrated p-values, bootstrap diagnostics, covariance matrix,
#'   and the fitted RCIF object.
#'
#' @details
#' The null hypothesis is that the cause-specific recurrent cumulative-incidence
#' functions are equal, equivalently `F_l(t) = F(t) / k` for every cause under
#' the formulation used here.
#'
#' This implementation should be treated as a **reconstruction**, not as an
#' exact reproduction of unpublished author software. The published appendix
#' motivates integrated contrasts of `F_l - F/k`, while the covariance here is
#' estimated transparently by subject-level bootstrap.
#'
#' @references
#' Sivadasan SM, Sankaran PG (2022). A nonparametric test for comparing recurrent
#' cumulative incidence functions. *METRON*. \doi{10.1007/s40300-022-00228-x}.
#'
#' Sivadasan SM, Sankaran PG (2023). Nonparametric estimation of cumulative
#' incidence functions of recurrent events. *Statistica*, 83(1), 3-25.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   id = rep(1:4, each = 3),
#'   episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4),
#'   cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
#' )
#' ss_rcif_equal_causes_test(dat, episode = "episode", B = 999, seed = 1)
#' }
#'
#' @export
ss_rcif_equal_causes_test <- function(
    data,
    id = "id",
    gap = "gap",
    status = "status",
    cause = "cause",
    episode = NULL,
    subject_weight = c("one", "followup"),
    cause_weighting = c("shared", "literal_eq15"),
    survival_side = c("right", "left"),
    tau = Inf,
    weight = NULL,
    B = 999,
    seed = NULL) {
  if (B < 20) {
    stop("Use B >= 20; B=999 or more is recommended")
  }

  cause_weighting <- match.arg(cause_weighting)
  survival_side <- match.arg(survival_side)

  d <- .validate_gap_data(data, id, gap, status, episode)
  if (!cause %in% names(d)) {
    stop("Missing cause column: ", cause)
  }

  fit <- ss_rcif(
    d,
    id = id,
    gap = gap,
    status = status,
    cause = cause,
    episode = episode,
    subject_weight = subject_weight,
    cause_weighting = cause_weighting,
    survival_side = survival_side,
    tau = tau
  )

  v <- .integrated_rcif_contrast(fit, weight)
  levels <- fit$causes
  k <- length(levels)
  C <- .helmert_rows(k)
  q <- as.vector(C %*% v)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  boot_v <- matrix(NA_real_, B, k)
  failures <- 0L

  for (b in seq_len(B)) {
    db <- .cluster_bootstrap_df(d, id)
    ans <- try({
      fb <- ss_rcif(
        db,
        id = id,
        gap = gap,
        status = status,
        cause = cause,
        episode = episode,
        subject_weight = subject_weight,
        cause_weighting = cause_weighting,
        survival_side = survival_side,
        tau = tau
      )

      if (!setequal(as.character(fb$causes), as.character(levels))) {
        stop("cause absent")
      }

      if (!is.null(weight) && !is.function(weight) &&
          length(weight) != nrow(fb$curve)) {
        stop("numeric weight incompatible with bootstrap grid; use callable weight")
      }

      vb0 <- .integrated_rcif_contrast(fb, weight)
      names(vb0) <- as.character(fb$causes)
      as.numeric(vb0[as.character(levels)])
    }, silent = TRUE)

    if (inherits(ans, "try-error")) {
      failures <- failures + 1L
    } else {
      boot_v[b, ] <- ans
    }
  }

  good <- apply(boot_v, 1, function(x) all(is.finite(x)))
  boot_v <- boot_v[good, , drop = FALSE]

  if (nrow(boot_v) < max(20, k + 5)) {
    stop("Too few successful bootstrap replicates")
  }

  Sigma_v <- stats::cov(boot_v)
  Sigma_q <- C %*% Sigma_v %*% t(C)
  pi <- .pinv_rank(Sigma_q)

  if (!pi$rank) {
    stat <- NA_real_
    pchi <- NA_real_
    pboot <- NA_real_
  } else {
    stat <- as.numeric(t(q) %*% pi$inv %*% q)
    pchi <- stats::pchisq(stat, df = pi$rank, lower.tail = FALSE)

    center <- sweep(boot_v, 2, colMeans(boot_v), "-")
    qb <- center %*% t(C)
    bstat <- rowSums((qb %*% pi$inv) * qb)
    pboot <- (1 + sum(bstat >= stat)) / (length(bstat) + 1)
  }

  list(
    method = paste(
      "Sivadasan/Sankaran METRON RCIF equal-causes test",
      "(bootstrap reconstruction)"
    ),
    null = "F_1(t)=...=F_k(t)=F(t)/k",
    logrank_analogue = FALSE,
    v = setNames(v, as.character(levels)),
    contrast = q,
    chisq = stat,
    df = pi$rank,
    p.value.chisq = pchi,
    p.value.bootstrap.wald = pboot,
    bootstrap_success = nrow(boot_v),
    bootstrap_failures = failures,
    cov_v = Sigma_v,
    fit = fit
  )
}
