#' Luo-Huang weighted recurrent gap-time rank test
#'
#' Implements the two-sample `G_rho*` test of Luo and Huang using the
#' weighted-risk-set representation of recurrent gap-time data. With `rho = 0`,
#' the statistic is the recurrent gap-time analogue of the ordinary log-rank
#' test. With `rho = 1`, it is an analogue of the Peto-Prentice generalized
#' Wilcoxon test.
#'
#' The variance is estimated from subject-level martingale-style contributions,
#' retaining the within-subject dependence among recurrent gaps.
#'
#' @param data A data frame in long recurrent gap-time format.
#' @param group Character scalar naming a subject-level two-group variable.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column.
#' @param episode Optional character scalar naming the within-subject episode
#'   order column.
#' @param rho Nonnegative rank-weight parameter. `rho = 0` gives the recurrent
#'   log-rank analogue; `rho = 1` gives the recurrent Peto-Prentice analogue.
#' @param tau Optional upper truncation time.
#'
#' @return A list containing the test statistic, standard error, `z` statistic,
#'   chi-square statistic, p-value, subject-level influence contributions, and
#'   event-time details.
#'
#' @details
#' The current implementation requires exactly two groups. Group ordering is
#' determined by first appearance among subjects; `group0` and `group1` are
#' returned explicitly so the score direction is interpretable.
#'
#' @references
#' Luo X, Huang CY (2011). Analysis of recurrent gap time data using the weighted
#' risk-set method and the modified within-cluster resampling method. *Statistics
#' in Medicine*, 30, 301-311. \doi{10.1002/sim.4074}.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:4, each = 3),
#'   episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4),
#'   arm = rep(c("A", "A", "B", "B"), each = 3)
#' )
#' wc_logrank(dat, group = "arm", episode = "episode", rho = 0)
#'
#' @export
wc_logrank <- function(data, group, id = "id", gap = "gap", status = "status",
                       episode = NULL, rho = 0, tau = Inf) {
  if (rho < 0) {
    stop("rho must be >= 0")
  }

  d <- .validate_gap_data(data, id, gap, status, episode)
  if (!group %in% names(d)) {
    stop("Missing group column: ", group)
  }

  recs <- .subject_records(d, id, gap, status, group = group)
  lev <- unique(vapply(recs, function(r) as.character(r$group), character(1)))
  if (length(lev) != 2) {
    stop("wc_logrank currently requires exactly two groups")
  }

  control <- lev[1]
  treatment <- lev[2]
  Z <- as.numeric(
    vapply(recs, function(r) as.character(r$group), character(1)) == treatment
  )

  times <- .event_grid(recs, tau)
  if (!length(times)) {
    stop("No completed recurrent gaps within tau")
  }

  cc <- .wrs_components(recs, times)
  risk <- colSums(cc$R_i)
  dN <- colSums(cc$dN_i)
  dLambda0 <- ifelse(risk > 0, dN / risk, 0)

  surv_left <- numeric(length(times))
  s <- 1
  for (j in seq_along(times)) {
    surv_left[j] <- s
    s <- s * max(0, 1 - dLambda0[j])
  }

  w <- surv_left^rho
  risk1 <- colSums(cc$R_i * Z)
  e <- ifelse(risk > 0, risk1 / risk, 0)

  zmat <- matrix(Z, nrow = length(Z), ncol = length(times))
  event_score <- rowSums(
    sweep(cc$dN_i, 2, w, "*") * sweep(zmat, 2, e, "-")
  )

  n <- length(recs)
  G <- mean(event_score)

  dM <- cc$dN_i - sweep(cc$R_i, 2, dLambda0, "*")
  psi <- rowSums(
    sweep(dM, 2, w, "*") * sweep(zmat, 2, e, "-")
  )

  var_sqrtn_G <- mean(psi^2)
  se <- sqrt(var_sqrtn_G / n)
  z <- if (se > 0) G / se else NA_real_
  chisq <- z^2
  p <- stats::pchisq(chisq, df = 1, lower.tail = FALSE)

  detail <- data.frame(
    time = times,
    risk = risk,
    risk_group1 = risk1,
    expected_group_fraction = e,
    dN = dN,
    dLambda0 = dLambda0,
    surv_left = surv_left,
    weight = w
  )

  list(
    method = "Luo-Huang WRS G_rho*",
    rho = rho,
    interpretation = if (rho == 0) {
      "recurrent gap-time log-rank analogue"
    } else if (rho == 1) {
      "recurrent Peto-Prentice analogue"
    } else {
      "weighted recurrent gap-time G_rho* test"
    },
    group0 = control,
    group1 = treatment,
    n_subjects = n,
    G = G,
    se = se,
    z = z,
    chisq = chisq,
    df = 1L,
    p.value = p,
    subject_influence = data.frame(
      id = vapply(recs, function(r) as.character(r$id), character(1)),
      psi = psi,
      event_score = event_score
    ),
    detail = detail
  )
}
