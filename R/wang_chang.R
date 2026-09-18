#' Wang-Chang weighted-risk-set recurrent gap-time survival estimator
#'
#' Estimates the marginal survival distribution of recurrent gap times using the
#' subject-weighted risk-set representation that is algebraically equivalent to
#' the Wang-Chang estimator and used by Luo and Huang.
#'
#' Each subject contributes weight `1 / m_i*` to each eligible recurrent gap,
#' where `m_i*` is the number of completed gaps used by the weighted-risk-set
#' construction. This avoids allowing subjects with many recurrent events to
#' dominate the marginal gap-time survival estimate.
#'
#' @param data A data frame in long gap-time format with one row per observed
#'   gap.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column. Completed
#'   recurrent gaps are coded `1`; the final censored gap is coded `0`.
#' @param episode Optional character scalar naming the within-subject episode
#'   order column.
#' @param tau Optional upper truncation time for the recurrent gap distribution.
#'
#' @return A data frame with one row per distinct completed gap time and columns
#'   `time`, `risk`, `dN`, `dLambda`, `surv_left`, and `surv`.
#'
#' @references
#' Wang MC, Chang SH (1999). Nonparametric estimation of a recurrent survival
#' function. *Journal of the American Statistical Association*, 94, 146-153.
#'
#' Luo X, Huang CY (2011). Analysis of recurrent gap time data using the weighted
#' risk-set method and the modified within-cluster resampling method. *Statistics
#' in Medicine*, 30, 301-311. \doi{10.1002/sim.4074}.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:3, c(3, 2, 3)),
#'   gap = c(2, 3, 4, 1, 5, 2, 4, 3),
#'   status = c(1, 1, 0, 1, 0, 1, 1, 0),
#'   episode = c(1, 2, 3, 1, 2, 1, 2, 3)
#' )
#' wc_surv(dat, episode = "episode")
#'
#' @export
wc_surv <- function(data, id = "id", gap = "gap", status = "status",
                    episode = NULL, tau = Inf) {
  d <- .validate_gap_data(data, id, gap, status, episode)
  recs <- .subject_records(d, id, gap, status)
  times <- .event_grid(recs, tau)

  if (!length(times)) {
    return(data.frame(
      time = numeric(),
      risk = numeric(),
      dN = numeric(),
      dLambda = numeric(),
      surv_left = numeric(),
      surv = numeric()
    ))
  }

  cc <- .wrs_components(recs, times)
  risk <- colSums(cc$R_i)
  dN <- colSums(cc$dN_i)
  dLambda <- ifelse(risk > 0, dN / risk, 0)

  surv_left <- numeric(length(times))
  surv <- numeric(length(times))
  s <- 1
  for (j in seq_along(times)) {
    surv_left[j] <- s
    s <- s * max(0, 1 - dLambda[j])
    surv[j] <- s
  }

  data.frame(
    time = times,
    risk = risk,
    dN = dN,
    dLambda = dLambda,
    surv_left = surv_left,
    surv = surv
  )
}
