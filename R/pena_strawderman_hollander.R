#' Pena--Strawderman--Hollander recurrent survival estimator
#'
#' Estimates the recurrent gap-time survival distribution using the generalized
#' product-limit estimator of Pena, Strawderman, and Hollander (PSH). With the
#' long gap-time data convention used in this package, PSH is the ordinary
#' Kaplan--Meier product limit formed from all observed gap records. Thus each
#' recurrent gap has equal weight, unlike [wc_surv()], which weights each
#' subject's eligible gaps by `1 / m_i*`.
#'
#' @param data A data frame in long gap-time format with one row per observed
#'   gap.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column.
#' @param episode Optional character scalar naming the within-subject episode
#'   order column.
#' @param tau Optional upper truncation time for the recurrent gap distribution.
#'
#' @return A data frame with one row per distinct completed gap time and columns
#'   `time`, `risk`, `dN`, `dLambda`, `surv_left`, and `surv`.
#'
#' @details
#' At a gap time `t`, PSH uses `R(t) = sum(i,j) I(Y_ij >= t)` and
#' `dN(t) = sum(i,j) I(Y_ij = t, delta_ij = 1)`, then applies
#' `prod(u <= t) {1 - dN(u) / R(u)}`. The renewal/IID interpretation of this
#' estimator differs from the marginal recurrent gap-time estimand of
#' Wang--Chang.
#'
#' @references
#' Pena EA, Strawderman RL, Hollander M (2001). Nonparametric estimation with
#' recurrent event data. *Journal of the American Statistical Association*,
#' 96, 1299-1315. \doi{10.1198/016214501753381887}.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:3, each = 3), episode = rep(1:3, 3),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5),
#'   status = rep(c(1, 1, 0), 3)
#' )
#' psh_surv(dat, episode = "episode")
#'
#' @export
psh_surv <- function(data, id = "id", gap = "gap", status = "status",
                     episode = NULL, tau = Inf) {
  d <- .validate_gap_data(data, id, gap, status, episode)
  times <- sort(unique(as.numeric(d[[gap]][d[[status]] == 1 & d[[gap]] <= tau])))
  if (!length(times)) {
    return(data.frame(
      time = numeric(), risk = numeric(), dN = numeric(), dLambda = numeric(),
      surv_left = numeric(), surv = numeric()
    ))
  }

  gap_value <- as.numeric(d[[gap]])
  event <- as.integer(d[[status]])
  risk <- vapply(times, function(tt) sum(gap_value >= tt), numeric(1))
  dN <- vapply(times, function(tt) sum(gap_value == tt & event == 1), numeric(1))
  dLambda <- ifelse(risk > 0, dN / risk, 0)
  surv_left <- c(1, utils::head(cumprod(pmax(0, 1 - dLambda)), -1))
  surv <- cumprod(pmax(0, 1 - dLambda))
  data.frame(
    time = times, risk = risk, dN = dN, dLambda = dLambda,
    surv_left = surv_left, surv = surv
  )
}
