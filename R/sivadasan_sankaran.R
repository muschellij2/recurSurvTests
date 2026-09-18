#' Sivadasan-Sankaran recurrent cumulative-incidence estimator
#'
#' Implements a transparent reference version of the recurrent cumulative
#' incidence function (RCIF) estimator proposed by Sivadasan and Sankaran for
#' recurrent competing-risk gap-time data.
#'
#' @param data A data frame in long recurrent gap-time format.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column.
#' @param cause Character scalar naming the recurrent-event cause column. Causes
#'   should be non-missing on completed-event rows.
#' @param episode Optional character scalar naming the within-subject episode
#'   order column.
#' @param subject_weight Either `"one"`, giving each subject unit weight,
#'   `"followup"`, weighting subjects by total follow-up, or a function taking a
#'   subject-specific data frame and returning one nonnegative scalar.
#' @param cause_weighting Either `"shared"` (default), which uses the same
#'   subject-level `m_i*` denominator for pooled and cause-specific event
#'   processes, or `"literal_eq15"`, which reproduces the printed cause-specific
#'   denominator in Equation 15 of the 2023 paper.
#' @param survival_side Whether RCIF increments use the right-continuous
#'   (`"right"`) or left-limit (`"left"`) survival estimate at each event time.
#' @param tau Optional upper truncation time.
#'
#' @return A list containing the pooled recurrent survival process, cause-specific
#'   recurrent cumulative-incidence curves, subject-weighting metadata, and an
#'   additivity diagnostic.
#'
#' @details
#' The printed Equation 15 in the 2023 *Statistica* paper uses a cause-specific
#' denominator `m_il` and an indicator requiring at least two events of a cause.
#' Taken literally, those cause-specific increments need not sum to the pooled
#' event process. Therefore the default is `cause_weighting = "shared"`, while
#' `"literal_eq15"` is provided for literal reproduction and sensitivity
#' analysis.
#'
#' @references
#' Sivadasan SM, Sankaran PG (2023). Nonparametric estimation of cumulative
#' incidence functions of recurrent events. *Statistica*, 83(1), 3-25.
#'
#' Sivadasan SM, Sankaran PG (2022). A nonparametric test for comparing recurrent
#' cumulative incidence functions. *METRON*. \doi{10.1007/s40300-022-00228-x}.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:4, each = 3),
#'   episode = rep(1:3, 4),
#'   gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
#'   status = rep(c(1, 1, 0), 4),
#'   cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
#' )
#' ss_rcif(dat, episode = "episode")
#'
#' @export
ss_rcif <- function(data, id = "id", gap = "gap", status = "status",
                    cause = "cause", episode = NULL,
                    subject_weight = c("one", "followup"),
                    cause_weighting = c("shared", "literal_eq15"),
                    survival_side = c("right", "left"), tau = Inf) {
  cause_weighting <- match.arg(cause_weighting)
  survival_side <- match.arg(survival_side)

  d <- .validate_gap_data(data, id, gap, status, episode)
  if (!cause %in% names(d)) {
    stop("Missing cause column: ", cause)
  }

  levels <- unique(d[[cause]][d[[status]] == 1 & !is.na(d[[cause]])])
  if (length(levels) < 2) {
    stop("At least two recurrent event causes are required")
  }

  recs <- .subject_records(d, id, gap, status, cause = cause)
  times <- .event_grid(recs, tau)
  if (!length(times)) {
    stop("No completed recurrent gaps within tau")
  }

  if (is.function(subject_weight)) {
    a <- vapply(
      recs,
      function(r) {
        as.numeric(subject_weight(d[d[[id]] == r$id, , drop = FALSE]))
      },
      numeric(1)
    )
    sw_name <- "callable"
  } else {
    subject_weight <- match.arg(subject_weight)
    a <- if (subject_weight == "one") {
      rep(1, length(recs))
    } else {
      vapply(recs, function(r) r$followup, numeric(1))
    }
    sw_name <- subject_weight
  }

  if (any(!is.finite(a)) || any(a < 0)) {
    stop("subject weights must be finite and >= 0")
  }

  n <- length(recs)
  nt <- length(times)
  R_i <- matrix(0, n, nt)
  dG_i <- matrix(0, n, nt)
  dGl_i <- stats::setNames(
    lapply(levels, function(x) matrix(0, n, nt)),
    as.character(levels)
  )

  for (i in seq_len(n)) {
    r <- recs[[i]]
    base <- a[i] / r$mstar

    for (j in seq_along(times)) {
      tt <- times[j]
      R_i[i, j] <- base * sum(r$used_gaps >= tt)
      ev <- r$used_status == 1
      dG_i[i, j] <- base * sum(ev & r$used_gaps == tt)

      for (ll in seq_along(levels)) {
        lev <- levels[[ll]]
        mask <- ev & !is.na(r$used_causes) & r$used_causes == lev

        if (cause_weighting == "shared") {
          dGl_i[[as.character(lev)]][i, j] <-
            base * sum(mask & r$used_gaps == tt)
        } else {
          mil <- sum(mask)
          if (mil >= 2) {
            dGl_i[[as.character(lev)]][i, j] <-
              (a[i] / mil) * sum(mask & r$used_gaps == tt)
          }
        }
      }
    }
  }

  risk <- colSums(R_i)
  dG <- colSums(dG_i)
  dLambda <- ifelse(risk > 0, dG / risk, 0)
  Lambda <- cumsum(dLambda)
  S <- exp(-Lambda)
  S_left <- c(1, utils::head(S, -1))
  S_use <- if (survival_side == "right") S else S_left

  out <- data.frame(
    time = times,
    risk_mass = risk,
    dG_mass = dG,
    dLambda = dLambda,
    Lambda = Lambda,
    S_left = S_left,
    S = S,
    F_overall = 1 - S
  )

  cause_curves <- list()
  for (lev in levels) {
    key <- as.character(lev)
    dGl <- colSums(dGl_i[[key]])
    dL <- ifelse(risk > 0, dGl / risk, 0)
    dF <- S_use * dL
    F <- cumsum(dF)

    cause_curves[[key]] <- list(
      level = lev,
      dG_mass = dGl,
      dLambda = dL,
      dF = dF,
      F = F
    )

    safe <- gsub("[^A-Za-z0-9_.]", "_", key)
    out[[paste0("dG_", safe)]] <- dGl
    out[[paste0("dLambda_", safe)]] <- dL
    out[[paste0("dF_", safe)]] <- dF
    out[[paste0("F_", safe)]] <- F
  }

  sumdG <- Reduce(`+`, lapply(cause_curves, `[[`, "dG_mass"))
  add_err <- max(abs(sumdG - dG))

  list(
    method = "Sivadasan/Sankaran RCIF",
    cause_weighting = cause_weighting,
    survival_side = survival_side,
    subject_weight = sw_name,
    causes = levels,
    n_subjects = n,
    curve = out,
    cause_curves = cause_curves,
    additivity_max_abs_dG = add_err,
    warning = if (cause_weighting == "literal_eq15") {
      paste(
        "literal Eq.(15) cause weighting generally does not add to",
        "the pooled event process"
      )
    } else {
      NULL
    }
  )
}

#' Compute integrated RCIF equal-cause contrasts
#'
#' @param fit An object returned by [ss_rcif()].
#' @param weight Either `NULL` for unit weights, a function of `fit$curve`, or a
#'   numeric vector with one value per RCIF event time.
#'
#' @return A numeric vector with one integrated contrast per cause.
#' @keywords internal
#' @noRd
.integrated_rcif_contrast <- function(fit, weight = NULL) {
  curve <- fit$curve
  levels <- fit$causes
  k <- length(levels)

  if (is.null(weight)) {
    w <- rep(1, nrow(curve))
  } else if (is.function(weight)) {
    w <- as.numeric(weight(curve))
  } else {
    w <- as.numeric(weight)
  }

  if (length(w) != nrow(curve)) {
    stop("weight must have one value per RCIF event time")
  }

  dF0 <- diff(c(0, curve$F_overall))
  vapply(
    levels,
    function(lev) {
      dFl <- fit$cause_curves[[as.character(lev)]]$dF
      sum(w * (dFl - dF0 / k))
    },
    numeric(1)
  )
}
