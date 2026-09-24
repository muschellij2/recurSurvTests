#' Validate recurrent gap-time input data
#'
#' Internal validator for the long-form gap-time data convention used throughout
#' the package.
#'
#' @param data A data frame with one row per observed gap.
#' @param id Character scalar naming the subject identifier column.
#' @param gap Character scalar naming the gap-duration column.
#' @param status Character scalar naming the event indicator column. Completed
#'   recurrent gaps are coded `1`; the terminal censored gap is coded `0`.
#' @param episode Optional character scalar naming a within-subject episode-order
#'   column.
#'
#' @return A reordered validated data frame.
#' @keywords internal
#' @noRd
.validate_gap_data <- function(data, id, gap, status, episode = NULL) {
  req <- c(id, gap, status)
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }

  d <- as.data.frame(data)
  d$.row_order_internal <- seq_len(nrow(d))
  if (!is.null(episode)) {
    if (!episode %in% names(d)) {
      stop("Missing episode column: ", episode)
    }
    d <- d[order(d[[id]], d[[episode]], d$.row_order_internal), , drop = FALSE]
  } else {
    d <- d[order(d[[id]], d$.row_order_internal), , drop = FALSE]
  }
  d$.row_order_internal <- NULL

  if (anyNA(d[[gap]]) || any(d[[gap]] < 0)) {
    stop("gap must be non-missing and >= 0")
  }
  if (!all(d[[status]] %in% c(0, 1))) {
    stop("status must be 0/1")
  }

  spl <- split(seq_len(nrow(d)), d[[id]], drop = TRUE)
  for (ii in spl) {
    s <- d[[status]][ii]
    sid <- d[[id]][ii[1]]
    if (utils::tail(s, 1) != 0) {
      stop("Subject ", sid, ": final row must have status=0")
    }
    if (length(s) > 1 && any(utils::head(s, -1) != 1)) {
      stop(
        "Subject ", sid,
        ": all rows before final censored gap must have status=1"
      )
    }
  }

  d
}

#' Convert long-form gap data to subject records
#'
#' Internal helper that creates one list element per subject and computes the
#' Wang-Chang `m_i*` quantity used by the reference implementations.
#'
#' @param d A validated data frame returned by `.validate_gap_data()`.
#' @param id,gap,status Column names.
#' @param group Optional grouping column name.
#' @param cause Optional recurrent-event cause column name.
#'
#' @return A list of subject-level records.
#' @keywords internal
#' @noRd
.subject_records <- function(d, id, gap, status, group = NULL, cause = NULL) {
  # `d` is already ordered by id.  Splitting row indices once avoids scanning
  # the full data frame again for every subject, which is quadratic in the
  # number of subjects for large Monte-Carlo reference populations.
  row_index <- split(seq_len(nrow(d)), d[[id]], drop = TRUE)

  lapply(row_index, function(ii) {
    g <- d[ii, , drop = FALSE]
    sid <- d[[id]][ii[1]]
    gaps <- as.numeric(g[[gap]])
    stat <- as.integer(g[[status]])
    m <- length(gaps)
    mstar <- if (m == 1) 1L else m - 1L

    z <- list(
      id = sid,
      m = m,
      mstar = mstar,
      gaps = gaps,
      status = stat,
      used_gaps = gaps[seq_len(mstar)],
      used_status = stat[seq_len(mstar)],
      followup = sum(gaps)
    )

    if (!is.null(group)) {
      vals <- unique(g[[group]])
      if (length(vals) != 1) {
        stop("group must be constant within subject ", sid)
      }
      z$group <- vals[[1]]
    }

    if (!is.null(cause)) {
      ca <- g[[cause]]
      z$causes <- ca
      z$used_causes <- ca[seq_len(mstar)]
    }

    z
  })
}

#' Construct the recurrent gap event-time grid
#'
#' @param recs Subject-level records.
#' @param tau Upper truncation time.
#'
#' @return Sorted unique completed recurrent gap times no greater than `tau`.
#' @keywords internal
#' @noRd
.event_grid <- function(recs, tau = Inf) {
  z <- unlist(
    lapply(recs, function(r) r$used_gaps[r$used_status == 1]),
    use.names = FALSE
  )
  z <- sort(unique(as.numeric(z)))
  z[z <= tau]
}

#' Compute weighted-risk-set components
#'
#' Computes subject-level weighted risk-set and event increments under the
#' Wang-Chang / Luo-Huang weighted-risk-set representation.
#'
#' @param recs Subject-level records.
#' @param times Event-time grid.
#'
#' @return A list with matrices `R_i` and `dN_i`.
#' @keywords internal
#' @noRd
.wrs_components <- function(recs, times) {
  n <- length(recs)
  nt <- length(times)
  R_i <- matrix(0, n, nt)
  dN_i <- matrix(0, n, nt)

  if (!nt) {
    return(list(R_i = R_i, dN_i = dN_i))
  }

  for (i in seq_len(n)) {
    r <- recs[[i]]
    ww <- 1 / r$mstar
    for (j in seq_along(times)) {
      tt <- times[j]
      R_i[i, j] <- ww * sum(r$used_gaps >= tt)
      dN_i[i, j] <- ww * sum(r$used_status == 1 & r$used_gaps == tt)
    }
  }

  list(R_i = R_i, dN_i = dN_i)
}

#' Helmert-style contrast rows
#'
#' @param k Number of causes.
#'
#' @return A `(k - 1) x k` orthonormalized contrast matrix.
#' @keywords internal
#' @noRd
.helmert_rows <- function(k) {
  C <- matrix(0, k - 1, k)
  for (j in 2:k) {
    C[j - 1, 1:(j - 1)] <- 1 / (j - 1)
    C[j - 1, j] <- -1
  }
  C / sqrt(rowSums(C^2))
}

#' Symmetric generalized inverse and numerical rank
#'
#' @param A Symmetric matrix.
#' @param rtol Relative eigenvalue tolerance.
#'
#' @return A list with `inv`, the generalized inverse, and `rank`.
#' @keywords internal
#' @noRd
.pinv_rank <- function(A, rtol = 1e-10) {
  ee <- eigen((A + t(A)) / 2, symmetric = TRUE)
  vmax <- max(abs(ee$values), 1)
  keep <- ee$values > rtol * vmax
  rk <- sum(keep)

  if (!rk) {
    return(list(inv = matrix(0, nrow(A), ncol(A)), rank = 0L))
  }

  V <- ee$vectors[, keep, drop = FALSE]
  list(
    inv = V %*% diag(1 / ee$values[keep], nrow = rk) %*% t(V),
    rank = as.integer(rk)
  )
}

#' Subject-level bootstrap resample
#'
#' Resamples subjects with replacement while assigning new bootstrap subject
#' identifiers.
#'
#' @param d Long-form data frame.
#' @param id Subject identifier column name.
#'
#' @return A bootstrap-resampled long-form data frame.
#' @keywords internal
#' @noRd
.cluster_bootstrap_df <- function(d, id) {
  ids <- unique(d[[id]])
  draw <- sample(ids, length(ids), replace = TRUE)
  pieces <- lapply(seq_along(draw), function(ii) {
    z <- d[d[[id]] == draw[ii], , drop = FALSE]
    z[[id]] <- ii
    z
  })
  do.call(rbind, pieces)
}
