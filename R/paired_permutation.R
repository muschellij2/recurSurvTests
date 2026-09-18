#' Paired label-permutation test for Wang--Chang recurrent gap-time curves
#'
#' Tests a two-visit paired comparison by repeatedly swapping the two visit
#' labels within each pair and recalculating [wc_logrank()]. This supplies
#' randomization-based inference when within-person visits are paired; it is not
#' an independent-arm rank test.
#'
#' @param data Long recurrent gap-time data.
#' @param pair Character scalar naming the paired-participant column. Each pair
#'   must have exactly two distinct visit-level `id` values.
#' @param group Character scalar naming the two visit labels.
#' @param id,gap,status,episode,rho,tau Passed to [wc_logrank()].
#' @param B Number of random within-pair label permutations.
#' @param seed Optional random seed.
#'
#' @return A list with the observed Wang--Chang rank fit, permutation z
#' statistics, and a two-sided randomization p-value.
#'
#' @details
#' The procedure requires one observation from each of the two labels within
#' every pair. It tests the sharp within-pair label-exchangeability null. The
#' usual [wc_logrank()] p-value is also returned for comparison, but it does not
#' itself account for pairing of visits.
#'
#' @examples
#' dat <- data.frame(
#'   id = rep(1:4, each = 2), pair = rep(1:2, each = 4),
#'   episode = rep(1:2, 4), gap = c(2, 4, 3, 2, 1, 5, 2, 3),
#'   status = rep(c(1, 0), 4), visit = rep(c("A", "B"), each = 2)
#' )
#' wc_paired_permutation(dat, pair = "pair", group = "visit", episode = "episode",
#'   B = 99, seed = 1)
#'
#' @export
wc_paired_permutation <- function(data, pair, group, id = "id", gap = "gap",
                                  status = "status", episode = NULL, rho = 0,
                                  tau = Inf, B = 999, seed = NULL) {
  if (!pair %in% names(data)) stop("Missing pair column: ", pair)
  if (!group %in% names(data)) stop("Missing group column: ", group)
  if (length(B) != 1 || B < 1 || B != as.integer(B)) stop("B must be a positive integer")
  if (!is.null(seed)) set.seed(seed)
  d <- .validate_gap_data(data, id, gap, status, episode)
  visit <- unique(d[, c(id, pair, group), drop = FALSE])
  if (anyDuplicated(visit[[id]])) stop("pair and group must be constant within each id")
  labels <- unique(as.character(visit[[group]]))
  if (length(labels) != 2L) stop("wc_paired_permutation requires exactly two groups")
  by_pair <- split(seq_len(nrow(visit)), visit[[pair]], drop = TRUE)
  if (any(vapply(by_pair, length, integer(1)) != 2L) || any(vapply(by_pair, function(ii) {
    length(unique(as.character(visit[[group]][ii]))) != 2L
  }, logical(1)))) {
    stop("each pair must have exactly two ids, one in each group")
  }

  observed <- wc_logrank(d, group = group, id = id, gap = gap, status = status,
    episode = episode, rho = rho, tau = tau)
  perm_z <- numeric(B)
  for (b in seq_len(B)) {
    perm_group <- as.character(visit[[group]])
    for (ii in by_pair) if (stats::runif(1) < 0.5) perm_group[ii] <- rev(perm_group[ii])
    dd <- d
    dd$.paired_group_internal <- perm_group[match(dd[[id]], visit[[id]])]
    perm_z[b] <- wc_logrank(dd, group = ".paired_group_internal", id = id,
      gap = gap, status = status, episode = episode, rho = rho, tau = tau)$z
  }
  list(
    method = "paired within-pair label permutation for Luo-Huang WRS G_rho*",
    rho = rho, B = as.integer(B), observed = observed,
    permutation_z = perm_z,
    p.value = (1 + sum(abs(perm_z) >= abs(observed$z), na.rm = TRUE)) / (B + 1)
  )
}
