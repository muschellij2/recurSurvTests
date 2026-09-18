# Bootstrap equal-cause test for recurrent cumulative-incidence functions

Constructs integrated cause-specific RCIF contrasts and estimates their
covariance by resampling subjects with replacement. The resulting
quadratic form is a transparent Wald reconstruction of the equal-cause
test described by Sivadasan and Sankaran.

## Usage

``` r
ss_rcif_equal_causes_test(
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
  seed = NULL
)
```

## Arguments

- data:

  A data frame in long recurrent gap-time format.

- id:

  Character scalar naming the subject identifier column.

- gap:

  Character scalar naming the gap-duration column.

- status:

  Character scalar naming the event indicator column.

- cause:

  Character scalar naming the recurrent-event cause column.

- episode:

  Optional character scalar naming the within-subject episode order
  column.

- subject_weight:

  Passed to
  [`ss_rcif()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif.md).

- cause_weighting:

  Passed to
  [`ss_rcif()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif.md).

- survival_side:

  Passed to
  [`ss_rcif()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif.md).

- tau:

  Optional upper truncation time.

- weight:

  Either `NULL` for unit weights, a function of the fitted RCIF curve
  returning one weight per event time, or a numeric vector with one
  value per event time. A function is recommended for bootstrap analyses
  because bootstrap event-time grids can vary.

- B:

  Number of subject-level bootstrap replicates. At least 20 are
  required; 999 or more are recommended for analysis.

- seed:

  Optional random-number seed.

## Value

A list containing the observed integrated contrasts, Helmert contrast
vector, Wald chi-square statistic, effective degrees of freedom,
asymptotic and bootstrap-calibrated p-values, bootstrap diagnostics,
covariance matrix, and the fitted RCIF object.

## Details

The null hypothesis is that the cause-specific recurrent
cumulative-incidence functions are equal, equivalently
`F_l(t) = F(t) / k` for every cause under the formulation used here.

This implementation should be treated as a **reconstruction**, not as an
exact reproduction of unpublished author software. The published
appendix motivates integrated contrasts of `F_l - F/k`, while the
covariance here is estimated transparently by subject-level bootstrap.

## References

Sivadasan SM, Sankaran PG (2022). A nonparametric test for comparing
recurrent cumulative incidence functions. *METRON*.
[doi:10.1007/s40300-022-00228-x](https://doi.org/10.1007/s40300-022-00228-x)
.

Sivadasan SM, Sankaran PG (2023). Nonparametric estimation of cumulative
incidence functions of recurrent events. *Statistica*, 83(1), 3-25.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
)
ss_rcif_equal_causes_test(dat, episode = "episode", B = 999, seed = 1)
} # }
```
