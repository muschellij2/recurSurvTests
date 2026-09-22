# Subject-bootstrap confidence intervals for recurrent survival curves

Forms confidence intervals for a
[`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
Wang–Chang marginal gap-time survival curve or a
[`psh_surv()`](https://muschellij2.github.io/recurSurvTests/reference/psh_surv.md)
Pena–Strawderman–Hollander (pooled Kaplan–Meier) curve. Entire subject
histories, rather than individual gap rows, are resampled. This
preserves the within-subject dependence and the event-count structure
that row-level resampling would destroy.

## Usage

``` r
recurrent_surv_ci(
  data,
  estimator = c("wc", "psh"),
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  tau = Inf,
  times = NULL,
  conf.level = 0.95,
  interval = c("pointwise", "simultaneous"),
  ci_method = c("bootstrap", "normal"),
  B = 999,
  seed = NULL,
  keep_bootstrap = FALSE
)
```

## Arguments

- data:

  A data frame in the long gap-time format used by
  [`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
  and
  [`psh_surv()`](https://muschellij2.github.io/recurSurvTests/reference/psh_surv.md).

- estimator:

  Either `"wc"` for Wang–Chang or `"psh"` for PSH/pooled KM.

- id, gap, status, episode, tau:

  Passed to the selected curve estimator.

- times:

  Times at which intervals are evaluated. The observed curve's event
  times are used by default. The fitted curve is evaluated as a
  right-continuous step function at supplied times.

- conf.level:

  Confidence level in `(0, 1)`.

- interval:

  Either `"pointwise"` (the standard choice for a survival curve) or
  `"simultaneous"`. Pointwise intervals are cluster-bootstrap percentile
  intervals and have the stated coverage separately at each fixed time.
  Simultaneous intervals are max-t bootstrap bands intended to cover all
  supplied times jointly.

- ci_method:

  Either `"bootstrap"` (default) or `"normal"`. The normal option forms
  pointwise intervals from the independence-style Nelson–Aalen/Greenwood
  standard-error approximation. It is included to reproduce the type of
  pointwise interval obtainable from legacy fitters; it is not generally
  appropriate for dependent recurrent gaps. Simultaneous bands require
  `ci_method = "bootstrap"`.

- B:

  Number of subject-level bootstrap resamples. At least 100 are required
  for `ci_method = "bootstrap"`; 999 or more are recommended for final
  analysis. The default of 999 follows common bootstrap and permutation
  practice. For stable two-sided 95 percent percentile limits or max-t
  band critical values, use 2,000 or more resamples when computation
  permits. It is ignored for `ci_method = "normal"`.

- seed:

  Optional random-number seed.

- keep_bootstrap:

  Logical; retain the bootstrap survival matrix in the returned object.
  This can be large when many time points are requested.

## Value

A list with the observed `fit`, a `curve` data frame containing `time`,
`surv`, `std.error`, `lower`, and `upper`, bootstrap diagnostics, and
(for simultaneous bands) the max-t critical value. `bootstrap_surv` is
included only when `keep_bootstrap = TRUE`.

## Details

Pointwise and simultaneous coverage answer different questions. A 95
percent pointwise interval has approximately 95 percent coverage at one
pre-specified gap time; it is not a 95 percent band for the whole curve.
Use `interval = "simultaneous"` when the inferential statement concerns
all supplied times at once. The max-t band uses the bootstrap standard
deviation at each time and is clipped to `[0, 1]`.

The 2.5th and 97.5th percentile limits of a two-sided 95 percent
interval are estimated from roughly 25 tail resamples at `B = 999` and
50 at `B = 2000`. Thus 999 is a useful default, while 2,000 or more
reduces Monte-Carlo variability in tail limits. A confidence interval
for a survival quantile or a between-group curve difference is a
different inferential target and is not returned by this function.

The resampling unit is the subject. This is important for both
estimators: PSH's pooled product limit has correlated gap rows, and
Wang–Chang's subject-level weights depend on the complete recurrence
history. Bootstrap intervals are not a substitute for checking
independent censoring and the estimand-specific assumptions of the
selected curve.

## References

Wang MC, Chang SH (1999). Nonparametric estimation of a recurrent
survival function. *Journal of the American Statistical Association*,
94, 146-153.

Pena EA, Strawderman RL, Hollander M (2001). Nonparametric estimation
with recurrent event data. *Journal of the American Statistical
Association*, 96, 1299-1315.
[doi:10.1198/016214501753381887](https://doi.org/10.1198/016214501753381887)
.

## Examples

``` r
# \donttest{
dat <- data.frame(
  id = rep(1:4, each = 3), episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4)
)
ci <- recurrent_surv_ci(dat, estimator = "wc", episode = "episode",
                        B = 999, seed = 1)
ci$curve
#>   time  surv std.error lower upper
#> 1    1 0.750 0.1254938 0.500 1.000
#> 2    2 0.375 0.1048812 0.125 0.500
#> 3    3 0.125 0.1071053 0.000 0.375
#> 4    4 0.000 0.0000000 0.000 0.000
# }
```
