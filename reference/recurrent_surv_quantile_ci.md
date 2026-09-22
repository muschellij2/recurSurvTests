# Confidence intervals for recurrent-survival quantiles

Estimates a quantile of a
[`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
Wang–Chang marginal gap-time survival curve or a
[`psh_surv()`](https://muschellij2.github.io/recurSurvTests/reference/psh_surv.md)
Pena–Strawderman–Hollander curve. For example, `survival_prob = 0.5`
estimates the median gap time, `inf {t: S(t) <= 0.5}`.

## Usage

``` r
recurrent_surv_quantile_ci(
  data,
  estimator = c("wc", "psh"),
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  tau = Inf,
  survival_prob = 0.5,
  conf.level = 0.95,
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

- survival_prob:

  Survival probability defining the quantile. `0.5` is the median;
  `0.25` is the 75th percentile of the gap-time distribution.

- conf.level:

  Confidence level in `(0, 1)`.

- B:

  Number of subject-bootstrap resamples. The default is 999; use 2,000
  or more for stable two-sided 95 percent tail limits when computation
  permits.

- seed:

  Optional random-number seed.

- keep_bootstrap:

  Logical; retain bootstrap quantiles when available.

## Value

A list containing the quantile `estimate`, its `lower` and `upper`
confidence limits, the fitted `curve`, method details, and bootstrap
diagnostics. An infinite estimate or upper limit means that the fitted
curve did not reach `survival_prob` before the available support.

## Details

The bootstrap resamples complete participant histories and refits the
curve. It therefore preserves recurrent-gap dependence and the
event-count weights of the WC estimator. Independence-style normal curve
limits are deliberately not inverted for a quantile interval because
they do not account for within-subject clustering.

A quantile interval, a pointwise curve interval, and a simultaneous
curve band have different coverage targets. Use
[`recurrent_surv_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_ci.md)
for the latter two curve-level quantities.

## Examples

``` r
# \donttest{
dat <- data.frame(
  id = rep(1:4, each = 3), episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4)
)
recurrent_surv_quantile_ci(dat, estimator = "wc", episode = "episode",
                          survival_prob = 0.5, B = 999, seed = 1)
#> $method
#> [1] "Wang-Chang subject-bootstrap percentile quantile interval"
#> 
#> $estimator
#> [1] "wc"
#> 
#> $survival_prob
#> [1] 0.5
#> 
#> $conf.level
#> [1] 0.95
#> 
#> $estimate
#> [1] 2
#> 
#> $lower
#> [1] 1
#> 
#> $upper
#> [1] 2
#> 
#> $fit
#>   time risk  dN   dLambda surv_left  surv
#> 1    1  4.0 1.0 0.2500000     1.000 0.750
#> 2    2  3.0 1.5 0.5000000     0.750 0.375
#> 3    3  1.5 1.0 0.6666667     0.375 0.125
#> 4    4  0.5 0.5 1.0000000     0.125 0.000
#> 
#> $bootstrap_success
#> [1] 999
#> 
#> $bootstrap_failures
#> [1] 0
#> 
#> $n_infinite
#> [1] 0
#> 
#> $prop_infinite
#> [1] 0
#> 
#> $warning
#> NULL
#> 
# }
```
