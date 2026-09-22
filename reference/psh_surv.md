# Pena–Strawderman–Hollander recurrent survival estimator

Estimates the recurrent gap-time survival distribution using the
generalized product-limit estimator of Pena, Strawderman, and Hollander
(PSH). With the long gap-time data convention used in this package, PSH
is the ordinary Kaplan–Meier product limit formed from all observed gap
records. Thus each recurrent gap has equal weight, unlike
[`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md),
which weights each subject's eligible gaps by `1 / m_i*`.

## Usage

``` r
psh_surv(
  data,
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  tau = Inf
)
```

## Arguments

- data:

  A data frame in long gap-time format with one row per observed gap.

- id:

  Character scalar naming the subject identifier column.

- gap:

  Character scalar naming the gap-duration column.

- status:

  Character scalar naming the event indicator column.

- episode:

  Optional character scalar naming the within-subject episode order
  column.

- tau:

  Optional upper truncation time for the recurrent gap distribution.

## Value

A data frame with one row per distinct completed gap time and columns
`time`, `risk`, `dN`, `dLambda`, `surv_left`, and `surv`.

## Details

At a gap time `t`, PSH uses `R(t) = sum(i,j) I(Y_ij >= t)` and
`dN(t) = sum(i,j) I(Y_ij = t, delta_ij = 1)`, then applies
`prod(u <= t) {1 - dN(u) / R(u)}`. The renewal/IID interpretation of
this estimator differs from the marginal recurrent gap-time estimand of
Wang–Chang.

## References

Pena EA, Strawderman RL, Hollander M (2001). Nonparametric estimation
with recurrent event data. *Journal of the American Statistical
Association*, 96, 1299-1315.
[doi:10.1198/016214501753381887](https://doi.org/10.1198/016214501753381887)
.

## Examples

``` r
dat <- data.frame(
  id = rep(1:3, each = 3), episode = rep(1:3, 3),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5),
  status = rep(c(1, 1, 0), 3)
)
psh_surv(dat, episode = "episode")
#>   time risk dN   dLambda surv_left      surv
#> 1    1    9  1 0.1111111 1.0000000 0.8888889
#> 2    2    8  3 0.3750000 0.8888889 0.5555556
#> 3    3    5  1 0.2000000 0.5555556 0.4444444
#> 4    4    3  1 0.3333333 0.4444444 0.2962963
```
