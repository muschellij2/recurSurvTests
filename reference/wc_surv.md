# Wang-Chang weighted-risk-set recurrent gap-time survival estimator

Estimates the marginal survival distribution of recurrent gap times
using the subject-weighted risk-set representation that is algebraically
equivalent to the Wang-Chang estimator and used by Luo and Huang.

## Usage

``` r
wc_surv(
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

  Character scalar naming the event indicator column. Completed
  recurrent gaps are coded `1`; the final censored gap is coded `0`.

- episode:

  Optional character scalar naming the within-subject episode order
  column.

- tau:

  Optional upper truncation time for the recurrent gap distribution.

## Value

A data frame with one row per distinct completed gap time and columns
`time`, `risk`, `dN`, `dLambda`, `surv_left`, and `surv`.

## Details

Each subject contributes weight `1 / m_i*` to each eligible recurrent
gap, where `m_i*` is the number of completed gaps used by the
weighted-risk-set construction. This avoids allowing subjects with many
recurrent events to dominate the marginal gap-time survival estimate.

## References

Wang MC, Chang SH (1999). Nonparametric estimation of a recurrent
survival function. *Journal of the American Statistical Association*,
94, 146-153.

Luo X, Huang CY (2011). Analysis of recurrent gap time data using the
weighted risk-set method and the modified within-cluster resampling
method. *Statistics in Medicine*, 30, 301-311.
[doi:10.1002/sim.4074](https://doi.org/10.1002/sim.4074) .

## Examples

``` r
dat <- data.frame(
  id = rep(1:3, c(3, 2, 3)),
  gap = c(2, 3, 4, 1, 5, 2, 4, 3),
  status = c(1, 1, 0, 1, 0, 1, 1, 0),
  episode = c(1, 2, 3, 1, 2, 1, 2, 3)
)
wc_surv(dat, episode = "episode")
#>   time risk  dN   dLambda surv_left      surv
#> 1    1  3.0 1.0 0.3333333 1.0000000 0.6666667
#> 2    2  2.0 1.0 0.5000000 0.6666667 0.3333333
#> 3    3  1.0 0.5 0.5000000 0.3333333 0.1666667
#> 4    4  0.5 0.5 1.0000000 0.1666667 0.0000000
```
