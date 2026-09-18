# Luo-Huang weighted recurrent gap-time rank test

Implements the two-sample `G_rho*` test of Luo and Huang using the
weighted-risk-set representation of recurrent gap-time data. With
`rho = 0`, the statistic is the recurrent gap-time analogue of the
ordinary log-rank test. With `rho = 1`, it is an analogue of the
Peto-Prentice generalized Wilcoxon test.

## Usage

``` r
wc_logrank(
  data,
  group,
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  rho = 0,
  tau = Inf
)
```

## Arguments

- data:

  A data frame in long recurrent gap-time format.

- group:

  Character scalar naming a subject-level two-group variable.

- id:

  Character scalar naming the subject identifier column.

- gap:

  Character scalar naming the gap-duration column.

- status:

  Character scalar naming the event indicator column.

- episode:

  Optional character scalar naming the within-subject episode order
  column.

- rho:

  Nonnegative rank-weight parameter. `rho = 0` gives the recurrent
  log-rank analogue; `rho = 1` gives the recurrent Peto-Prentice
  analogue.

- tau:

  Optional upper truncation time.

## Value

A list containing the test statistic, standard error, `z` statistic,
chi-square statistic, p-value, subject-level influence contributions,
and event-time details.

## Details

The variance is estimated from subject-level martingale-style
contributions, retaining the within-subject dependence among recurrent
gaps.

The current implementation requires exactly two groups. Group ordering
is determined by first appearance among subjects; `group0` and `group1`
are returned explicitly so the score direction is interpretable.

## References

Luo X, Huang CY (2011). Analysis of recurrent gap time data using the
weighted risk-set method and the modified within-cluster resampling
method. *Statistics in Medicine*, 30, 301-311.
[doi:10.1002/sim.4074](https://doi.org/10.1002/sim.4074) .

## Examples

``` r
dat <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  arm = rep(c("A", "A", "B", "B"), each = 3)
)
wc_logrank(dat, group = "arm", episode = "episode", rho = 0)
#> $method
#> [1] "Luo-Huang WRS G_rho*"
#> 
#> $rho
#> [1] 0
#> 
#> $interpretation
#> [1] "recurrent gap-time log-rank analogue"
#> 
#> $group0
#> [1] "A"
#> 
#> $group1
#> [1] "B"
#> 
#> $n_subjects
#> [1] 4
#> 
#> $G
#> [1] 0.1041667
#> 
#> $se
#> [1] 0.05311479
#> 
#> $z
#> [1] 1.961161
#> 
#> $chisq
#> [1] 3.846154
#> 
#> $df
#> [1] 1
#> 
#> $p.value
#> [1] 0.0498602
#> 
#> $subject_influence
#>   id        psi event_score
#> 1  1 0.06944444  -0.4166667
#> 2  2 0.11111111  -0.2500000
#> 3  3 0.12500000   0.5000000
#> 4  4 0.11111111   0.5833333
#> 
#> $detail
#>   time risk risk_group1 expected_group_fraction  dN  dLambda0 surv_left weight
#> 1    1  4.0         2.0               0.5000000 1.0 0.2500000     1.000      1
#> 2    2  3.0         1.5               0.5000000 1.5 0.5000000     0.750      1
#> 3    3  1.5         0.5               0.3333333 1.0 0.6666667     0.375      1
#> 4    4  0.5         0.0               0.0000000 0.5 1.0000000     0.125      1
#> 
```
