# Zhao et al. extended recurrent gap-time rank test

Tests equality of two Wang–Chang marginal recurrent gap-time survival
distributions using the extended rank-test construction of Zhao et al.
(2020). Unlike
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md),
this function uses Zhao et al.'s adjusted subject-level residuals (their
Equation 6) to estimate the robust variance.

## Usage

``` r
zhao_rank_test(
  data,
  group,
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  test = c("logrank", "gehan_breslow", "peto_prentice"),
  variance_method = c("pooled_risk", "zhao_eq6"),
  alternative = c("two.sided", "greater", "less"),
  tau = Inf
)
```

## Arguments

- data:

  A data frame in long recurrent gap-time format.

- group:

  Character scalar naming a subject-level two-group variable.

- id, gap, status, episode:

  See
  [`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md).

- test:

  Rank weight: `"logrank"` (LR), `"gehan_breslow"` (GB), or
  `"peto_prentice"` (PP).

- variance_method:

  Robust variance residual. The default `"pooled_risk"` uses a pooled
  risk-set denominator and has residuals that sum exactly to the score.
  `"zhao_eq6"` follows the group-specific denominator printed in Zhao et
  al.'s Equation 6 as a sensitivity/reproduction implementation.

- alternative:

  Alternative for the normal-score p-value.

- tau:

  Optional upper truncation time.

## Value

A list with the unstandardized score `W`, robust variance `variance`,
standard error, normal score `z`, p-value, subject residuals, and
event-time details. `group1` is the treatment group, determined by
second appearance among subjects; a positive score corresponds to more
events in `group1`.

## Details

The implementation follows Equations 3, 6, and 7 of Zhao et al. (2020).
The paper presents its calculation without tied observed times. Here,
tied events are handled by aggregating their Wang–Chang weighted
increments at a common time; this is a natural extension, but is not
separately derived in that paper. Peto–Prentice weights use the pooled
Wang–Chang survival estimate `exp(-cumulative hazard)`, as defined in
the paper.

## References

Zhao Q, Zhang B, LaValley MP, Massaro JM, Lunetta KL, Chang M (2020).
Extended rank tests for analyzing recurrent event data. *Statistics in
Biopharmaceutical Research*, 12(1), 90-98.
[doi:10.1080/19466315.2019.1601596](https://doi.org/10.1080/19466315.2019.1601596)
.

## Examples

``` r
dat <- data.frame(
  id = rep(1:4, each = 3), episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  arm = rep(c("A", "A", "B", "B"), each = 3)
)
zhao_rank_test(dat, group = "arm", episode = "episode")
#> $method
#> [1] "Zhao et al. extended logrank rank test"
#> 
#> $test
#> [1] "logrank"
#> 
#> $variance_method
#> [1] "pooled_risk"
#> 
#> $alternative
#> [1] "two.sided"
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
#> $W
#> [1] 0.4166667
#> 
#> $variance
#> [1] 0.04513889
#> 
#> $se
#> [1] 0.2124591
#> 
#> $z
#> [1] 1.961161
#> 
#> $p.value
#> [1] 0.0498602
#> 
#> $subject_residuals
#>   id group score_contribution adjustment   residual
#> 1  1     A         -0.4166667 -0.4861111 0.06944444
#> 2  2     A         -0.2500000 -0.3611111 0.11111111
#> 3  3     B          0.5000000  0.3750000 0.12500000
#> 4  4     B          0.5833333  0.4722222 0.11111111
#> 
#> $detail
#>   time risk risk_group1  dN   dLambda   survival weight
#> 1    1  4.0         2.0 1.0 0.2500000 0.77880078      1
#> 2    2  3.0         1.5 1.5 0.5000000 0.47236655      1
#> 3    3  1.5         0.5 1.0 0.6666667 0.24252107      1
#> 4    4  0.5         0.0 0.5 1.0000000 0.08921852      1
#> 
```
