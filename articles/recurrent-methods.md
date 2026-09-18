# Recurrent gap-time survival and competing-risk methods

``` r

library(recurSurvTests)
```

## Scope and data structure

Recurrent-event methods are often described with similar language while
targeting different quantities. This package is for **gap-time**
analyses: if subject $`i`$ has observed gap durations
$`Y_{i1}, \ldots, Y_{im_i}`$, each row of the input represents one such
gap. The final row is the right-censored gap (`status = 0`); preceding
rows are completed events (`status = 1`). The optional `episode` column
specifies their within-subject ordering.

This differs from a mean cumulative function, which accumulates events
on a calendar-time scale. The RCIF functions below also use recurrent
gap times, but target cause-specific incidence rather than a survival
curve. Thus a result from
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
should not be interpreted as a test of a calendar-time mean cumulative
function or an RCIF.

``` r

gaps <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  arm = rep(c("A", "A", "B", "B"), each = 3),
  cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
)
```

## Wang-Chang marginal gap-time survival

Let $`m_i^*`$ be the number of gap records used for subject $`i`$ in the
weighted-risk-set representation. At gap time $`t`$, the package
calculates

``` math
R(t) = \sum_i \frac{1}{m_i^*}\sum_j I(Y_{ij} \geq t), \qquad
dN(t) = \sum_i \frac{1}{m_i^*}\sum_j I(Y_{ij}=t,\ \delta_{ij}=1).
```

The estimated hazard increment is $`d\widehat\Lambda(t)=dN(t)/R(t)`$ and
the Wang-Chang-equivalent survival estimate is the product

``` math
\widehat S(t)=\prod_{u\leq t}\{1-d\widehat\Lambda(u)\}.
```

In code, `.subject_records()` validates and turns long rows into one
record per subject, `.event_grid()` finds distinct completed gap times,
and `.wrs_components()` constructs the subject-by-time risk and event
matrices.
[`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
sums those matrices and applies the product-limit update.

``` r

wc_surv(gaps, episode = "episode")
#>   time risk  dN   dLambda surv_left  surv
#> 1    1  4.0 1.0 0.2500000     1.000 0.750
#> 2    2  3.0 1.5 0.5000000     0.750 0.375
#> 3    3  1.5 1.0 0.6666667     0.375 0.125
#> 4    4  0.5 0.5 1.0000000     0.125 0.000
```

## Luo-Huang recurrent rank test

For two subject-level groups,
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
applies the Luo-Huang weighted risk-set $`G_\rho^*`$ statistic. At each
event time it contrasts the observed weighted events in group 1 with its
weighted risk-set expectation, using $`w(t)=\widehat S(t-)^\rho`$.
`rho = 0` is the recurrent gap-time log-rank analogue and `rho = 1` is
the Peto-Prentice analogue. The variance is formed from subject-level
martingale-style contributions, retaining within-subject dependence.

``` r

rank_fit <- wc_logrank(gaps, group = "arm", episode = "episode", rho = 0)
rank_fit[c("method", "interpretation", "G", "se", "p.value")]
#> $method
#> [1] "Luo-Huang WRS G_rho*"
#> 
#> $interpretation
#> [1] "recurrent gap-time log-rank analogue"
#> 
#> $G
#> [1] 0.1041667
#> 
#> $se
#> [1] 0.05311479
#> 
#> $p.value
#> [1] 0.0498602
```

This is Luo-Huang’s implementation, not an exact Zhao et al. (2020)
implementation; the available Zhao supplement does not supply enough
formulas or source code to reconstruct the latter.

## Zhao extended rank tests and simulation checks

[`zhao_rank_test()`](https://muschellij2.github.io/recurSurvTests/reference/zhao_rank_test.md)
is separate from
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
and implements the extended LR, GB, and PP score construction in Zhao et
al. (2020). Its default `variance_method = "pooled_risk"` has subject
residuals that sum exactly to the score and is covered by the package’s
null simulation tests. `"zhao_eq6"` follows the group-specific
denominator printed in Equation 6 as a sensitivity/reproduction
implementation. This distinction is explicit because the paper’s
published numerical tables cannot be regenerated exactly without its
original simulation software.

``` r

zhao_rank_test(gaps, group = "arm", episode = "episode",
               variance_method = "pooled_risk")
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
```

The PSG scenario treats a recording as 8 hours (480 minutes), with no
staggered entry, and scales a chosen gap distribution to target roughly
200 completed events per participant. Its default `Z ~ U(0.75, 1.25)`
avoids the extreme event-count dispersion in the paper’s `U(0.1, 1.9)`
stress scenario.

``` r

zhao_rank_simulation(n_sim = 1000, n_per_group = 100, scenario = "psg_200",
                     variance_method = "pooled_risk", seed = 1)
```

## Recurrent competing risks

[`ss_rcif()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif.md)
computes a pooled recurrent hazard and cause-specific increments. With
subject weight $`a_i`$, the default shared weighting uses $`a_i/m_i^*`$
for both pooled and cause-specific event processes. For cause $`l`$,

``` math
d\widehat F_l(t)=\widehat S(t)\,
\frac{dG_l(t)}{R(t)}, \qquad
\widehat F_l(t)=\sum_{u\leq t}d\widehat F_l(u).
```

The default `cause_weighting = "shared"` maintains pooled/cause
increment additivity. `"literal_eq15"` reproduces the printed Equation
15 denominator as a sensitivity option; its increments can be
non-additive and should be checked through `additivity_max_abs_dG`.

``` r

rcif <- ss_rcif(gaps, episode = "episode")
rcif$curve
#>   time risk_mass dG_mass   dLambda   Lambda    S_left          S F_overall dG_A
#> 1    1       4.0     1.0 0.2500000 0.250000 1.0000000 0.77880078 0.2211992  0.5
#> 2    2       3.0     1.5 0.5000000 0.750000 0.7788008 0.47236655 0.5276334  1.0
#> 3    3       1.5     1.0 0.6666667 1.416667 0.4723666 0.24252107 0.7574789  0.0
#> 4    4       0.5     0.5 1.0000000 2.416667 0.2425211 0.08921852 0.9107815  0.5
#>   dLambda_A       dF_A       F_A dG_B dLambda_B       dF_B       F_B
#> 1 0.1250000 0.09735010 0.0973501  0.5 0.1250000 0.09735010 0.0973501
#> 2 0.3333333 0.15745552 0.2548056  0.5 0.1666667 0.07872776 0.1760779
#> 3 0.0000000 0.00000000 0.2548056  1.0 0.6666667 0.16168072 0.3377586
#> 4 1.0000000 0.08921852 0.3440241  0.0 0.0000000 0.00000000 0.3377586
rcif$additivity_max_abs_dG
#> [1] 0
```

The optional
[`ss_rcif_equal_causes_test()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif_equal_causes_test.md)
integrates contrasts of $`F_l(t)-F(t)/k`$, uses Helmert contrasts, and
estimates covariance by resampling subjects. It is a bootstrap Wald
reconstruction of the METRON equal-cause RCIF test, not a log-rank test.
For real analyses use a substantially larger `B` (for example, 999 or
more) than would be appropriate for a quick example.
