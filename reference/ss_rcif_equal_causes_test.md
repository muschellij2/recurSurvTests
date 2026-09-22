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
# \donttest{
dat <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
)
ss_rcif_equal_causes_test(dat, episode = "episode", B = 999, seed = 1)
#> $method
#> [1] "Sivadasan/Sankaran METRON RCIF equal-causes test (bootstrap reconstruction)"
#> 
#> $null
#> [1] "F_1(t)=...=F_k(t)=F(t)/k"
#> 
#> $logrank_analogue
#> [1] FALSE
#> 
#> $v
#>          A          B 
#> -0.1113666 -0.1176322 
#> 
#> $contrast
#> [1] 0.00443042
#> 
#> $chisq
#> [1] 0.0009910507
#> 
#> $df
#> [1] 1
#> 
#> $p.value.chisq
#> [1] 0.974886
#> 
#> $p.value.bootstrap.wald
#> [1] 0.8709677
#> 
#> $bootstrap_success
#> [1] 991
#> 
#> $bootstrap_failures
#> [1] 8
#> 
#> $cov_v
#>              [,1]         [,2]
#> [1,]  0.010521834 -0.009751331
#> [2,] -0.009751331  0.009587243
#> 
#> $fit
#> $fit$method
#> [1] "Sivadasan/Sankaran RCIF"
#> 
#> $fit$cause_weighting
#> [1] "shared"
#> 
#> $fit$survival_side
#> [1] "right"
#> 
#> $fit$subject_weight
#> [1] "one"
#> 
#> $fit$causes
#> [1] "A" "B"
#> 
#> $fit$n_subjects
#> [1] 4
#> 
#> $fit$curve
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
#> 
#> $fit$cause_curves
#> $fit$cause_curves$A
#> $fit$cause_curves$A$level
#> [1] "A"
#> 
#> $fit$cause_curves$A$dG_mass
#> [1] 0.5 1.0 0.0 0.5
#> 
#> $fit$cause_curves$A$dLambda
#> [1] 0.1250000 0.3333333 0.0000000 1.0000000
#> 
#> $fit$cause_curves$A$dF
#> [1] 0.09735010 0.15745552 0.00000000 0.08921852
#> 
#> $fit$cause_curves$A$F
#> [1] 0.0973501 0.2548056 0.2548056 0.3440241
#> 
#> 
#> $fit$cause_curves$B
#> $fit$cause_curves$B$level
#> [1] "B"
#> 
#> $fit$cause_curves$B$dG_mass
#> [1] 0.5 0.5 1.0 0.0
#> 
#> $fit$cause_curves$B$dLambda
#> [1] 0.1250000 0.1666667 0.6666667 0.0000000
#> 
#> $fit$cause_curves$B$dF
#> [1] 0.09735010 0.07872776 0.16168072 0.00000000
#> 
#> $fit$cause_curves$B$F
#> [1] 0.0973501 0.1760779 0.3377586 0.3377586
#> 
#> 
#> 
#> $fit$additivity_max_abs_dG
#> [1] 0
#> 
#> $fit$warning
#> NULL
#> 
#> 
# }
```
