# Sivadasan-Sankaran recurrent cumulative-incidence estimator

Implements a transparent reference version of the recurrent cumulative
incidence function (RCIF) estimator proposed by Sivadasan and Sankaran
for recurrent competing-risk gap-time data.

## Usage

``` r
ss_rcif(
  data,
  id = "id",
  gap = "gap",
  status = "status",
  cause = "cause",
  episode = NULL,
  subject_weight = c("one", "followup"),
  cause_weighting = c("shared", "literal_eq15"),
  survival_side = c("right", "left"),
  tau = Inf
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

  Character scalar naming the recurrent-event cause column. Causes
  should be non-missing on completed-event rows.

- episode:

  Optional character scalar naming the within-subject episode order
  column.

- subject_weight:

  Either `"one"`, giving each subject unit weight, `"followup"`,
  weighting subjects by total follow-up, or a function taking a
  subject-specific data frame and returning one nonnegative scalar.

- cause_weighting:

  Either `"shared"` (default), which uses the same subject-level `m_i*`
  denominator for pooled and cause-specific event processes, or
  `"literal_eq15"`, which reproduces the printed cause-specific
  denominator in Equation 15 of the 2023 paper.

- survival_side:

  Whether RCIF increments use the right-continuous (`"right"`) or
  left-limit (`"left"`) survival estimate at each event time.

- tau:

  Optional upper truncation time.

## Value

A list containing the pooled recurrent survival process, cause-specific
recurrent cumulative-incidence curves, subject-weighting metadata, and
an additivity diagnostic.

## Details

The printed Equation 15 in the 2023 *Statistica* paper uses a
cause-specific denominator `m_il` and an indicator requiring at least
two events of a cause. Taken literally, those cause-specific increments
need not sum to the pooled event process. Therefore the default is
`cause_weighting = "shared"`, while `"literal_eq15"` is provided for
literal reproduction and sensitivity analysis.

## References

Sivadasan SM, Sankaran PG (2023). Nonparametric estimation of cumulative
incidence functions of recurrent events. *Statistica*, 83(1), 3-25.

Sivadasan SM, Sankaran PG (2022). A nonparametric test for comparing
recurrent cumulative incidence functions. *METRON*.
[doi:10.1007/s40300-022-00228-x](https://doi.org/10.1007/s40300-022-00228-x)
.

## Examples

``` r
dat <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
)
ss_rcif(dat, episode = "episode")
#> $method
#> [1] "Sivadasan/Sankaran RCIF"
#> 
#> $cause_weighting
#> [1] "shared"
#> 
#> $survival_side
#> [1] "right"
#> 
#> $subject_weight
#> [1] "one"
#> 
#> $causes
#> [1] "A" "B"
#> 
#> $n_subjects
#> [1] 4
#> 
#> $curve
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
#> $cause_curves
#> $cause_curves$A
#> $cause_curves$A$level
#> [1] "A"
#> 
#> $cause_curves$A$dG_mass
#> [1] 0.5 1.0 0.0 0.5
#> 
#> $cause_curves$A$dLambda
#> [1] 0.1250000 0.3333333 0.0000000 1.0000000
#> 
#> $cause_curves$A$dF
#> [1] 0.09735010 0.15745552 0.00000000 0.08921852
#> 
#> $cause_curves$A$F
#> [1] 0.0973501 0.2548056 0.2548056 0.3440241
#> 
#> 
#> $cause_curves$B
#> $cause_curves$B$level
#> [1] "B"
#> 
#> $cause_curves$B$dG_mass
#> [1] 0.5 0.5 1.0 0.0
#> 
#> $cause_curves$B$dLambda
#> [1] 0.1250000 0.1666667 0.6666667 0.0000000
#> 
#> $cause_curves$B$dF
#> [1] 0.09735010 0.07872776 0.16168072 0.00000000
#> 
#> $cause_curves$B$F
#> [1] 0.0973501 0.1760779 0.3377586 0.3377586
#> 
#> 
#> 
#> $additivity_max_abs_dG
#> [1] 0
#> 
#> $warning
#> NULL
#> 
```
