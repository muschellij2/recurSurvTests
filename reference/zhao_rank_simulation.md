# Simulate Zhao et al. null-size experiments

Generates recurrent gap-time data under the simulation design in Section
3 of Zhao et al. (2020) and repeatedly applies
[`zhao_rank_test()`](https://muschellij2.github.io/recurSurvTests/reference/zhao_rank_test.md).
This is a compact, reproducible diagnostic for nominal type-I-error
calibration; it is not intended to reproduce the paper's
100,000-replicate tables during a routine package test.

## Usage

``` r
zhao_rank_simulation(
  n_sim = 1000,
  n_per_group = 100,
  scenario = c("paper", "psg_200"),
  distribution = c("exponential", "weibull", "lognormal", "loglogistic"),
  heterogeneity = NULL,
  followup = NULL,
  target_events = 200,
  time_resolution = NULL,
  alpha = 0.05,
  test = c("logrank", "gehan_breslow", "peto_prentice"),
  variance_method = c("pooled_risk", "zhao_eq6"),
  alternative = c("two.sided", "greater", "less"),
  treatment_time_multiplier = 1,
  seed = NULL
)
```

## Arguments

- n_sim:

  Number of simulated datasets.

- n_per_group:

  Subjects per group.

- scenario:

  `"paper"` reproduces the broad design choices in Zhao et al.'s
  Section 3. `"psg_200"` represents a fixed 8-hour PSG recording in
  minutes, calibrated for approximately `target_events` completed events
  per participant.

- distribution:

  Baseline gap-time distribution used in Zhao et al.'s Table 1.

- heterogeneity:

  Two endpoints of the subject multiplier `Z ~ U(a, b)`.

- followup:

  Administrative study end time.

- target_events:

  Target mean number of completed events per participant for
  `scenario = "psg_200"`.

- time_resolution:

  Optional gap-time discretization. The PSG default is 0.5 minutes
  (30-second scoring epochs); leave `NULL` for continuous times.

- alpha:

  Nominal test level.

- test:

  Passed to
  [`zhao_rank_test()`](https://muschellij2.github.io/recurSurvTests/reference/zhao_rank_test.md).

- variance_method:

  Passed to
  [`zhao_rank_test()`](https://muschellij2.github.io/recurSurvTests/reference/zhao_rank_test.md).

- alternative:

  Passed to
  [`zhao_rank_test()`](https://muschellij2.github.io/recurSurvTests/reference/zhao_rank_test.md).
  Zhao et al.'s Tables 2–3 use the one-sided `"less"` alternative when
  `group = "treatment"` has longer gaps (a lower recurrent-event
  hazard).

- treatment_time_multiplier:

  Multiplicative factor applied to treatment gap times. The default `1`
  generates the null; values above one generate longer treatment gaps
  for power simulations.

- seed:

  Optional random seed.

## Value

A list containing the empirical rejection rate, its Monte-Carlo standard
error, normal-score moments, p-values, and simulation settings.

## References

Zhao et al. (2020), Section 3 and Tables 1–3.

## Examples

``` r
zhao_rank_simulation(n_sim = 100, n_per_group = 50, seed = 1)
#> $method
#> [1] "Zhao et al. recurrent gap-time null simulation"
#> 
#> $empirical_type1_error
#> [1] 0.03
#> 
#> $monte_carlo_se
#> [1] 0.01705872
#> 
#> $mean_z
#> [1] -0.01045476
#> 
#> $sd_z
#> [1] 0.9819088
#> 
#> $mean_completed_events
#> [1] 2.6862
#> 
#> $completed_events
#>   [1] 2.34 2.44 2.53 2.48 2.75 2.53 2.05 2.49 2.43 2.85 2.56 2.73 2.88 2.65 3.23
#>  [16] 3.16 2.08 2.74 2.74 2.18 2.74 2.65 2.46 2.71 2.53 2.38 2.67 3.37 2.48 2.36
#>  [31] 2.43 2.57 3.02 2.78 2.18 2.39 2.62 2.35 2.48 2.98 2.49 2.96 2.44 2.64 3.28
#>  [46] 2.90 2.21 3.63 2.32 2.68 2.29 2.97 3.24 3.15 2.65 2.84 2.75 2.95 2.49 3.01
#>  [61] 2.09 2.51 2.87 3.31 2.21 2.70 2.64 2.55 2.82 2.07 2.86 3.30 2.50 2.00 2.38
#>  [76] 2.72 2.29 2.48 2.77 3.05 2.94 2.65 3.06 3.32 2.42 2.54 3.10 2.80 2.63 2.95
#>  [91] 2.98 2.35 3.37 2.55 3.64 2.33 2.57 3.24 2.36 2.82
#> 
#> $p.value
#>   [1] 0.42759826 0.25455788 0.88657363 0.72044762 0.76990219 0.42291364
#>   [7] 0.03564831 0.30734748 0.56168774 0.43743168 0.46007301 0.91368288
#>  [13] 0.13255998 0.38680933 0.06117539 0.48748937 0.44475941 0.09147498
#>  [19] 0.27840008 0.31896408 0.31770888 0.17215255 0.65637182 0.10167055
#>  [25] 0.78239610 0.76576748 0.92122523 0.39338447 0.28304607 0.40954796
#>  [31] 0.97933062 0.52854265 0.83516933 0.49150491 0.70495140 0.32600236
#>  [37] 0.02671187 0.05079297 0.23788046 0.84307990 0.75755375 0.59470986
#>  [43] 0.38702529 0.96826004 0.64943607 0.48295417 0.81848982 0.23198970
#>  [49] 0.68445096 0.79619158 0.29856677 0.73565656 0.15222865 0.32623095
#>  [55] 0.67324821 0.64633998 0.31758363 0.60934800 0.44093859 0.32370875
#>  [61] 0.76768741 0.14817340 0.22457104 0.78048597 0.01045672 0.37225657
#>  [67] 0.95581496 0.63129835 0.61502456 0.35783904 0.10302232 0.91703033
#>  [73] 0.55684891 0.12026785 0.17976007 0.06142051 0.63670218 0.05441719
#>  [79] 0.99178546 0.57249701 0.49240766 0.77996101 0.22442126 0.22192111
#>  [85] 0.49250383 0.99442298 0.21816216 0.77187737 0.26988574 0.74088635
#>  [91] 0.84649554 0.93331085 0.14207898 0.47694984 0.48445774 0.48256955
#>  [97] 0.91763827 0.33514423 0.50146228 0.14627453
#> 
#> $z
#>   [1]  0.793308229 -1.139348372 -0.142641112 -0.357860630 -0.292502844
#>   [6] -0.801376623  2.100916318  1.020803631  0.580336473  0.776537519
#>  [11] -0.738726631 -0.108394349 -1.504083073 -0.865418144 -1.872225677
#>  [16] -0.694307537  0.764180958 -1.687667992  1.083920429  0.996588932
#>  [21]  0.999177152 -1.365319818  0.444927956 -1.636808217  0.276197863
#>  [26] -0.297915738 -0.098890477  0.853495708  1.073501704 -0.824689381
#>  [31]  0.025908121 -0.630232233 -0.208076380 -0.687917228  0.378645130
#>  [36]  0.982197905  2.215702712  1.953224770  1.180301154 -0.197955471
#>  [41]  0.308694740  0.532023303 -0.865024587 -0.039790637 -0.454545750
#>  [46]  0.701559098 -0.229487777  1.195249150 -0.406396952  0.258279011
#>  [51]  1.039511817 -0.337610728 -1.431703670 -0.981733913  0.421694235
#>  [56] -0.458852658  0.999435799  0.511004362 -0.770609235  0.986865142
#>  [61]  0.295401214 -1.446013568 -1.214462806  0.278685787  2.560348773
#>  [66] -0.892254441 -0.055406069 -0.479900208 -0.502914260 -0.919490566
#>  [71]  1.630374719 -0.104175178 -0.587528358  1.553650305  1.341494129
#>  [76]  1.870456048 -0.472314695 -1.923500797  0.010295578 -0.564377851
#>  [81] -0.686484468  0.279369844 -1.214855375 -1.221435665 -0.686331917
#>  [86] -0.006989815 -1.231429795 -0.289920096  1.103325723  0.330679818
#>  [91] -0.193591774  0.083680011  1.468092920  0.711216888  0.699150885
#>  [96] -0.702175780 -0.103409122 -0.963803946 -0.672190736 -1.452817164
#> 
#> $settings
#> $settings$n_sim
#> [1] 100
#> 
#> $settings$n_per_group
#> [1] 50
#> 
#> $settings$scenario
#> [1] "paper"
#> 
#> $settings$distribution
#> [1] "exponential"
#> 
#> $settings$heterogeneity
#> [1] 0.1 1.9
#> 
#> $settings$followup
#> [1] 180
#> 
#> $settings$target_events
#> NULL
#> 
#> $settings$mean_gap
#> NULL
#> 
#> $settings$time_resolution
#> NULL
#> 
#> $settings$alpha
#> [1] 0.05
#> 
#> $settings$test
#> [1] "logrank"
#> 
#> $settings$variance_method
#> [1] "pooled_risk"
#> 
#> $settings$alternative
#> [1] "two.sided"
#> 
#> $settings$treatment_time_multiplier
#> [1] 1
#> 
#> 
```
