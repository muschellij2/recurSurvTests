# Subject-bootstrap confidence intervals for recurrent survival curves

[`recurrent_surv_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_ci.md)
supplies subject-level bootstrap uncertainty intervals for the
Wang–Chang (WC) marginal gap-time curve and the
Pena–Strawderman–Hollander (PSH) curve. The resampling unit is the
complete participant history. Resampling individual rows would
incorrectly regard a participant’s recurrent gaps as independent and,
for WC, would break the subject-specific event-count weights.

## Relationship to earlier R implementations

Earlier R software provides valuable, but different, uncertainty tools.
`survrec::wc.fit()` and `survrec::psh.fit()` return pointwise standard
errors; the associated plotting method uses them to draw
normal-approximation pointwise limits. `survrec::survdiffr()`
additionally creates WC or PSH bootstrap replicates for a selected
survival quantile and returns a `boot` object, from which users can
request a quantile confidence interval. Therefore the present function
does not introduce subject bootstrap resampling as a new statistical
idea. Instead, it makes full-curve subject bootstrap intervals available
with the same call for WC and PSH, states whether the result is
pointwise or simultaneous, and returns the fitted curve and limits in a
common format. The inspected `newTestSurvRec` WC and PSH fitters provide
standard errors but no exported bootstrap interface.

Those legacy standard errors use the independence-style
Nelson–Aalen/Greenwood calculation

``` math
\widehat{\mathrm{SE}}\{\widehat S(t)\} =
\widehat S(t)\left\{\sum_{u\le t}dN(u)/R(u)^2\right\}^{1/2}.
```

They can produce a normal pointwise curve CI, but they do not produce a
simultaneous band and are not generally adequate for correlated
recurrent gaps. `recurrent_surv_ci(..., ci_method = "normal")` exposes
this calculation using this package’s own risk increments. It reproduces
the PSH calculation on the same long records; for WC it is an analogous
weighted-risk-set approximation rather than a claim to reproduce another
package’s variance code. Its default subject bootstrap is preferable
when participant-level dependence is relevant.

## Pointwise intervals are the usual curve summary

The default is a percentile interval at each event time. A 95%
*pointwise* interval has approximately 95% repeated-sampling coverage at
one pre-specified time. This is the usual uncertainty display alongside
a survival curve; it does not mean that the full curve is covered with
probability 0.95.

``` r

gaps <- data.frame(
  id = rep(1:12, each = 3),
  episode = rep(1:3, 12),
  gap = rep(c(1, 2, 4), 12) + rep(0:11, each = 3) / 10,
  status = rep(c(1, 1, 0), 12)
)

wc_ci <- recurrent_surv_ci(
  gaps, estimator = "wc", episode = "episode", B = 199, seed = 1
)
head(wc_ci$curve)
#>   time      surv  std.error     lower     upper
#> 1  1.0 0.9583333 0.03994660 0.8750000 1.0000000
#> 2  1.1 0.9166667 0.05233246 0.8333333 1.0000000
#> 3  1.2 0.8750000 0.06121266 0.7500000 1.0000000
#> 4  1.3 0.8333333 0.06914942 0.6798611 0.9583333
#> 5  1.4 0.7916667 0.07317706 0.6666667 0.9166667
#> 6  1.5 0.7500000 0.07373072 0.6250000 0.8750000
```

The returned `curve` has the estimated survival, bootstrap standard
error, and lower and upper limits. `B = 999` or larger is recommended
for analysis; `B = 199` only keeps this vignette quick. The default
`B = 999` follows common bootstrap and permutation practice. For a
two-sided 95% interval, however, its 2.5% and 97.5% limits are
determined by only about 25 bootstrap draws in each tail. Use `B = 2000`
or more when stable 95% percentile limits or max-t band critical values
are important and runtime permits.

For a legacy-style normal pointwise interval, use the standard-error
option. It is shown here for comparison, not as the default
recurrent-event analysis.

``` r

psh_normal_ci <- recurrent_surv_ci(
  gaps, estimator = "psh", episode = "episode", ci_method = "normal"
)
head(psh_normal_ci$curve)
#>   time      surv  std.error     lower     upper
#> 1  1.0 0.9722222 0.02700617 0.9192911 1.0000000
#> 2  1.1 0.9444444 0.03763503 0.8706811 1.0000000
#> 3  1.2 0.9166667 0.04540030 0.8276837 1.0000000
#> 4  1.3 0.8888889 0.05161113 0.7877329 0.9900449
#> 5  1.4 0.8611111 0.05677994 0.7498245 0.9723978
#> 6  1.5 0.8333333 0.06117145 0.7134395 0.9532272
```

## Simultaneous bands answer a different question

If an analysis makes a statement about the whole displayed curve,
request a simultaneous band. The implementation uses the bootstrap
distribution of the largest absolute, standard-error-scaled curve
deviation over the requested times. It therefore controls coverage
jointly over that finite time grid, and is generally wider than
pointwise intervals.

``` r

psh_band <- recurrent_surv_ci(
  gaps, estimator = "psh", episode = "episode", interval = "simultaneous",
  B = 199, seed = 1
)
psh_band$critical.value
#> [1] 2.827631
head(psh_band$curve)
#>   time      surv  std.error     lower     upper
#> 1  1.0 0.9722222 0.02663107 0.8969194 1.0000000
#> 2  1.1 0.9444444 0.03488831 0.8457932 1.0000000
#> 3  1.2 0.9166667 0.04080844 0.8012755 1.0000000
#> 4  1.3 0.8888889 0.04609961 0.7585362 1.0000000
#> 5  1.4 0.8611111 0.04878471 0.7231660 0.9990562
#> 6  1.5 0.8333333 0.04915381 0.6943445 0.9723222
```

For a band restricted to clinically meaningful gap times, pass them
explicitly via `times`. This both specifies the scientific claim and
avoids creating a wide band over poorly supported tail times.

## Which confidence statement is being made?

These intervals are not interchangeable.

| Quantity | Typical output | Coverage statement |
|----|----|----|
| Curve ordinate at one pre-specified time, $`S(t_0)`$ | Default bootstrap percentile interval | About 95% coverage at $`t_0`$ |
| Several separately reported curve ordinates | Pointwise intervals at each time | About 95% at each time separately; not 95% for all times together |
| Entire curve over supplied grid $`mathcal T`$ | `interval = "simultaneous"` max-t band | About 95% coverage of all $`S(t)`$ for $`t \in \mathcal T`$ jointly |
| Survival quantile, such as a median gap | [`recurrent_surv_quantile_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_quantile_ci.md) subject-bootstrap interval | Coverage for an inverse-curve quantity, not for $`S(t)`$ |
| Difference between two group curves | Not returned by [`recurrent_surv_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_ci.md) | Requires a two-group, cluster-resampling design and a specified contrast |

The normal-SE option is only a pointwise approximation. It cannot become
a simultaneous band by drawing many normal intervals, and it can
understate uncertainty when recurrent gaps are correlated. A
simultaneous band is useful for a global visual or scientific claim, but
its width means pointwise intervals are usually preferable for reporting
a small number of pre-specified times.

## Quantile statistics, including the median

The median gap time is the inverse-curve statistic
$`\inf\{t:S(t)\le 0.5\}`$. It has a different coverage target from a CI
for the survival ordinate $`S(t)`$.
[`recurrent_surv_quantile_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_quantile_ci.md)
resamples complete participant histories, refits the selected curve, and
takes percentile limits of the resulting median (or other quantile)
estimates.

``` r

median_ci <- recurrent_surv_quantile_ci(
  gaps, estimator = "wc", episode = "episode", survival_prob = 0.5,
  B = 199, seed = 1
)
median_ci[c("estimate", "lower", "upper")]
#> $estimate
#> [1] 2.1
#> 
#> $lower
#> [1] 1.9
#> 
#> $upper
#> [1] 2.2
```

The function deliberately does not invert the normal/Greenwood curve
limits for this purpose: those limits do not account for clustered
recurrent gaps.

## Restricted means and two-group contrasts

When the curve does not reach 0.5, the restricted mean gap time
$`\int_0^\tau S(t)dt`$ remains defined. It and its subject-bootstrap
interval are available through
[`recurrent_surv_rmst_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_rmst_ci.md).

``` r

recurrent_surv_rmst_ci(gaps, estimator = "wc", episode = "episode",
                       tau = 3, B = 199, seed = 1)
#> $method
#> [1] "Wang-Chang subject-bootstrap restricted mean gap time"
#> 
#> $estimator
#> [1] "wc"
#> 
#> $tau
#> [1] 3
#> 
#> $conf.level
#> [1] 0.95
#> 
#> $estimate
#> [1] 2.045833
#> 
#> $lower
#> [1] 1.830278
#> 
#> $upper
#> [1] 2.240556
#> 
#> $fit
#>    time risk  dN    dLambda  surv_left       surv
#> 1   1.0 12.0 0.5 0.04166667 1.00000000 0.95833333
#> 2   1.1 11.5 0.5 0.04347826 0.95833333 0.91666667
#> 3   1.2 11.0 0.5 0.04545455 0.91666667 0.87500000
#> 4   1.3 10.5 0.5 0.04761905 0.87500000 0.83333333
#> 5   1.4 10.0 0.5 0.05000000 0.83333333 0.79166667
#> 6   1.5  9.5 0.5 0.05263158 0.79166667 0.75000000
#> 7   1.6  9.0 0.5 0.05555556 0.75000000 0.70833333
#> 8   1.7  8.5 0.5 0.05882353 0.70833333 0.66666667
#> 9   1.8  8.0 0.5 0.06250000 0.66666667 0.62500000
#> 10  1.9  7.5 0.5 0.06666667 0.62500000 0.58333333
#> 11  2.0  7.0 1.0 0.14285714 0.58333333 0.50000000
#> 12  2.1  6.0 1.0 0.16666667 0.50000000 0.41666667
#> 13  2.2  5.0 0.5 0.10000000 0.41666667 0.37500000
#> 14  2.3  4.5 0.5 0.11111111 0.37500000 0.33333333
#> 15  2.4  4.0 0.5 0.12500000 0.33333333 0.29166667
#> 16  2.5  3.5 0.5 0.14285714 0.29166667 0.25000000
#> 17  2.6  3.0 0.5 0.16666667 0.25000000 0.20833333
#> 18  2.7  2.5 0.5 0.20000000 0.20833333 0.16666667
#> 19  2.8  2.0 0.5 0.25000000 0.16666667 0.12500000
#> 20  2.9  1.5 0.5 0.33333333 0.12500000 0.08333333
#> 21  3.0  1.0 0.5 0.50000000 0.08333333 0.04166667
#> 
#> $bootstrap_success
#> [1] 199
#> 
#> $bootstrap_failures
#> [1] 0
```

For independent treatment arms,
[`recurrent_surv_contrast_ci()`](https://muschellij2.github.io/recurSurvTests/reference/recurrent_surv_contrast_ci.md)
resamples subjects within arm; supplying `pair` instead resamples
complete paired visit clusters. It provides pointwise differences,
ratios, or log-ratios, and a simultaneous band when requested. These are
curve-difference estimands, not the rank-test statistic.

``` r

gaps$arm <- rep(c("control", "treatment"), each = 18)
contrast_ci <- recurrent_surv_contrast_ci(
  gaps, group = "arm", estimator = "wc", episode = "episode",
  contrast = "difference", B = 199, seed = 1
)
plot_recurrent_surv_ci(wc_ci)
```

![](recurrent-survival-confidence-intervals_files/figure-html/contrast-1.png)

## Coverage check

`data-raw/curve-ci-coverage-simulation.R` is an opt-in simulation
runner. It generates recurrent histories with administrative terminal
censoring and unequal observed event counts, obtains a large-population
reference curve for each estimator, and evaluates pointwise and joint
coverage separately. The large-population reference is intentional:
under terminal-gap recurrent data, the pooled PSH curve need not target
an IID exponential gap survival function.

``` r

Sys.setenv(RUN_CURVE_CI_COVERAGE = "true", CURVE_CI_N_SIM = "1000",
           CURVE_CI_B = "999")
source("data-raw/curve-ci-coverage-simulation.R")

out <- readRDS("data-raw/curve-ci-coverage-results.rds")
out$pointwise
out$simultaneous
```

The standard analysis should report pointwise intervals for fixed-time
curve estimates and use simultaneous bands only when making a joint,
whole-curve claim. The simulation output must be reviewed with
Monte-Carlo standard errors before claiming nominal coverage for a
chosen setting.
