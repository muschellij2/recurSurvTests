# recurSurvTests

[![Codecov test coverage](https://codecov.io/gh/muschellij2/recurSurvTests/graph/badge.svg)](https://app.codecov.io/gh/muschellij2/recurSurvTests)

`recurSurvTests` provides transparent, research-oriented reference
implementations for recurrent **gap-time** survival and recurrent competing-risk
analyses.

It deliberately separates marginal recurrent gap-time survival, mean cumulative
recurrent-event burden over calendar time (not implemented here), and recurrent
cause-specific cumulative incidence.

## Methods

- `wc_surv()` estimates Wang-Chang marginal recurrent gap-time survival using
  the equivalent Luo-Huang weighted risk set.
- `psh_surv()` estimates the Pena-Strawderman-Hollander generalized
  product-limit curve (the pooled gap-record Kaplan-Meier curve).
- `wc_logrank()` implements Luo-Huang's two-sample `G_rho*` statistic. With
  `rho = 0`, it is a recurrent gap-time log-rank analogue; it is not an exact
  Zhao et al. (2020) implementation.
- `zhao_rank_test()` separately implements Zhao et al.'s extended LR/GB/PP
  score construction. `zhao_rank_simulation()` provides paper and PSG-like
  null-calibration scenarios.
- `ss_rcif()` estimates Sivadasan-Sankaran recurrent cause-specific cumulative
  incidence functions (RCIFs).
- `ss_rcif_equal_causes_test()` is a subject-bootstrap Wald reconstruction of
  the Sivadasan-Sankaran/METRON equal-cause RCIF test. It is a weighted CIF
  contrast, not a log-rank test.

## Installation

```r
# install.packages("pak")
pak::pak("muschellij2/recurSurvTests")
```

## Data convention

Supply one row per observed gap, ordered within subject by `episode` (or already
in episode order). The final row for each subject is a right-censored gap with
`status = 0`; all prior rows have `status = 1`. `gap` is a duration, not a
calendar time. For RCIF analyses, completed events additionally have `cause`.

See `vignette("recurrent-methods", package = "recurSurvTests")` for the
estimands, formulas, implementation details, and a worked example.

## Relationship to existing software

This package is a transparent reference collection, not a claim that every
component is a new algorithm. `newTestSurvRec` already provides Wang--Chang and
Pena--Strawderman--Hollander recurrent-curve machinery, so `wc_surv()` and
`psh_surv()` are compatibility/reference implementations. `reda` and `mets`
focus on calendar-time mean cumulative recurrent-event functions; those are a
different estimand and are deliberately not implemented here. The package's
distinct scope is the Luo--Huang `G_rho*` two- and K-group tests, the separately
labeled Zhao score/variance variants and simulations, and the
Sivadasan--Sankaran recurrent cause-specific incidence methods. See the
[full scope audit](docs/related-software.md) for the comparison and sources.

In practical terms, use `reda::mcfDiff()` or
`mets::recurrent_marginal()`/`mets::test_logrankRecurrent()` for calendar-time
marginal event means and their AUC/logrank-type comparisons. Use this package
when the target is a marginal **gap-time survival distribution**, a
Luo--Huang/Zhao recurrent-gap rank contrast, or a recurrent competing-cause
incidence curve. The package is intended to make those distinctions and their
variance choices inspectable, not to duplicate the mean-function ecosystem.

The Zhao implementation is particularly useful because it applies LR,
Gehan--Breslow, and Peto--Prentice weighting to Wang--Chang survival curves
while retaining subject-level recurrent-gap dependence through robust residual
covariance. It therefore answers a different question from a mean-function
test, and exposes the paper's competing variance conventions for sensitivity
analysis rather than hiding them behind a generic “recurrent log-rank” label.

## Example

```r
library(recurSurvTests)

gaps <- data.frame(
  id = rep(1:4, each = 3),
  episode = rep(1:3, 4),
  gap = c(2, 3, 4, 1, 4, 3, 2, 2, 5, 1, 3, 4),
  status = rep(c(1, 1, 0), 4),
  arm = rep(c("A", "A", "B", "B"), each = 3),
  cause = c("A", "B", NA, "A", "A", NA, "B", "A", NA, "B", "B", NA)
)

wc_surv(gaps, episode = "episode")
wc_logrank(gaps, group = "arm", episode = "episode", rho = 0)$p.value
ss_rcif(gaps, episode = "episode")$curve
```

## References

- Wang MC, Chang SH (1999). *JASA* 94, 146-153.
- Luo X, Huang CY (2011). *Statistics in Medicine* 30, 301-311.
  <https://doi.org/10.1002/sim.4074>
- Sivadasan SM, Sankaran PG (2023). *Statistica* 83, 3-25.
- Sivadasan SM, Sankaran PG (2022). *METRON*.
  <https://doi.org/10.1007/s40300-022-00228-x>
