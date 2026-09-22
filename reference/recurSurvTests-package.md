# recurSurvTests: recurrent survival and recurrent competing-risk reference methods

`recurSurvTests` contains transparent reference implementations for
recurrent gap-time survival and recurrent competing-risk analyses
discussed in the Wang-Chang, Luo-Huang, and Sivadasan-Sankaran papers.

## Details

The package currently exposes recurrent survival, rank-test, and
competing-risk reference functions:

- [`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
  for a Wang-Chang / weighted-risk-set marginal recurrent gap-time
  survival curve.

- [`psh_surv()`](https://muschellij2.github.io/recurSurvTests/reference/psh_surv.md)
  for the Pena–Strawderman–Hollander generalized product-limit recurrent
  survival curve (the pooled gap-record Kaplan–Meier curve).

- [`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
  for the Luo-Huang `G_rho*` two-sample rank test; `rho = 0` is the
  recurrent gap-time log-rank analogue.

- [`ss_rcif()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif.md)
  for the Sivadasan-Sankaran recurrent cumulative-incidence estimator.

- [`ss_rcif_equal_causes_test()`](https://muschellij2.github.io/recurSurvTests/reference/ss_rcif_equal_causes_test.md)
  for a subject-bootstrap Wald reconstruction of the equal-cause
  recurrent cumulative-incidence test.

## Data convention

Input data are in long format with one row per observed gap. Within each
subject, rows are in episode order. The final row is the right-censored
gap (`status = 0`), while all preceding rows are completed recurrent
events (`status = 1`). `gap` is a gap duration rather than calendar
time. For recurrent competing-risk analyses, `cause` is defined on
completed-event rows.

## Important implementation note

These are research reference implementations assembled from published
formulas. In particular, the equal-cause RCIF test is explicitly a
bootstrap reconstruction rather than claimed byte-for-byte reproduction
of unpublished author software.

## References

Wang MC, Chang SH (1999). Nonparametric estimation of a recurrent
survival function. *Journal of the American Statistical Association*,
94, 146-153.

Pena EA, Strawderman RL, Hollander M (2001). Nonparametric estimation
with recurrent event data. *Journal of the American Statistical
Association*, 96, 1299-1315.
[doi:10.1198/016214501753381887](https://doi.org/10.1198/016214501753381887)
.

Luo X, Huang CY (2011). Analysis of recurrent gap time data using the
weighted risk-set method and the modified within-cluster resampling
method. *Statistics in Medicine*, 30, 301-311.
[doi:10.1002/sim.4074](https://doi.org/10.1002/sim.4074) .

Sivadasan SM, Sankaran PG (2023). Nonparametric estimation of cumulative
incidence functions of recurrent events. *Statistica*, 83(1), 3-25.

Sivadasan SM, Sankaran PG (2022). A nonparametric test for comparing
recurrent cumulative incidence functions. *METRON*.
[doi:10.1007/s40300-022-00228-x](https://doi.org/10.1007/s40300-022-00228-x)
.

## See also

Useful links:

- <https://muschellij2.github.io/recurSurvTests/>

- <https://github.com/muschellij2/recurSurvTests>

- Report bugs at <https://github.com/muschellij2/recurSurvTests/issues>

## Author

**Maintainer**: John Muschelli <muschellij2@gmail.com>

Authors:

- John Muschelli <muschellij2@gmail.com>
