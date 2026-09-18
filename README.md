# recurSurvTests

A small research-oriented R package skeleton containing reference implementations
for recurrent gap-time survival and recurrent competing-risk methods discussed in
our September 2026 analysis.

## Included methods

- `wc_surv()` — Wang-Chang / Luo-Huang weighted-risk-set estimator of the
  marginal recurrent gap-time survival curve.
- `wc_logrank()` — Luo-Huang `G_rho*` two-sample rank test. `rho = 0` is the
  recurrent gap-time **log-rank analogue**; `rho = 1` is the recurrent
  Peto-Prentice analogue.
- `ss_rcif()` — Sivadasan-Sankaran recurrent cumulative-incidence estimator.
- `ss_rcif_equal_causes_test()` — subject-bootstrap Wald reconstruction of the
  Sivadasan-Sankaran / METRON equal-cause RCIF test.

## Important status note about Zhao et al. (2020)

The paper *Extended Rank Tests for Analyzing Recurrent Event Data* explicitly
extends rank tests to Wang-Chang survival curves and uses robust subject-level
variance estimation. The supplementary DOCX supplied in the source chat contains
simulation details and residual diagnostics, **not implementation code or the
full test formulas**. Therefore this package does not label any function as an
exact Zhao et al. implementation. `wc_logrank()` implements the independently
published Luo-Huang Wang-Chang-equivalent recurrent rank test.

If the full Zhao et al. paper/formulas are added later, implement them in a
separate `R/zhao_rank_tests.R` file and validate against the paper's simulation
settings rather than silently treating `wc_logrank()` as identical.

## Install locally

From the parent directory containing `recurSurvTests/`:

```r
install.packages(c("roxygen2", "devtools"))
roxygen2::roxygenise("recurSurvTests")
devtools::install("recurSurvTests")
```

Or open the package directory as an RStudio project and run:

```r
roxygen2::roxygenise()
devtools::document()
devtools::check()
devtools::install()
```

Before sharing publicly, change the placeholder email, GitHub URL, and
BugReports URL in `DESCRIPTION`.

## Data convention

One row per observed gap. Within each subject, rows are in episode order. The
final row is a right-censored gap (`status = 0`); preceding rows are completed
recurrent events (`status = 1`). `gap` is a gap duration, not calendar time.
For RCIF analyses, `cause` is present on completed-event rows.

## References

- Wang MC, Chang SH (1999). *JASA* 94:146-153.
- Luo X, Huang CY (2011). *Statistics in Medicine* 30:301-311.
  DOI: 10.1002/sim.4074.
- Zhao Q, Zhang B, LaValley MP, Massaro JM, Lunetta KL, Chang M (2020).
  *Statistics in Biopharmaceutical Research* 12(1):90-98.
  DOI: 10.1080/19466315.2019.1601596.
- Sivadasan SM, Sankaran PG (2023). *Statistica* 83(1):3-25.
- Sivadasan SM, Sankaran PG (2022). *METRON*.
  DOI: 10.1007/s40300-022-00228-x.

## Validation status

The code has been statically reviewed against the formulas used to construct the
reference implementation, but the current execution environment did not contain
R, so `R CMD check` was not run here. Treat this as a package starting point and
run `devtools::check()` before analysis or publication use.
