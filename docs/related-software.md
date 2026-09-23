# Related software and scope audit

This audit records the boundary between `recurSurvTests` and established R
packages. It is intentionally explicit: the package does not claim that the
Wang--Chang or Pena--Strawderman--Hollander estimators are new algorithms.
Those estimators are included as transparent reference implementations because
they are the estimands needed by the rank tests and confidence intervals here.

| Software | Main target | Relationship to `recurSurvTests` |
| --- | --- | --- |
| `newTestSurvRec` / `TestSurvRec` | Recurrent survival curves and several two-sample recurrent-curve tests, including Wang--Chang and generalized Pena product-limit machinery | Overlaps with `wc_surv()` and `psh_surv()`. Its `Dif.Surv.Rec()` source path calls the PSH/GPLE-style fit and uses event-time risk-set variance; it is not Zhao's Wang--Chang subject-residual Equation 6--7 construction. `wc_logrank()` is the Luo--Huang weighted-risk-set `G_rho*` construction, a different implementation and estimand description. |
| `reda` | Calendar-time mean cumulative functions, `mcfDiff()` and Lawless--Nadeau/Cook--Lawless--Nadeau pseudo-score tests, plus frailty/rate models | Not a duplicate of the gap-time survival or recurrent-CIF targets in this package. Mean cumulative event burden is deliberately out of scope here. The implementation is available from CRAN, R-universe, and [GitHub](https://github.com/wenjie2wang/reda). |
| `mets` | Marginal mean recurrent-event functions and logrank-type tests for those means, including AUC and terminal-event/cumulative-incidence settings | Same calendar-time mean-function target distinction as `reda`; these are not replacements for `wc_logrank()` or `ss_rcif_equal_causes_test()`. The current development source is [GitHub](https://github.com/kkholst/mets) and its package/vignettes are published through R-universe. |
| `survrec` (archived) | Earlier recurrent-survival estimators and bootstrap comparisons | Useful historical cross-check; not a source of the Luo--Huang K-group test, Zhao diagnostics, or Sivadasan--Sankaran RCIF implementation here. |

The package-specific contribution is therefore the documented, estimand-aware
reference collection: Luo--Huang two- and K-group weighted-risk-set rank tests,
the separately labeled Zhao et al. score/variance variants and simulations, and
Sivadasan--Sankaran recurrent cause-specific incidence (including the explicit
`cause_weighting` sensitivity option). Generic bootstrap intervals and the
Wang--Chang/PSH curves should be read as compatibility/reference components,
not as claims of algorithmic novelty.

## Ecosystem fit

For a calendar-time question such as "how many events does a participant
experience by time `t`, on average?", use `reda::mcfDiff()` or
`mets::recurrent_marginal()`/`mets::test_logrankRecurrent()`. These methods are
the appropriate complements when the estimand is a marginal mean or its AUC.
For a gap-time question such as "does the distribution of the next recurrent
gap differ between groups after accounting for unequal numbers of observed
gaps?", use `wc_surv()` and `wc_logrank()` here. For competing recurrent causes,
use `ss_rcif()` and its equal-cause contrast. This separation prevents calling
every recurrent-event comparison a log-rank test when the underlying target is a
mean function or a cause-specific incidence curve.

The Zhao implementation is useful because it occupies a narrow gap between
those ecosystems. It applies familiar LR, Gehan--Breslow, and Peto--Prentice
weighting ideas to Wang--Chang marginal recurrent gap-time survival curves,
then uses subject-level residuals and a robust covariance to retain dependence
among a participant's recurrent gaps. Thus it tests a survival-distribution
contrast, not an event-count mean. It also makes the practical variance choice
visible: this package exposes the pooled-risk residual convention and the
group-specific denominator printed in Zhao's Equation 6 as separate options,
with simulations and validation diagnostics. That transparency is valuable
because the paper's formulas and published tables do not uniquely settle every
implementation detail, and an apparently ordinary "recurrent log-rank" label
can otherwise conceal a different estimand or variance calculation.

Sources checked: [newTestSurvRec CRAN documentation](https://cran.r-project.org/package=newTestSurvRec),
[the CRAN `newTestSurvRec` mirror on GitHub](https://github.com/cran/newTestSurvRec),
[reda CRAN documentation](https://cran.r-project.org/package=reda),
[reda R-universe](https://wenjie2wang.r-universe.dev/reda),
[the reda GitHub repository](https://github.com/wenjie2wang/reda),
[the mets recurrent-events vignette](https://kkholst.github.io/mets/articles/recurrent-events.html),
[mets `test_logrankRecurrent()` documentation](https://kkholst.github.io/mets/reference/test_logrankRecurrent.html),
and [the mets GitHub repository](https://github.com/kkholst/mets).

A targeted search of CRAN, R-universe, and GitHub found maintained software for
the Wang--Chang/PSH curves and calendar-time mean-function analyses, but did
not identify a separate maintained implementation of the Zhao extended
Wang--Chang rank tests or the Sivadasan--Sankaran recurrent-CIF estimator. This
is a practical ecosystem finding rather than a proof that no unpublished or
private code exists.
