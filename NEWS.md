# recurSurvTests 0.0.1

* Added an explicit related-software scope audit. The documentation now
  identifies the overlap of `wc_surv()` and `psh_surv()` with `newTestSurvRec`
  and distinguishes the calendar-time mean-function methods in `reda` and
  `mets` from this package's gap-time and recurrent-CIF targets.
* Corrected the paper-scenario log-logistic gap generator to match Zhao
  Table 1 (scale exp(4), inverse-survival exponent 1 before the square root).
  Earlier log-logistic simulation results must be regenerated. Other baseline
  distributions and the PSG generator are unchanged. Array subset submissions
  now retain the original 120-task partition, and outputs identify the corrected
  generator for validation reporting.
* Initial CRAN submission.
* added a number of methods and compared those methods to `newTestSurvRec` and `survrec`.
