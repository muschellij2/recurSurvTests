# recurSurvTests 0.0.1

* Corrected the paper-scenario log-logistic gap generator to match Zhao
  Table 1 (scale exp(4), inverse-survival exponent 1 before the square root).
  Earlier log-logistic simulation results must be regenerated. Other baseline
  distributions and the PSG generator are unchanged. Array subset submissions
  now retain the original 120-task partition, and outputs identify the corrected
  generator for validation reporting.
* Initial CRAN submission.
* added a number of methods and compared those methods to `newTestSurvRec` and `survrec`.
