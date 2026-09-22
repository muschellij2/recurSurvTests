# Subject-bootstrap confidence intervals for two recurrent survival curves

Compares two WC or PSH curves at requested gap times using a cluster
bootstrap. Independent arms are resampled within arm; when `pair` is
given, complete paired clusters are resampled, preserving the paired
visits.

## Usage

``` r
recurrent_surv_contrast_ci(
  data,
  group,
  pair = NULL,
  estimator = c("wc", "psh"),
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  tau = Inf,
  times = NULL,
  contrast = c("difference", "ratio", "log_ratio"),
  interval = c("pointwise", "simultaneous"),
  conf.level = 0.95,
  B = 999,
  seed = NULL
)
```

## Arguments

- data:

  Long recurrent gap-time data.

- group:

  Subject-level two-group column.

- pair:

  Optional paired-cluster column.

- estimator, id, gap, status, episode, tau:

  Passed to curve estimators.

- times:

  Evaluation times; pooled observed event times by default.

- contrast:

  `"difference"`, `"ratio"`, or `"log_ratio"`, defined as group 1
  minus/divided by group 0.

- interval:

  Pointwise percentile intervals or simultaneous max-t bands.

- conf.level, B, seed:

  Usual bootstrap controls.

## Value

Observed group curves, a contrast curve with limits, and bootstrap
diagnostics.
