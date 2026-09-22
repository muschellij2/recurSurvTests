# Restricted mean recurrent gap time with a subject-bootstrap confidence interval

Estimates the area under a Wang–Chang or PSH recurrent gap-time survival
curve through `tau`, and obtains a percentile confidence interval by
resampling complete subject histories.

## Usage

``` r
recurrent_surv_rmst_ci(
  data,
  estimator = c("wc", "psh"),
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  tau,
  conf.level = 0.95,
  B = 999,
  seed = NULL
)
```

## Arguments

- data, id, gap, status, episode:

  Passed to
  [`wc_surv()`](https://muschellij2.github.io/recurSurvTests/reference/wc_surv.md)
  or
  [`psh_surv()`](https://muschellij2.github.io/recurSurvTests/reference/psh_surv.md).

- estimator:

  Either `"wc"` or `"psh"`.

- tau:

  Positive restriction time.

- conf.level:

  Confidence level.

- B:

  Number of subject bootstrap resamples; 999 is the default and 2000 or
  more is recommended for stable two-sided 95 percent tails.

- seed:

  Optional random seed.

## Value

A list with the restricted mean `estimate`, percentile limits, fitted
curve, and bootstrap diagnostics.

## Details

The restricted mean gap time is `integral_0^tau S(u) du`. It remains
estimable when a survival median is not reached. Resampling is by
subject, preserving recurrent-gap clustering and WC weights.
