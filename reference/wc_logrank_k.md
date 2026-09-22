# K-group Luo-Huang weighted-risk-set recurrent gap-time rank test

Extends the two-group
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
score-vector construction to `K >= 2` groups. It uses subject-level
influence vectors and their empirical covariance so repeated gaps remain
clustered. For two groups it is an omnibus version of
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md).

## Usage

``` r
wc_logrank_k(
  data,
  group,
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  rho = 0,
  tau = Inf
)
```

## Arguments

- data, group, id, gap, status, episode, rho, tau:

  As for
  [`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md).

## Value

An omnibus chi-square test, score vector, covariance, and event-time
details.

## References

Luo X, Huang CY (2011). Analysis of recurrent gap time data using the
weighted risk-set method. *Statistics in Medicine*, 30, 301-311.
