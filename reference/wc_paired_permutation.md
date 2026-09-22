# Paired label-permutation test for Wang–Chang recurrent gap-time curves

Tests a two-visit paired comparison by repeatedly swapping the two visit
labels within each pair and recalculating
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md).
This supplies randomization-based inference when within-person visits
are paired; it is not an independent-arm rank test.

## Usage

``` r
wc_paired_permutation(
  data,
  pair,
  group,
  id = "id",
  gap = "gap",
  status = "status",
  episode = NULL,
  rho = 0,
  tau = Inf,
  B = 999,
  seed = NULL
)
```

## Arguments

- data:

  Long recurrent gap-time data.

- pair:

  Character scalar naming the paired-participant column. Each pair must
  have exactly two distinct visit-level `id` values.

- group:

  Character scalar naming the two visit labels.

- id, gap, status, episode, rho, tau:

  Passed to
  [`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md).

- B:

  Number of random within-pair label permutations.

- seed:

  Optional random seed.

## Value

A list with the observed Wang–Chang rank fit, permutation z statistics,
and a two-sided randomization p-value.

## Details

The procedure requires one observation from each of the two labels
within every pair. It tests the sharp within-pair label-exchangeability
null. The usual
[`wc_logrank()`](https://muschellij2.github.io/recurSurvTests/reference/wc_logrank.md)
p-value is also returned for comparison, but it does not itself account
for pairing of visits.

## Examples

``` r
dat <- data.frame(
  id = rep(1:4, each = 2), pair = rep(1:2, each = 4),
  episode = rep(1:2, 4), gap = c(2, 4, 3, 2, 1, 5, 2, 3),
  status = rep(c(1, 0), 4), visit = rep(c("A", "B"), each = 2)
)
wc_paired_permutation(dat, pair = "pair", group = "visit", episode = "episode",
  B = 99, seed = 1)
#> $method
#> [1] "paired within-pair label permutation for Luo-Huang WRS G_rho*"
#> 
#> $rho
#> [1] 0
#> 
#> $B
#> [1] 99
#> 
#> $observed
#> $observed$method
#> [1] "Luo-Huang WRS G_rho*"
#> 
#> $observed$rho
#> [1] 0
#> 
#> $observed$interpretation
#> [1] "recurrent gap-time log-rank analogue"
#> 
#> $observed$group0
#> [1] "A"
#> 
#> $observed$group1
#> [1] "B"
#> 
#> $observed$n_subjects
#> [1] 4
#> 
#> $observed$G
#> [1] -0.2083333
#> 
#> $observed$se
#> [1] 0.1301041
#> 
#> $observed$z
#> [1] -1.601282
#> 
#> $observed$chisq
#> [1] 2.564103
#> 
#> $observed$df
#> [1] 1
#> 
#> $observed$p.value
#> [1] 0.1093146
#> 
#> $observed$subject_influence
#>   id         psi event_score
#> 1  1 -0.09722222  -0.6666667
#> 2  2 -0.34722222   0.0000000
#> 3  3 -0.37500000  -0.5000000
#> 4  4 -0.01388889   0.3333333
#> 
#> $observed$detail
#>   time risk risk_group1 expected_group_fraction dN  dLambda0 surv_left weight
#> 1    1    4           2               0.5000000  1 0.2500000      1.00      1
#> 2    2    3           2               0.6666667  2 0.6666667      0.75      1
#> 3    3    1           1               1.0000000  1 1.0000000      0.25      1
#> 
#> 
#> $permutation_z
#>  [1] -1.6012815 -1.6012815 -0.2443389 -1.6012815 -0.2443389 -1.6012815
#>  [7] -0.2443389 -0.2443389 -1.6012815 -0.2443389 -0.2443389 -0.2443389
#> [13] -1.6012815 -1.6012815 -0.2443389 -0.2443389 -1.6012815 -1.6012815
#> [19] -0.2443389 -0.2443389 -1.6012815 -1.6012815 -1.6012815 -1.6012815
#> [25] -1.6012815 -0.2443389 -1.6012815 -1.6012815 -0.2443389 -0.2443389
#> [31] -0.2443389 -1.6012815 -0.2443389 -0.2443389 -0.2443389 -0.2443389
#> [37] -1.6012815 -0.2443389 -0.2443389 -1.6012815 -0.2443389 -1.6012815
#> [43] -0.2443389 -0.2443389 -1.6012815 -1.6012815 -1.6012815 -1.6012815
#> [49] -1.6012815 -1.6012815 -0.2443389 -0.2443389 -0.2443389 -1.6012815
#> [55] -1.6012815 -1.6012815 -1.6012815 -1.6012815 -0.2443389 -0.2443389
#> [61] -0.2443389 -1.6012815 -0.2443389 -0.2443389 -0.2443389 -0.2443389
#> [67] -0.2443389 -1.6012815 -1.6012815 -1.6012815 -1.6012815 -1.6012815
#> [73] -0.2443389 -0.2443389 -0.2443389 -1.6012815 -1.6012815 -0.2443389
#> [79] -0.2443389 -1.6012815 -0.2443389 -0.2443389 -0.2443389 -1.6012815
#> [85] -0.2443389 -1.6012815 -0.2443389 -0.2443389 -1.6012815 -1.6012815
#> [91] -1.6012815 -1.6012815 -0.2443389 -1.6012815 -1.6012815 -0.2443389
#> [97] -0.2443389 -0.2443389 -0.2443389
#> 
#> $p.value
#> [1] 0.49
#> 
```
