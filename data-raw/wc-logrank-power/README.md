# Two-group recurrent gap-time power study

This study estimates the size and power of the Luo-Huang weighted-risk-set
tests `wc_logrank(rho = 0)` and `wc_logrank(rho = 1)` as functions of sample
size and treatment gap-time effect. It is separate from the PSG validation and
does not claim to reproduce Zhao et al.'s simulation tables.

Each data set has equal independent groups, subject-specific multiplicative
heterogeneity, staggered entry, and administrative censoring at 180 time units.
Baseline gaps are exponential or Weibull. The subject multiplier is uniform on
`(0.5, 1.5)` (low heterogeneity) or `(0.01, 1.99)` (high heterogeneity).
Treatment gaps are multiplied by `exp(delta)` for
`delta = -0.35, -0.20, 0, 0.20, 0.35`; `delta = 0` is the size cell, positive
values mean longer treatment gaps, and negative values mean shorter gaps.
Each arm has 25, 50, 100, or 200 subjects. Both tests use the same generated
data set in every replicate. The default is 5,000 replicates per cell, giving a
Monte Carlo standard error of about 0.003 at a true rejection probability of
0.05 and about 0.007 at power 0.5.

Run the 20-task array from the repository root:

```sh
sbatch data-raw/wc-logrank-power/simulation.sbatch
```

After all tasks finish, aggregate them:

```sh
Rscript data-raw/wc-logrank-power/aggregate.R
```

Task RDS files and logs stay under this directory. The aggregate writes
`results.csv`, `results.rds`, and `power-curves.png` here. The task IDs are
assigned to design cells deterministically; each task output records its
simulation settings and raw p-values so the summaries can be independently
recomputed.
