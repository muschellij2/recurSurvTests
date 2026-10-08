# K-group null calibration diagnostic

The existing 20-task study gives a precise rejection estimate above 0.05 for
the 30/45/60-group PSG-shaped null. This diagnostic checks whether the excess
shrinks with more independent subjects or changes when the same total sample
size is allocated evenly. It also records an F-reference calculation as a
finite-sample diagnostic; this is not proposed as a replacement test unless
its behavior is supported by broader validation.

The four designs are the original unbalanced N=135, balanced N=135, unbalanced
N=270, and unbalanced N=540. Each uses the original frailty, burst/recovery,
recording duration, and terminal censoring generator under the global null.
Each design has 5,000 independent simulation replicates. The primary output
contains the documented chi-square p-value, test statistic, and degrees of
freedom; it also contains the exploratory finite-sample F-reference p-value.

Run from the repository root with:

```sh
sbatch data-raw/wc-logrank-k-diagnostics/simulation.sbatch
```

Then summarize after all four tasks finish:

```sh
Rscript data-raw/wc-logrank-k-diagnostics/aggregate.R
```

The run is intentionally separate from the previous K-group output. Task RDS
files are written under `results/`; no summary is updated until all four tasks
are present. Compare empirical size and confidence intervals across the four
designs. If the size approaches 0.05 as N grows, finite-sample Wald calibration
is implicated. If it remains elevated, investigate the score/influence
covariance derivation and the effect of the within-subject state process.
