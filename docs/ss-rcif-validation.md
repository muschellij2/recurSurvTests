# RCIF finite-sample validation pilot

The [Sivadasan--Sankaran 2023 paper](https://rivista-statistica.unibo.it/article/download/11727/18151/81363) specifies exponential frailty with mean one, two cause-specific Weibull hazards with common shape 2 and coefficients 0.1 and 0.12, and both fixed and uniform censoring limits of 3 or 5. Under this model, the marginal cause-specific CIF has a closed form:

\[
F_l(t) = \frac{\lambda_l}{\lambda_1+\lambda_2}
\left(1-\frac{1}{1+(\lambda_1+\lambda_2)t^2}\right).
\]

Run the small pilot from the package root:

```sh
Rscript data-raw/ss-rcif-quick-simulation.R 100 100 3
```

Arguments are replications, subjects, and censoring limit (3 or 5). The script reports mean estimate, signed bias, and MSE for each cause at gap times 1 and 2 under fixed and uniform censoring. It uses the same generated dataset for all three estimators: literal Equation 15 with paper-style integration, shared weighting with paper-style integration, and shared weighting with exact exponential-decrement allocation. For more stable results, increase the replications and run each of the paper's sample sizes 50, 100, 200, and 500 under both limits.

In a 100-replication pilot with 100 subjects and censoring limit 3, the cause A bias at time 2 was -0.137 for literal Equation 15, -0.0017 for shared weighting with paper-style integration, and -0.0006 for shared weighting with exact decrement under fixed censoring. Under uniform censoring, the corresponding biases were -0.182, 0.0019, and 0.0044. These are exploratory Monte Carlo results, not a reproduction of the paper's tables. They indicate that the Equation 15 singleton exclusion warrants investigation; they do not resolve whether the printed formula is a typo or establish general consistency of the shared-weight estimator.
