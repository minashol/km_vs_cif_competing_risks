# Bias from Treating Competing Events as Censoring

This repository contains a small simulation study comparing the naive Kaplan-Meier estimator with the cumulative incidence function in a competing risks setting.

The event of interest is relapse, while death before relapse is treated as a competing event.

## Main idea

In competing risks settings, treating competing events as ordinary censoring changes the target estimand. The naive Kaplan-Meier estimator may overestimate the real-world cumulative incidence of the event of interest.

## Files

- `KM vs. CIF.R`: R code used to generate the simulation results and figures.
- `KM-vs-CIF.pdf`: short report summarising the simulation design, results, and interpretation.

## Simulation setup

For each simulated individual, two event times are generated:

- \(T_1\): time to relapse,
- \(T_2\): time to death before relapse.

The relapse rate is fixed at

$$
\lambda_1 = 0.08,
$$

while the competing event rate varies over

```{r}
lambda2_values <- c(0, 0.04, 0.08, 0.16)
