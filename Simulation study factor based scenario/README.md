## Simulation Study: Factor-Based Scenario

This simulation study evaluates correlation estimation under approximate factor model settings, following the design described in Section titled `Approximate factor model simulations` of the manuscript. Data are generated from a factor model with a low-rank common component and an idiosyncratic covariance matrix, and methods are compared under varying assumptions on factor observability and idiosyncratic structure. The study benchmarks classical and robust POET-based estimators against GLASD-based correlation estimation when the number of factors is unknown.

Simulations are conducted for dimensions \( p = 50 \) and \( p = 100 \), under two regimes for the idiosyncratic covariance matrix \( \Sigma_u \): diagonal and dense random. Estimation accuracy is summarized as reported in the manuscript.

Please follow the steps below to reproduce the results.

1) **Data generation**
Run `Data_generate.m` for `p = 50, 100`. This script generates factor-model data under the specified idiosyncratic covariance regimes.

2) **Naive estimator (unknown factors)**
Run `Method_Naive.m` for `p = 50, 100`, with `Sigma_u_regime = 'diagonal'` and `'random'`.

3) **POET (unknown factors)**
Run `Method_FanEtAl2013.m` for `p = 50, 100`, with `Sigma_u_regime = 'diagonal'` and `'random'`.

4) **Robust POET (known factors)**
Run `Method_FanEtAl2019.m` for `p = 50, 100`, with `Sigma_u_regime = 'diagonal'` and `'random'`.

5) **Robust POET (unknown factors)**
Run `Method_FanEtAl2019_unobservedF.m` for `p = 50, 100`, with `Sigma_u_regime = 'diagonal'` and `'random'`.

6) **GLASD (unknown factors)**
Run `Method_GLASD.m` for `p = 50, 100`, with `Sigma_u_regime = 'diagonal'` and `'random'`.

7) **Summarize results**
After completing all runs, execute `summary_output.m`, which aggregates the results and produces the summary tables reported in the manuscript for this scenario.

**Note.** Originally, the `Data/` folder contained more than 400 generated data files and the `Output/` folder contained more than 600 generated output files. Due to their large size, these files are not included in this repository. Only the code files required to fully reproduce the results are provided.

