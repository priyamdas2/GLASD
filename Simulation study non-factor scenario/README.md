## Simulation Study: Non–Factor-Based Scenario

This simulation study evaluates robust correlation estimation methods in settings without latent factor structure, following the design described in Section 6.1 of the manuscript. The study compares the classical Gaussian estimator, Tyler’s M-estimator, and several loss-based estimators optimized via GLASD under multiple contamination mechanisms, correlation structures, and dimensional regimes.

Data are generated under two correlation structures—Block-Toeplitz (sparse) and non-sparse—and four distributional settings: asymmetric Gaussian contamination (row-wise, column-wise, and random elementwise) and heavy-tailed sampling from a multivariate t-distribution. Performance is assessed across three dimensional regimes: \((p,n) = (20,100), (50,500), (100,500)\).

Please follow the steps below to reproduce the results reported in the manuscript.

1) **Data generation**
Run `Data_generation.m` for \((p,n) = (20,100), (50,500), (100,500)\), with - `structure_options = 2` (BlockToeplitz) and `structure_options = 4` (non-sparse), and - `dist_type = 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'`.

2) **Naive (Gaussian) estimator**
Run `Simulation_study_Naive.m` for \((p,n) = (20,100), (50,500), (100,500)\), with - `Ctype = 'BlockToeplitz'` and `'random'` (the `'random'` option corresponds to the non-sparse scenario), and - `dist_type = 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'`.

3) **Tyler’s M-estimator**
Run`Simulation_study_Tyler.m` for \((p,n) = (20,100), (50,500), (100,500)\), with - `Ctype = 'BlockToeplitz'` and `'random'`, and - `dist_type = 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'`.

4) **GLASD–Huber estimator**
Run`Simulation_study.m` for \((p,n) = (20,100), (50,500), (100,500)\), with - `Ctype = 'BlockToeplitz'` and `'random'`, and - `dist_type = 'gaussian_asymmetric_row', 'gaussian_asymmetric_col', 'gaussian_asymmetric_random', 't'`.

5) **Summarize simulation outputs**

After completing all simulation runs, execute `Simulation_output/Summary.m` and `Simulation_output/Summary_FINAL_TABLE.m` to generate the final summary table reported in the manuscript.

**Note.** The full simulation outputs are large in size (`Simulation_data/` contains more than 300 files and `Simulation_output/` contains more than 800 files) and are therefore not included in this repository. Only the code files required to fully reproduce the results are provided.
