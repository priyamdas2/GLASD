## Manopt Toolbox Requirement

For copyright reasons, the **Manopt** source files are not redistributed in this repository.  
Please download the Manopt MATLAB toolbox directly from the official website:

👉 https://www.manopt.org/

After downloading and extracting the toolbox, the `manopt` directory should have the following structure:

```text
manopt/
├── autodiff
├── core
├── lifts
├── manifolds
├── solvers
└── tools
```
After downloading, please replace the existing `manopt` folder in this repository (which contains only instructions and no source files) with the downloaded `manopt` folder.

---

## Code Reproducibility (for benchmark study)

To reproduce all numerical results reported in the main manuscript, please follow the steps below.

## Reproducing the Benchmark Results

This repository provides fully reproducible scripts to generate all benchmark results reported in the manuscript, including summary tables and boxplots.

### 1. Reproduce the Benchmark Summary Table

To reproduce the benchmark comparison table reported in the manuscript, proceed as follows.

**1.1 GLASD benchmarks**

Run `Benchmark_study_GLASD.m` for \( M \in \{5, 10, 20, 50\} \) and `input_vals = [1, 2, 3, 4]`. All output files are automatically saved in the `ALL_Outputs/` directory.

**1.2 Competing methods benchmarks**

Run  
`Benchmark_study_others.m` for \( M \in \{5, 10, 20, 50\} \) and `input_vals = [1, 2, 3, 4]`. All output files are automatically saved in the `ALL_Outputs/` directory.

**1.3 Summarize results**

After completing the above runs, execute  
`ALL_Outputs/Summary_ALL.m`, which aggregates all generated outputs and produces the benchmark summary table reported in the manuscript.

---

### 2. Reproduce the Benchmark Boxplots

To reproduce the boxplots based on repeated simulations, proceed as follows.

**2.1 GLASD (100 replications)**

Run `Benchmark_study_100_reps_GLASD.m` for \( M = 10 \) and `input_vals = [1, 2, 3, 4]`. All output files are automatically saved in the `ALL_Outputs/` directory.

**2.2 Competing methods (100 replications)**

Run `Benchmark_study_100_reps_others.m` for \( M = 10 \) and `input_vals = [1, 2, 3, 4]`. All output files are automatically saved in the `ALL_Outputs/` directory.

**2.3 Generate boxplots**

Finally, run `ALL_Outputs/Summary_ALL_100_reps_boxplot.R`, which constructs the benchmark boxplots from the aggregated outputs and reproduces the figures reported in the manuscript.

