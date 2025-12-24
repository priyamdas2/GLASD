# GLASD: Global Adaptive Stochastic Descent for Robust Correlation Estimation

This repository provides an implementation of **GLASD (Global Adaptive Stochastic Descent)**, a loss-agnostic, black-box optimization framework for **robust estimation of high-dimensional correlation matrices** under general Mahalanobis-type loss functions. The methods and experiments in this repository correspond to the manuscript:

> **Robust Correlation Estimation under General Loss Functions**  
> Priyam Das, Department of Biostatistics, Virginia Commonwealth University

GLASD is designed to reliably optimize **non-convex, non-smooth, and even discontinuous objectives** over the space of positive definite correlation matrices, without requiring gradient information or smoothness assumptions.

---

## Key Features

- **Loss-agnostic framework**  
  Supports Gaussian, Huber, truncated quadratic, Tukey’s biweight, and other user-defined Mahalanobis losses.

- **Geometry-aware optimization**  
  Uses a bijective hyperspherical (angular) parameterization of the correlation matrix space to enforce constraints by construction.

- **Derivative-free global optimization**  
  Combines adaptive stochastic coordinate descent with controlled forced exploration, enabling reliable global search.

- **Theoretical guarantees**  
  Includes global convergence guarantees and statistical consistency results under mild regularity conditions.

- **Practical applicability**  
  Demonstrated through extensive simulations and real-data applications, including cancer proteomic and microbiome datasets.

---

### Notation note

Throughout the code, comments, and documentation, the symbols **`M`** and **`p`** may occasionally be used interchangeably.  
Both symbols refer to the **same quantity**: the **dimension (number of rows/columns) of the correlation matrix of interest**.

---

## 📂 Repository Structure

```text

GLASD/
├── Benchmark GLASD Box/                      # Empirical convergence benchmarks for GLASD and competing optimizers on box-constrained test functions (Ackley, Griewank, Rastrigin, Rosenbrock)
├── Benchmark/                                # Benchmark comparison studies for robust correlation estimation methods reported in the manuscript
├── Case study high dim/                      # High-dimensional real-data case study with CRC microbiome analyses
├── Case study/                               # Real-data case study with cancer proteomics data
├── Simulation study factor based scenario/   # Simulation studies under approximate factor model settings 
├── Simulation study non-factor scenario/     # Simulation studies without factor structure 
├── images/                                   # Figures used in the manuscript and supplementary material
├── README.md                                 # Project overview and usage instructions


```




