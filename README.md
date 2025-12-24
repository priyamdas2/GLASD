# GLASD: Global Adaptive Stochastic Descent for Robust Correlation Estimation

This repository provides an implementation of **GLASD (Global Adaptive Stochastic Descent)**, a loss-agnostic, black-box optimization framework for **robust estimation of high-dimensional correlation matrices** under general Mahalanobis-type loss functions. The methods and experiments in this repository correspond to the manuscript:

> **Robust Correlation Estimation under General Loss Functions**  
> Priyam Das, Department of Biostatistics, Virginia Commonwealth University

GLASD is designed to reliably optimize **non-convex, non-smooth, and even discontinuous objectives** over the space of positive definite correlation matrices, without requiring gradient information or smoothness assumptions.

---

## 🔑 Key Features

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

## 📌 Problem Setting

GLASD is developed to address **robust estimation of correlation matrices** in high-dimensional settings where observations may exhibit outliers, heavy tails, or other departures from Gaussian assumptions. In such regimes, the sample correlation matrix is both statistically unstable and highly sensitive to extreme observations, motivating the use of alternative loss functions that control the influence of large Mahalanobis distances.

Rather than committing to a single likelihood or contamination model, GLASD formulates correlation estimation as a **general loss-based optimization problem** over the constrained space of positive definite correlation matrices with unit diagonal. Specifically, the target estimator is defined as the solution to

<pre>
minimize    f(C)
subject to  C ∈ 𝒞<sub>M</sub>
</pre>

where  
- \( 𝒞_M \) denotes the space of full-rank \( M \times M \) correlation matrices, and  
- \( f(C) \) is a user-specified Mahalanobis-type objective function.

This formulation accommodates a wide class of loss functions, including Gaussian likelihoods, Huber-type losses, truncated quadratic losses, and redescending M-estimators (e.g., Tukey’s biweight). Importantly, the resulting optimization problem is typically **non-convex, non-smooth, and sometimes discontinuous**, and is defined over a geometrically constrained parameter space.

GLASD is designed to solve this problem directly—without requiring gradient information, smoothness assumptions, or problem-specific algorithmic tailoring—by combining a geometry-aware reparameterization of the correlation matrix space with a derivative-free global optimization strategy.

---

## 🧠 Method Overview

GLASD operates by reparameterizing the space of correlation matrices using a **bijective hyperspherical (angular) mapping**, which transforms the constrained matrix optimization problem into a **box-constrained optimization problem in Euclidean space**. This mapping guarantees feasibility by construction: every point in the angular domain corresponds to a valid correlation matrix.

Once reformulated, optimization is performed using **Global Adaptive Stochastic Descent**, a derivative-free algorithm that combines adaptive coordinate-wise descent with controlled stochastic exploration. The method requires only objective evaluations and does not rely on gradients, smoothness, or continuity assumptions.

The overall geometry-aware optimization pipeline is illustrated in the following figures (available in the `images/` folder):

**Correlation matrix geometry and angular mapping** 
<p align="center">
  <img src="images/GLASD_mapping_v6.jpg" width="70%">
</p


**GLASD algorithmic flowchart**
<p align="center">
  <img src="images/GLASD_flowchart_v8.jpg" width="70%">
</p

Together, these figures summarize how GLASD enforces correlation constraints, adapts step sizes and directional probabilities, and balances exploitation with global exploration.

---

## 📐 Theoretical Guarantees

GLASD is equipped with formal **global convergence guarantees** under mild regularity conditions on the objective function, requiring only boundedness and measurability. In addition, when applied to loss-based correlation estimation, the resulting estimators are shown to be **statistically consistent** under a uniform law of large numbers. Importantly, these results hold without assuming convexity, differentiability, or smoothness of the loss function, making GLASD applicable to truncated and redescending objectives that fall outside the scope of classical optimization theory.

---

## 📊 Benchmark Experiments

To assess optimization performance on structured constrained domains, GLASD is benchmarked against classical local and global optimization methods using four standard non-convex test functions:

- Ackley  
- Griewank  
- Rastrigin  
- Rosenbrock  

These functions are embedded into the correlation matrix space by evaluating them on the off-diagonal entries of valid correlation matrices. Summary comparisons and variability across dimensions are reported in the manuscript and supplementary material.

Representative benchmark visualizations are provided as boxplots below, summarizing optimization performance.

<p align="center">
  <img src="images/Boxplots_Fun_and_Time_all_functions.png" width="80%">
</p
---

## 🧪 Simulation Studies

GLASD is evaluated through two complementary simulation regimes:

- **Non-factor-based scenarios**  
  Dense correlation structures under heavy-tailed and contaminated distributions, assessing robustness without latent structure assumptions.

- **Factor-based scenarios**  
  Approximate factor models with structured idiosyncratic covariance, comparing GLASD to POET-based and robust factor methods under known and unknown factor settings.

These studies demonstrate that GLASD produces stable and competitive correlation estimates across a wide range of dimensions, contamination mechanisms, and structural assumptions.

---

## 🧬 Case Study: CRC Microbiome Correlation Analysis

GLASD is applied to a colorectal cancer (CRC) microbiome dataset to estimate robust correlation structures among microbial taxa. Using a Huber-loss-based objective, GLASD yields a **spectrally stable and interpretable correlation estimate**, in contrast to classical estimators that exhibit numerical instability under heavy-tailed behavior.

A representative heatmap from this analysis is available in the `images/` directory and is reported in the manuscript and supplementary material.

---

This repository is intended to serve both as a **research artifact accompanying the manuscript** and as a **general-purpose implementation** for robust, geometry-aware correlation matrix optimization.








