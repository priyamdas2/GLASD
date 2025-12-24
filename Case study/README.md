## Case Study: Cancer Proteomic Data Correlation Analysis

This case study illustrates robust correlation estimation for high-dimensional cancer proteomic data using GLASD with a Huber loss, highlighting stable pathway-level association patterns in the presence of outliers. 

To reproduce the results reported in the manuscript and supplementary material, follow the steps below.

1) **Data extraction**
Run `BRCA_data_extraction.R` to extract and preprocess the breast cancer proteomic dataset used in this analysis.

2) **Correlation estimation**
Run `BRCA_analysis.m` to estimate correlation matrices by optimizing the Huber-loss-based objective using the GLASD algorithm.

3) **Visualization**
Run `BRCA_plot.R` to generate the heatmaps and summary plots reported in this section of the manuscript and supplementary material.
