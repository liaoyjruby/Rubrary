
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rubrary

<!-- badges: start -->
<!-- badges: end -->

Personal functions, mostly for Graeber Lab analysis and to test R
package building. All functions subject to sudden changes.

Installing this package will install many packages used throughout
Graeber lab analyses.

## Installation

You can install the development version of Rubrary from
[GitHub](https://github.com/) with:

``` r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("liaoyjruby/Rubrary")
```

## Features

### Exploratory Data Analysis

- `check_normal`
- `plot_distribution`
- `plot_density`
- `plot_scatter`
- `plot_scatter_compare`
- `plot_boxplot_paired`

### Principal Components Analysis

- `run_PCA`
  - `plot_screeplot`
- `plot_PCA`

### Differential Gene / Pathway Analysis

- `run_RRHO`

### Survival Analysis
