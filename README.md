
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epievolve

<!-- badges: start -->

<!-- badges: end -->

The goal of epievolve is to estimate time-varying Rt(t) using the
odin-monty framework. Estimates Rt(t) from time series of case counts
using a mechanistic SIR model. Enables interpretable inference on
epidemic dynamics and intervention effects.

## Installation

You can install the development version of epievolve from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("meaganjharrington/epievolve")
```

## Example

``` r
# not sure if this works as of 7/01
library(epievolve)
df <- read.csv(incidence.csv)

cases <- c(3,2,2,2,1,3)
time <- c(1,2,3,4,5,6)

out <- final_estimate_Rt_step(
incidence = df,
N = 1e6,
gamma = 1/6,
beta_breaks = c(1, 17, 35)
)

plot(out$Rt_series)
```
