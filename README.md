
<!-- badges: start -->

![status](https://img.shields.io/badge/status-under%20construction-yellow)
<!-- badges: end -->

## Introduction

The `ReplicateMST2018` package…

## Installation

The `ReplicateMST2018` package is hosted on GitHub at
<https://github.com/ian-xu-economics/ReplicateMST2018/>. It can be
installed using the `remotes::install_github()` function:

``` r
# install.packages("remotes")
remotes::install_github("ian-xu-economics/ReplicateMST2018")
```

## Using `ReplicateMST2018`

After installing `ReplicateMST2018`, we can attach the package to our
session using the base `library()` function. We’ll also attach some
other packages needed to run the examples in the vignette.

``` r
library(ReplicateMST2018)
library(tidyverse)
library(glue)
library(latex2exp)
```

## Replicating Figures

### The DGP

In order to produce the figures in Mogstad, Santos, and Torgovitsky
(2018), we first need to duplicate the data generating process (DGP)
used in MST 2018. They consider a simple example with a trinary
instrument, $Z \in \{0,1,2\}$, with $P(Z = 0) = 0.5$, $P(Z = 1) = 0.4$,
and $P(Z = 2) = 0.1$. The propensity score is specified as
$p(0) = 0.35$, $p(1) = 0.6$, and $p(2) = 0.7$. They take the outcome
$Y \in \{0,1\}$ to be binary and restrict $\mathcal{M}$ to contain only
MTR pairs that are bounded between 0 and 1. The data are generated using
the MTR functions

$$
m_0(u) = 0.6b^2_0(u) + 0.4b^2_1(u) + 0.3b^2_2(u) \quad \text{ and}
$$ $$
m_1(u) = 0.75b^2_0(u) + 0.5b^2_1(u) + 0.25b^2_2(u)\text{,}
$$ where $b^2_k$ is the $k$th Bernstein basis polynoial of degree 2.

This DGP has already been pre-programmed into `dgp_MST2018()`, so we’ll
just save the DGP used in the paper into `dgp`.

``` r
dgp <- dgp_MST2018()
```

Because of the complexity of the code used to produce the figures, the
code is quite long. Therefore, we do not include the code here; we only
output the figures. To see the code, open the `README.Rmd` file.

### Figure 1

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Figure 2

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

### Figure 3

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Figure 4

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Figure 5

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Figure 6

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### Figure 7

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Figure 8

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
