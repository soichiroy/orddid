
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orddid

<!-- badges: start -->

<!-- badges: end -->

  - Author: [Soichiro Yamauchi](https://soichiroy.github.io/)

  - For a detailed description of the method see:
    
    Difference-in-Differences for Ordinal Outcome

## Installation Instructions

You can install the released version of orddid from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("orddid")
```

And the development version from [GitHub](https://github.com/) with:

``` r
require("devtools")
install_github("soichiroy/orddid", dependencies=TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## load package 
library(orddid)

## load example data 
data(gun_twowave)


## Estimate causal effects 
set.seed(1234)
fit <- ord_did(
  Ynew = gun_twowave$guns12,
  Yold = gun_twowave$guns10,
  treat = gun_twowave$treat_100mi,
  cut = c(0, 1),
  n_boot = 500,
  pre = FALSE,
  verbose = TRUE
)
#> 
 50 out of 500 bootstrap iterations
 100 out of 500 bootstrap iterations
 150 out of 500 bootstrap iterations
 200 out of 500 bootstrap iterations
 250 out of 500 bootstrap iterations
 300 out of 500 bootstrap iterations
 350 out of 500 bootstrap iterations
 400 out of 500 bootstrap iterations
 450 out of 500 bootstrap iterations
 500 out of 500 bootstrap iterations

## view summary 
summary(fit)
#> ── Effect Estimates ────────────────────────────────────
#>              Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> Delta[2-3] -0.01385 0.00724  -0.02619  -0.00153   -0.0285  -0.00039
#> Delta[3]    0.00693 0.00889  -0.00816   0.02076   -0.0114   0.02481
```
