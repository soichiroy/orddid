
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orddid

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/soichiroy/orddid.svg?branch=master)](https://travis-ci.org/soichiroy/orddid)
<!-- badges: end -->

  - Author: [Soichiro Yamauchi](https://soichiroy.github.io/)

  - For a detailed description of the method see:
    
    Difference-in-Differences for Ordinal Outcome

## Installation Instructions

You can install the development version from
[GitHub](https://github.com/) with:

``` r
require("devtools")
install_github("soichiroy/orddid", dependencies=TRUE)
```

## Example: Two Time Periods

This is a basic example which shows you how to solve a common problem:

``` r
## load package
library(orddid)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

## load example data
data("gun_twowave")

## Estimate causal effects
set.seed(1234)
fit <- ord_did(
  Ynew = gun_twowave %>% filter(year == 2012) %>% pull(guns),
  Yold = gun_twowave %>% filter(year == 2010) %>% pull(guns),
  treat = gun_twowave %>% filter(year == 2012) %>% pull(treat_100mi),
  id_cluster = gun_twowave %>% filter(year == 2010) %>% pull(reszip),
  cut = c(0, 1),
  n_boot = 100,
  pre = FALSE,
  verbose = FALSE
)

## view summary
summary(fit)
#> ── Effect Estimates ────────────────────────────────────────────────────
#>              Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> Delta[2-3] -0.00562 0.00485  -0.01314   0.00253  -0.01426    0.0035
#> Delta[3]    0.00431 0.00610  -0.00589   0.01302  -0.00679    0.0146

## non-cumulative effect
summary(fit, cumulative = FALSE)
#> ── Effect Estimates ────────────────────────────────────────────────────
#>            Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> Delta[1]  0.00562 0.00485  -0.00253   0.01314  -0.00350   0.01426
#> Delta[2] -0.00993 0.00775  -0.02247   0.00277  -0.02401   0.00422
#> Delta[3]  0.00431 0.00610  -0.00589   0.01302  -0.00679   0.01461
```

## Example: Additional Pre-treatment Period is Available

``` r
## load data
data("gun_threewave")
gun_threewave <- na.omit(gun_threewave)

## further subset to no-treated people through 2012
case_use <- gun_threewave %>%
  filter(year == 2012) %>%
  filter(pds_100mi == "Untreated in Previous Decade" & t_100mi == 0) %>%
  pull(caseid)
dat_14   <- gun_threewave %>% filter(caseid %in% case_use)

## check if subsetting is success full
## there should be no one treated until 2014
dat_14 %>% group_by(year, t_100mi) %>% summarize(n = n())
#> # A tibble: 4 x 3
#> # Groups:   year [3]
#>    year t_100mi     n
#>   <dbl>   <dbl> <int>
#> 1  2010       0  2811
#> 2  2012       0  2813
#> 3  2014       0  2143
#> 4  2014       1   664


## subset to comple-cases (exist from 2010 through 2014)
case14    <- dat_14 %>% filter(year == 2014) %>% pull(caseid)
case12    <- dat_14 %>% filter(year == 2012) %>% pull(caseid)
case10    <- dat_14 %>% filter(year == 2010) %>% pull(caseid)
case_full <- intersect(intersect(case14, case12), case10)

## treat Y2012 as "post" and Y2010 as "pre"
Ynew  <- dat_14 %>% filter(caseid %in% case_full & year == 2012) %>%
          pull(guns)
Yold  <- dat_14 %>% filter(caseid %in% case_full & year == 2010) %>%
          pull(guns)
treat <- dat_14 %>% filter(caseid %in% case_full & year == 2014) %>%
          pull(t_100mi)
zip   <- dat_14 %>% filter(caseid %in% case_full & year == 2014) %>%
          pull(reszip)

## estimate parameters
fit <- ord_did(Ynew, Yold, treat, zip, cut = c(0, 1),
               n_boot = 100, pre = TRUE, verbose = FALSE)

## equivalence test
equiv_test <- equivalence_test(
  object = fit, alpha = 0.05
)

## view result
summary(equiv_test)
#> ── Equivalence Test ────────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>        0.020332       -0.036161        0.021210        0.000204 
#> 
#> [1] H0 of no-equivalence is REJECTED with threshold 0.054

## plot result
plot(equiv_test, ylim = c(-0.1, 0.1), fill = FALSE)
```

![](man/figures/README-example2-1.png)<!-- -->

``` r

## test with different threshold
equiv_test2 <- equivalence_test(
  object = fit, alpha = 0.05, threshold = 0.01
)

## view result
summary(equiv_test2)
#> ── Equivalence Test ────────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>          0.0203         -0.0362          0.0212          0.8589 
#> 
#> [1] H0 of no-equivalence is NOT REJECTED with threshold 0.01
```
