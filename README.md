
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orddid

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/soichiroy/orddid.svg?branch=master)](https://travis-ci.org/soichiroy/orddid)
<!-- badges: end -->

  - Author: [Soichiro Yamauchi](https://soichiroy.github.io/)

## Installation Instructions

You can install the development version from
[GitHub](https://github.com/) with:

``` r
require("devtools")
install_github("soichiroy/orddid", dependencies=TRUE)
```

## Example: Two Time Periods


``` r
## Estimate causal effects
set.seed(1234)
fit <- ord_did(
  Ynew = gun_twowave %>% filter(year == 2012) %>% pull(guns),
  Yold = gun_twowave %>% filter(year == 2010) %>% pull(guns),
  treat = gun_twowave %>% filter(year == 2012) %>% pull(treat_100mi),
  id_cluster = gun_twowave %>% filter(year == 2010) %>% pull(reszip),
  n_boot = 100,
  pre = FALSE,
  verbose = FALSE
)

## view summary of the output 
## non-cumulative effects
summary(fit, cumulative = FALSE)
#> ── Effect Estimates ──────────────────────────────────────────────────
#>           Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> zeta[1]  0.00559 0.00476  -0.00237   0.01233  -0.00544   0.01365
#> zeta[2] -0.00991 0.00716  -0.02138   0.00119  -0.02282   0.00404
#> zeta[3]  0.00432 0.00634  -0.00590   0.01536  -0.00882   0.01647

## cumulative effects
summary(fit)
#> ── Effect Estimates ──────────────────────────────────────────────────
#>              Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> Delta[2-3] -0.00559 0.00476   -0.0123   0.00237  -0.01365   0.00544
#> Delta[3]    0.00432 0.00634   -0.0059   0.01536  -0.00882   0.01647
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
fit <- ord_did(Ynew, Yold, treat, zip,
               n_boot = 100, pre = TRUE, verbose = FALSE)

## equivalence test
equiv_test <- equivalence_test(
  object = fit, alpha = 0.05
)

## view result
summary(equiv_test)
#> ── Equivalence Test ──────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>         0.02033        -0.03847         0.02321         0.00102 
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
#> ── Equivalence Test ──────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>          0.0203         -0.0385          0.0232          0.8267 
#> 
#> [1] H0 of no-equivalence is NOT REJECTED with threshold 0.01
```
