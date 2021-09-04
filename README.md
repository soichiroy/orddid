
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orddid

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/soichiroy/orddid)
[![Codecov test
coverage](https://codecov.io/gh/soichiroy/orddid/branch/master/graph/badge.svg)](https://codecov.io/gh/soichiroy/orddid?branch=master)
[![R build
status](https://github.com/soichiroy/orddid/workflows/R-CMD-check/badge.svg)](https://github.com/soichiroy/orddid/actions)
<!-- badges: end -->

  - Author: [Soichiro Yamauchi](https://soichiroy.github.io/)
  - Package manual:
    [PDF](https://soichiroy.github.io/files/packages/orddid_0_1_0.pdf)

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
#> ── Effect Estimates ────────────────────────────────────────────────────────────
#>           Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> zeta[1]  0.00559 0.00485  -0.00252   0.01314  -0.00347   0.01426
#> zeta[2] -0.00991 0.00775  -0.02245   0.00276  -0.02400   0.00421
#> zeta[3]  0.00432 0.00611  -0.00587   0.01301  -0.00687   0.01459

## cumulative effects
summary(fit)
#> ── Effect Estimates ────────────────────────────────────────────────────────────
#>              Effect      SE 90% Lower 90% Upper 95% Lower 95% Upper
#> Delta[2-3] -0.00559 0.00485  -0.01314   0.00252  -0.01426   0.00347
#> Delta[3]    0.00432 0.00611  -0.00587   0.01301  -0.00687   0.01459
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
#> `summarise()` has grouped output by 'year'. You can override using the `.groups` argument.
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

## wald test 
wald_test <- wald_test(fit) 
round(wald_test$p_value, 3)
#>      [,1]
#> [1,] 0.28


pnorm_00 <- function(x) pnorm(x, mean = fit$fit$theta$mu00, sd = fit$fit$theta$sd00)
pnorm_10 <- function(x) pnorm(x, mean = fit$fit$theta$mu10, sd = fit$fit$theta$sd10)
pnorm_01 <- function(x) pnorm(x, mean = fit$fit$theta$mu01, sd = fit$fit$theta$sd01)
pnorm_11 <- function(x) pnorm(x, mean = fit$fit$theta$mu11, sd = fit$fit$theta$sd11)

# par(mfrow = c(1, 2)) 
curve(pnorm_00, from = -3, to = 4, xlab = "", ylab = "")
curve(pnorm_10, from = -3, to = 4, col = 'red', add = TRUE)
curve(pnorm_01, from = -3, to = 4, add = TRUE, lwd = 1.5, lty = 2)
curve(pnorm_11, from = -3, to = 4, col = 'red', add = TRUE, lwd = 1.5, lty = 2)
```

![](man/figures/README-example2-1.png)<!-- -->

``` r

## equivalence test
equiv_test <- equivalence_test(
  object = fit, alpha = 0.05
)

## view result
summary(equiv_test)
#> ── Equivalence Test ────────────────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>         0.02033        -0.04272         0.03611         0.00632 
#> 
#> [1] H0 of non-equivalence is REJECTED with threshold 0.054

## plot result
plot(equiv_test, ylim = c(-0.1, 0.1), fill = FALSE)
```

![](man/figures/README-example2-2.png)<!-- -->

``` r

## test with different threshold
equiv_test2 <- equivalence_test(
  object = fit, alpha = 0.05, threshold = 0.01
)

## view result
summary(equiv_test2)
#> ── Equivalence Test ────────────────────────────────────────────────────────────
#> Estimate (tmax)           Lower           Upper          pvalue 
#>          0.0203         -0.0427          0.0361          0.7835 
#> 
#> [1] H0 of non-equivalence is NOT REJECTED with threshold 0.01
```
