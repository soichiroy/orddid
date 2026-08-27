
<!-- README.md is generated from README.Rmd. Please edit that file -->

# orddid

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-1.2.0-blue.svg)](https://github.com/soichiroy/orddid)
[![Codecov test
coverage](https://codecov.io/gh/soichiroy/orddid/branch/master/graph/badge.svg)](https://codecov.io/gh/soichiroy/orddid?branch=master)
[![R build
status](https://github.com/soichiroy/orddid/workflows/R-CMD-check/badge.svg)](https://github.com/soichiroy/orddid/actions)
<!-- badges: end -->

R package for difference-in-differences with ordinal outcomes (Yamauchi
2026+).

- Author: [Soichiro Yamauchi](https://soichiroy.github.io/)

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
require("devtools")
install_github("soichiroy/orddid", dependencies = TRUE)
```

## Example: Parametric Estimation

The `gun_twowave` dataset contains a two-wave panel (2010 and 2012) from
the CCES, merged with mass-shooting treatment data from Barney and
Schaffner (2019). The outcome `guns` is a three-category ordinal measure
of attitudes toward gun control.

``` r
library(orddid)
library(dplyr)

data("gun_twowave")
```

The `ord_did()` function takes a data frame with column names for the
outcome, a post-period indicator, a treatment indicator, and a cluster
variable. Because treatment assignment is determined by whether a mass
shooting occurred nearby between 2010 and 2012, we construct a
time-invariant treatment indicator from the post-period observations.

``` r
## identify treated individuals (ever treated in post period)
treated_ids <- gun_twowave %>%
  filter(treat_100mi == 1) %>%
  distinct(caseid)

dat <- gun_twowave %>%
  mutate(
    post  = (year == 2012),
    treat = caseid %in% treated_ids$caseid,
    guns  = as.numeric(guns)
  )

## estimate effects
set.seed(1234)
fit <- ord_did(
  data    = dat,
  outcome = "guns",
  post    = "post",
  treat   = "treat",
  cluster = "reszip",
  n_boot  = 50
)

## category-specific effects
fit$estimate_effects
#> # A tibble: 3 × 4
#>   category   effect lower_ci upper_ci
#>      <int>    <dbl>    <dbl>    <dbl>
#> 1        1  0.00562 -0.00264  0.0129 
#> 2        2 -0.00993 -0.0246   0.00341
#> 3        3  0.00431 -0.00650  0.0150

## bounds on the relative treatment effect
fit$relative_effects
#> # A tibble: 2 × 8
#>   tau_assumption       effect_lb effect_ub status lower_ci upper_ci c_crit
#>   <chr>                    <dbl>     <dbl> <chr>     <dbl>    <dbl>  <dbl>
#> 1 none                   -0.155     0.158  ok      -0.168    0.173    1.64
#> 2 stochastic_dominance   -0.0608    0.0601 ok      -0.0714   0.0708   1.64
#> # ℹ 1 more variable: n_boot_valid <int>
```

## Example: Semiparametric Estimation with Covariates

The semiparametric estimator adjusts for covariates using an NPMLE-based
approach. The `gun_twowave_cov` dataset is a smaller subset of the panel
augmented with IRT ideology scores and demographics.

The covariate selected by `fixed_index` (the first by default) must be
continuous, nonconstant, and oriented so its coefficient in the
package’s internal index is positive. Here `ideology` (an IRT score)
satisfies this requirement. The default semiparametric regularization is
an I-spline; the ridge version is selected with
`semiparametric_regularization = "ridge"`. Both use the paper’s two
alternating first-stage updates followed by a final NPMLE refit;
fixed-point and optimizer diagnostics are returned in
`point_diagnostics`.

The three estimators used in the paper’s simulation are available
through the same interface: use `method = "parametric_covariates"` for
probit, `method = "semiparametric"` for NPMLE plus I-spline, and add
`semiparametric_regularization = "ridge"` for NPMLE plus ridge. For the
Appendix design `covariates = c("intercept", "X1", "X2")`, use
`fixed_index = 2`; the intercept is held fixed and the continuous `X1`
coefficient is normalized.

``` r
data("gun_twowave_cov")

treated_ids_cov <- gun_twowave_cov %>%
  filter(treat_100mi == 1) %>%
  distinct(caseid)

dat_cov <- gun_twowave_cov %>%
  mutate(
    post  = (year == 2012),
    treat = caseid %in% treated_ids_cov$caseid
  )

set.seed(1234)
fit_semi <- ord_did(
  data       = dat_cov,
  outcome    = "guns",
  post       = "post",
  treat      = "treat",
  cluster    = "reszip",
  n_boot     = 0, # keep the rendered example fast; use 500+ for inference
  method     = "semiparametric",
  covariates = c("ideology", "ideology2", "female", "educ")
)

## category-specific effects
fit_semi$estimate_effects
#> # A tibble: 3 × 4
#>   category   effect lower_ci upper_ci
#>      <int>    <dbl>    <dbl>    <dbl>
#> 1        1 -0.0612        NA       NA
#> 2        2 -0.00138       NA       NA
#> 3        3  0.0626        NA       NA

## bounds on the relative treatment effect
fit_semi$relative_effects
#> # A tibble: 2 × 8
#>   tau_assumption       effect_lb effect_ub status lower_ci upper_ci c_crit
#>   <chr>                    <dbl>     <dbl> <chr>     <dbl>    <dbl>  <dbl>
#> 1 none                   -0.0944     0.281 ok           NA       NA     NA
#> 2 stochastic_dominance    0.0362     0.171 ok           NA       NA     NA
#> # ℹ 1 more variable: n_boot_valid <int>
```

## Example: Equivalence Test

When an additional pre-treatment period is available, the distributional
parallel trends assumption can be assessed using `run_equivalence()`.
The `gun_threewave` dataset contains a three-wave panel (2010, 2012,
2014). We use the 2010–2012 waves to test the assumption among
individuals who remain untreated through 2012.

``` r
data("gun_threewave")
gun_threewave <- na.omit(gun_threewave)

## keep individuals untreated through 2012
case_use <- gun_threewave %>%
  filter(year == 2012,
         pds_100mi == "Untreated in Previous Decade",
         t_100mi == 0) %>%
  pull(caseid)

## restrict to complete 2010-2012 cases
case10 <- gun_threewave %>%
  filter(caseid %in% case_use, year == 2010) %>% pull(caseid)
case12 <- gun_threewave %>%
  filter(caseid %in% case_use, year == 2012) %>% pull(caseid)
case_full <- intersect(case10, case12)

## individuals treated between 2012 and 2014
future_treated <- gun_threewave %>%
  filter(year == 2014, t_100mi == 1) %>%
  pull(caseid)

dat_eq <- gun_threewave %>%
  filter(caseid %in% case_full, year %in% c(2010, 2012)) %>%
  mutate(
    post  = (year == 2012),
    treat = caseid %in% future_treated,
    guns  = as.numeric(guns)
  )
```

``` r
set.seed(1234)
eq_result <- run_equivalence(
  data    = dat_eq,
  outcome = "guns",
  post    = "post",
  treat   = "treat",
  cluster = "reszip",
  n_boot  = 50
)

## worst-case bias inflation factor (1/M)
eq_result$inf_factor
#> [1] 1.135909

## plot the equivalence test result
plot_equivalence(eq_result)
```

![](man/figures/README-equiv-fit-1.png)<!-- -->

## References

Yamauchi, Soichiro (2026+). “[Difference-in-Differences for Ordinal
Outcomes: Application to the Effect of Mass Shootings on Attitudes
toward Gun
Control](https://soichiroy.github.io/files/papers/did_ordinal_paper2.pdf).”
