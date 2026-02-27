## ------------------------------------------------------------------
## Create gun_twowave_cov: two-wave panel with individual covariates
## ------------------------------------------------------------------
##
## Source:
##   gun_twowave_with_covariates.rds produced by 02-construct-covariates.R
##   in the paper replication files.  That script:
##     1. Loads the CCES 2010-12 panel and the mass-shooting treatment data
##     2. Fits a 2PL IRT model to CCES policy items and extracts ideology
##        scores (ideology, ideology2 = second dimension)
##     3. Extracts demographics (female, educ, race) from raw CCES
##
## This script subsamples ~2000 individuals (both waves, ~4000 rows)
## to create a small dataset suitable for package examples.
## It is included for reproducibility reference only and is not meant
## to be re-run during package installation.
## ------------------------------------------------------------------

library(dplyr)

## -- load augmented data ------------------------------------------------
dat_full <- readRDS(
  "path/to/gun_twowave_with_covariates.rds"
)

## -- keep relevant columns ----------------------------------------------
cols_keep <- c(
  "caseid", "year", "guns",
  "treat_25mi", "treat_100mi", "reszip",
  "ideology", "ideology2", "female", "educ"
)
dat_full <- dat_full[, cols_keep]

## -- subsample individuals ----------------------------------------------
set.seed(20250220)
ids_all <- unique(dat_full$caseid)
ids_use <- sample(ids_all, size = 2000)

gun_twowave_cov <- dat_full %>%
  filter(caseid %in% ids_use) %>%
  as.data.frame()

## -- basic checks -------------------------------------------------------
stopifnot(
  length(unique(gun_twowave_cov$caseid)) == 2000,
  all(c("ideology", "female", "educ") %in% names(gun_twowave_cov))
)

## -- save ---------------------------------------------------------------
save(gun_twowave_cov, file = "data/gun_twowave_cov.rda", compress = "xz")
