## load packages
require(dplyr)
require(ggplot2)
require(haven)

## load data
dat <- read_dta("CCES_10_12_14_panel_2year_subset.dta")

## subset data
dat_use <- dat %>%
  filter(sameres == 1) %>%
  select(guns12, guns10, treat_100mi, treat10_12_2, partyid3_10, pds_100mi) %>%
  na.omit()

## save object
barney_schaffner_18 <- dat_use
usethis::use_data(barney_schaffner_18)
