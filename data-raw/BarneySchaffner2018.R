#
#
# title: R script to generate "gun_twowave" and "gun_threewave"
#
# author: Soichiro Yamauchi
#
# note:
#   Original data can be found at
#     Barney, David; Schaffner, Brian, 2018, "Replication Data for:
#        Reexamining the Effect of Mass Shootings on Public Support for Gun Control",
#     https://doi.org/10.7910/DVN/YJQIXP, Harvard Dataverse
#
#

##
## load packages
##
require(dplyr)
require(ggplot2)
require(haven)
require(readr)
require(tidyr)

## ------------------------------------------------------------------------- ##
##                  create two wave data: "gun_twowave_sub"                  ##
##                  two year subset of 2012-2014 panel                       ##
## ------------------------------------------------------------------------- ##
##

## load data
dat <- read_dta("data-raw/CCES_10_12_14_panel_2year_subset.dta")

## subset variables
dat %>%
  filter(sameres == 1) %>% ## only no-movers
  select(guns12, guns10, treat_100mi, treat10_12_2, partyid3_10, pds_100mi, zip) %>%
  na.omit() -> gun_twowave_sub

## save object
usethis::use_data(gun_twowave_sub, overwrite = TRUE)


## ------------------------------------------------------------------------- ##
##                  create two wave data: "gun_twowave_full"                  ##
## ------------------------------------------------------------------------- ##
## load data
dat1 <- read_dta("data-raw/final_longform_10_12_merged.dta")


dat1 %>% filter(year == 2010) %>%
  pull(CC320) %>%
  summary(.)
dat1 %>% filter(year == 2012) %>%
  pull(CC320) %>%
  summary(.)

table(
  dat1 %>% filter(year == 2010) %>%
    pull(CC320),
  dat1 %>% filter(year == 2012) %>%
    pull(CC320)
)


dat1 %>%
  filter(no_move == 1) %>%
  mutate(rezip = regzip_post) %>%
  select(caseid, year, CC320, pds_100mi, pds_25mi) %>%
  pivot_wider(
    id_cols = c(caseid, pds_100mi, pds_25mi),
    names_from = year, names_prefix = "guns",
    values_from = CC320) -> gun_twowave_full

  summary(gun_twowave_full)


    , treat_25mi, treat_100mi, pds_100mi, pds_25mi,
      regzip_post) %>%

## ------------------------------------------------------------------------- ##
##                create three wave data: "gun_threewave"                    ##
## ------------------------------------------------------------------------- ##

## The original data is in .dta format. However, I was not able to run the
##    the following R code:
##        dat <- haven::read_dta("final_longform_10_14_merged.dta")
##    which generates the following error message:
##        Error message: Failed to parse final_longform_10_14_merged.dta:
##        Unable to convert string to the requested encoding (invalid byte sequence).
##
## Therefore, I manually opned the filed in State and export the file to .csv format
##


## load data
dat2 <- read_csv("data-raw/final_longform_10_14_merged.csv")

## subset variables
dat2 %>%
  filter(no_move_allwaves==1) %>%
  select(caseid, year, guns, t_100mi, t_25mi, pds_100mi, pds_25mi, reszip) %>%
  na.omit() -> gun_threewave

## save object
usethis::use_data(gun_threewave, overwrite = TRUE)
