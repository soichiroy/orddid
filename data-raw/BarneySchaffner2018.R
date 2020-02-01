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
require(labelled)

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
  rename(caseid = caseid_1) %>%
  select(caseid, guns12, guns10, treat_100mi, treat10_12_2, partyid3_10, pds_100mi, zip) %>%
  na.omit() -> gun_twowave_sub

## save object
usethis::use_data(gun_twowave_sub, overwrite = TRUE)


## ------------------------------------------------------------------------- ##
##                  create two wave data: "gun_twowave_full"                 ##
## ------------------------------------------------------------------------- ##
## there are some observations whos `regzip_post` differ across waves,
## even though `no_move == 1`
##

## load data
dat1 <- read_dta("data-raw/final_longform_10_12_merged.dta")

## some zip changes even after taking `no_move == 1`
## keep only ones no changes in zip
dat1 %>%
  filter(no_move == 1) %>%
  pivot_wider(id_cols = caseid,
    names_from = year, names_prefix = "zip",
    values_from = regzip_post
  ) %>%
  filter(zip2010 != zip2012) %>%
  pull(caseid) -> remove_id

## create a group (treatment / control)
# dat1 %>%
#   filter(no_move == 1) %>%
#   pivot_wider(id_cols = caseid,
#     names_from = year, names_prefix = "treat",
#     values_from = treat_100mi
#   ) %>%
#   mutate(treat_100mi = pmax(treat2010, treat2012)) %>%
#   select(caseid, treat_100mi) -> treat_100mi
#
# dat1 %>%
#   filter(no_move == 1) %>%
#   pivot_wider(id_cols = caseid,
#     names_from = year, names_prefix = "treat",
#     values_from = treat_25mi
#   ) %>%
#   mutate(treat_25mi = pmax(treat2010, treat2012)) %>%
#   select(caseid, treat_25mi) -> treat_25mi

## final data
gun_twowave <- dat1 %>%
  filter(no_move == 1) %>%
  filter(!caseid %in% remove_id) %>%
  mutate(guns = ifelse(CC320 == 2, 1, ifelse(CC320 == 3, 2, 3))) %>%
  select(caseid, year, guns, pds_100mi, pds_25mi, regzip_post,
    party2010, pid3, pid7,
    treat_25mi, treat_100mi) %>%
  rename(reszip = regzip_post) %>%
  add_value_labels(
    treat_100mi = c(Untreated = 0, Treated = 1),
    treat_25mi = c(Untreated = 0, Treated = 1),
    guns = c("Less Strict" = 1, "Kept As They Are" = 2, "More Strict" = 3)
  ) %>%
  na.omit()

## keep only the balanced observations
gun_twowave %>%
  pivot_wider(id_cols = c(caseid),
    names_from = year, names_prefix = 'year',
    values_from = guns
) %>%
  na.omit() %>%
  pull(caseid) -> id_keep


## save to R
gun_twowave <- gun_twowave %>%
  filter(caseid %in% id_keep)
usethis::use_data(gun_twowave, overwrite = TRUE)


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
gun_threewave <- dat2 %>%
  filter(no_move_allwaves==1) %>%
  select(caseid, year, guns, t_100mi, t_25mi, pds_100mi, pds_25mi, reszip,
    pid3, pid7, party2010)

## save object
usethis::use_data(gun_threewave, overwrite = TRUE)
