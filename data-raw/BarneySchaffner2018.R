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

## ------------------------------------------------------------------------- ##
##                  create two wave data: "gun_twowave"                      ## 
## ------------------------------------------------------------------------- ##

## load data
dat <- read_dta("CCES_10_12_14_panel_2year_subset.dta")

## subset variables 
dat %>%
  filter(sameres == 1) %>% ## only no-movers 
  select(guns12, guns10, treat_100mi, treat10_12_2, partyid3_10, pds_100mi) %>%
  na.omit() -> gun_twowave

## save object
usethis::use_data(gun_twowave)


## ------------------------------------------------------------------------- ##
##                create three wave data: "gun_threewave"                    ## 
## ------------------------------------------------------------------------- ##

## The original data is in .dta format. However, I was not able to run the 
##    the following R code:
##        dat <- read_dta("final_longform_10_14_merged.dta")
##    which generates the following error message:
##        Error message: Failed to parse final_longform_10_14_merged.dta: 
##        Unable to convert string to the requested encoding (invalid byte sequence).
##
## Therefore, I manually opned the filed in State and export the file to .csv format 
##


## load data 
dat2 <- read_csv("final_longform_10_14_merged.csv")

## subset variables 
dat2 %>% 
  filter(no_move_allwaves==1) %>%
  select(caseid, year, guns, t_100mi, t_25mi, pds_100mi, pds_25mi) %>%
  na.omit() -> gun_threewave

## save object
usethis::use_data(gun_threewave)
