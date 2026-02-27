#' Two-wave panel from Barney and Schaffner (2019)
#' 
#' \code{gun_twowave} is a part of the replication data from Barney and Schaffner (2019).
#' The data set is constructed combining the two-wave panel from the Cooperative Congressional Election Study (CCES) 2010-12 study 
#'  with the mass shooting data originally analyzed by Newman and Hartman (2019).
#' 
#' @name gun_twowave
#' @docType data 
#' @format An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 33106 rows and 11 columns.
#'   \describe{
#'      \item{caseid}{Unique id of a respondent.}
#'      \item{year}{Year of the response.}
#'      \item{guns}{Attitudes towards gun control, corresponding to a variable `CC320' in CCES. 
#'                  (1) "Less Strict", (2) "Kept As They Are" and (3) "More Strict".}
#'      \item{pds_25mi}{``Prior exposure'' variable, taking (1) if the mass shootings happened within 25 miles between 2000 and 2010.}
#'      \item{pds_100mi}{``Prior exposure'' variable, taking (1) if the mass shootings happened within 100 miles between 2000 and 2010.}
#'      \item{reszip}{Zipcode of a respondent.}
#'      \item{party2010}{3-point scale party id constructed by Barney and Schaffner (2019) based on `pid7' in 2010.}
#'      \item{pid3}{3-point scale party id, measured in 2010.}
#'      \item{pid7}{7-point scale party id, measured in 2010.}
#'      \item{treat_25mi}{Treatment indicator which takes 1 if the mass shootings happened within 25 miles.}
#'      \item{treat_100mi}{Treatment indicator which takes 1 if the mass shootings happened within 100 miles.}
#'    }
#' @source 
#' Barney, David; Schaffner, Brian. 
#'  ``Replication Data for: Reexamining the Effect of Mass Shootings on Public Support for Gun Control.''
#'  \url{https://doi.org/10.7910/DVN/YJQIXP}, Harvard Dataverse, V1, 2018.
#'
#'  Ansolabehere, Stephen; Schaffner, Brian. 
#'   ``2010 - 2012 CCES Panel Study.''
#'   \url{https://doi.org/10.7910/DVN/24416}, Harvard Dataverse, V4, 2014.
#'
#' @references
#' Barney, David J., and Brian F. Schaffner. 
#'  "Reexamining the Effect of Mass Shootings on Public Support for Gun Control." 
#'   British Journal of Political Science 49.4 (2019): 1555-1565. 
#'   \url{https://doi.org/10.1017/S0007123418000352}.
#'
#' Newman, B. J. and Hartman, T. K.
#'  ``Mass shootings and public support for gun control.''
#'  British Journal of Political Science, 49.4 (2019): 1527–1553.
#'  \url{https://doi.org/10.1017/S0007123417000333}.
#' @keywords dataset 
#' @examples
#'  require(dplyr)
#'  require(haven)
#'  data("gun_twowave")
#'  gun_twowave
"gun_twowave"


#' Three-wave panel from Barney and Schaffner (2019)
#'
#' \code{gun_threewave} is a part of the replication data from Barney and Schaffner (2019).
#' The data set is constructed combining the three-wave panel from the Cooperative Congressional Election Study (CCES) 2010-12-14 study 
#'  with the mass shooting data originally analyzed by Newman and Hartman (2019).
#'
#' @name gun_threewave
#' @docType data
#' @format An object of class spec_tbl_df (inherits from tbl_df, tbl, data.frame) with 23832 rows and 11 columns.
#' @source 
#' Barney, David; Schaffner, Brian. 
#'  ``Replication Data for: Reexamining the Effect of Mass Shootings on Public Support for Gun Control.''
#'  \url{https://doi.org/10.7910/DVN/YJQIXP}, Harvard Dataverse, V1, 2018.
#'
#' Schaffner, Brian; Ansolabehere, Stephen.
#'  ``2010-2014 Cooperative Congressional Election Study Panel Survey.''
#'  \url{https://doi.org/10.7910/DVN/TOE8I1}, Harvard Dataverse, V11, 2015.
#'
#' @references
#' Barney, David J., and Brian F. Schaffner. 
#'  "Reexamining the Effect of Mass Shootings on Public Support for Gun Control." 
#'   British Journal of Political Science 49.4 (2019): 1555-1565. 
#'   \url{https://doi.org/10.1017/S0007123418000352}.
#'
#' Newman, B. J. and Hartman, T. K.
#'  ``Mass shootings and public support for gun control.''
#'  British Journal of Political Science, 49.4 (2019): 1527–1553.
#'  \url{https://doi.org/10.1017/S0007123417000333}.
#' @keywords dataset 
#' @examples
#'  require(dplyr)
#'  require(haven)
#'  data("gun_threewave")
#'  gun_threewave
"gun_threewave"


#' Two-wave panel with covariates
#'
#' A subset of approximately 2000 individuals from the CCES 2010-12 panel,
#' augmented with IRT ideology scores and demographic variables.
#' The ideology scores are estimated via a two-parameter logistic (2PL) IRT model
#' fitted to CCES policy items.  This dataset is designed for demonstrating the
#' semiparametric estimator, which requires at least one continuous covariate.
#'
#' @name gun_twowave_cov
#' @docType data
#' @format A data frame with 4000 rows and 10 columns.
#'   \describe{
#'      \item{caseid}{Unique id of a respondent.}
#'      \item{year}{Year of the response (2010 or 2012).}
#'      \item{guns}{Attitudes towards gun control: (1) "Less Strict", (2) "Kept As They Are", (3) "More Strict".}
#'      \item{treat_25mi}{Treatment indicator: 1 if a mass shooting happened within 25 miles.}
#'      \item{treat_100mi}{Treatment indicator: 1 if a mass shooting happened within 100 miles.}
#'      \item{reszip}{Zipcode of a respondent.}
#'      \item{ideology}{First-dimension IRT ideology score (continuous).}
#'      \item{ideology2}{Second-dimension IRT ideology score (continuous).}
#'      \item{female}{Binary indicator for female respondent.}
#'      \item{educ}{Education level (integer scale).}
#'   }
#' @source
#' Barney, David; Schaffner, Brian.
#'  ``Replication Data for: Reexamining the Effect of Mass Shootings on Public Support for Gun Control.''
#'  \url{https://doi.org/10.7910/DVN/YJQIXP}, Harvard Dataverse, V1, 2018.
#'
#'  Ansolabehere, Stephen; Schaffner, Brian.
#'   ``2010 - 2012 CCES Panel Study.''
#'   \url{https://doi.org/10.7910/DVN/24416}, Harvard Dataverse, V4, 2014.
#'
#' @references
#' Barney, David J., and Brian F. Schaffner.
#'  "Reexamining the Effect of Mass Shootings on Public Support for Gun Control."
#'   British Journal of Political Science 49.4 (2019): 1555-1565.
#'   \url{https://doi.org/10.1017/S0007123418000352}.
#' @keywords dataset
#' @examples
#'  data("gun_twowave_cov")
#'  head(gun_twowave_cov)
"gun_twowave_cov"
