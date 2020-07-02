#' The CTtools package
#'
#' This package is meant to replicate some of the useful SAS macros used in the Cancer Center to make R more friendly and easier to use
#'
#' @section Functions:
#'
#' Below are a list of some of the most useful functions in \code{CTtools}.
#'
#' \code{\link{accgph}}: Accrual graph
#'
#' \code{\link{acc.quarter}}: Quarterly accrual table
#'
#' \code{\link{eot.table}}: End of treatment table for Alliance DSMB report
#'
#' \code{\link{aerpt.summary}}: AE summary table for Alliance DSMB report
#'
#' \code{\link{aerpt.listing}}: AE listing table for Alliance DSMB report
#'
#' \code{\link{irbtox_summary}}: AE summary table for Mayo DSMB report
#'
#' \code{\link{irbtox_listing}}: AE listing table for Mayo DSMB report
#'
#' \code{\link{survival.table}}: Table for survival analysis output
#'
#' \code{\link{forest_plot}}: Creates a forestplot for survival meta-analysis
#'
#' \code{\link{km.plot}}: Creates a publication quality KM plot
#'
#' \code{\link{tidy_tableby}}: Creates a nicely formatted tableby() summary when knitting to .docx using flextable
#'
#' @docType package
#' @name CTtools
NULL

#### commands to build the package using devtools
# devtools::check_man()
# devtools::test()
# devtools::check()
# withr::with_libpaths("/projects/ccsstat/across/R/R_packages/", devtools::install(), action = "prefix")
# devtools::build() # builds in home directory as well so the package can be used in shiny apps
# devtools::install() # installs in local directory
# library(CTtools, lib.loc = "/projects/ccsstat/across/R/R_packages/")



#' @importFrom magrittr `%>%`
#'
#' @export
magrittr::`%>%`

#' @importFrom lubridate `%m+%`
#'
#' @export
lubridate::`%m+%`

#' @importFrom lubridate `%m-%`
#'
#' @export
lubridate::`%m-%`

#' @import ggplot2
#' @import dplyr
#' @import flextable
NULL
