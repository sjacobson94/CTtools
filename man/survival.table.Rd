% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival.table.R
\name{survival.table}
\alias{survival.table}
\title{Survival table}
\usage{
survival.table(
  data,
  stat,
  time,
  label = "Rate",
  arm = "arm",
  ref = as.character(levels(data2$arm)[1]),
  rate_span = c(12),
  time_level = "month",
  rate_label = "PFS",
  one.sided = FALSE
)
}
\arguments{
\item{data}{A dataset containing columns necessary for survival analysis.}

\item{stat}{Survival event variable}

\item{time}{Survival time variable, currently only supports months}

\item{label}{Label for survival variable}

\item{arm}{Arm variable}

\item{ref}{Option to specify the reference arm. Must be entered how it is in the data.}

\item{rate_span}{Vector of time levels desired for survival rate. Used in conjunction with time_level.}

\item{time_level}{Time level desired. Currently supports "months" and "years".}

\item{rate_label}{Label for survival rate in the table (eg. rate_label = "PFS" gives "12-month PFS rate).}

\item{one.sided}{Logical. If TRUE returns one-sided p-values. Default is FALSE}
}
\value{
This function returns a nicely formatted table of survival analysis statistics.
}
\description{
This function creates a table/dataset with number of events/N, 95\% median (months) survival CI, 95\% survival rate CI (eg. 12-month PFS rate), and 95\% HR CI.
}
