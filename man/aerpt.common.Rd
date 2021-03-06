% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aerpt.common.R
\name{aerpt.common}
\alias{aerpt.common}
\title{Create listing of commomly occurring AEs for Alliance DSMB reports}
\usage{
aerpt.common(
  cytox,
  id,
  threshold = 10,
  grade.cutoff = 3,
  rel.cutoff = 0,
  rel.include = 1,
  byvar = "arm",
  arm.labels = NULL,
  toxvar = "toxicity",
  relvar = "rel_smed",
  gradevar = "grade"
)
}
\arguments{
\item{cytox}{The toxicity dataset containing multiple observations per patient with one observation per AE > grade 0}

\item{id}{The patient ID variable used in the dataset}

\item{threshold}{Numeric threshold level for commonly occurring Grade 3+ AEs. Default to 10}

\item{grade.cutoff}{Grade cutoff for minimum AE grade. Inclusive}

\item{rel.cutoff}{Attribution cutoff for minimum AE attribution. Inclusive}

\item{rel.include}{Value of 0 or 1. If 0, "Regardless of Attribution" printed in header. If 1, "At least possibly related" printed in header}

\item{byvar}{Variable used to separate summary by arm. If not present, "All Patients" printed for arm. Set to NULL if calculations are desired for all patients together.}

\item{arm.labels}{Option to include additional Arm labels in the table header. Add as vector of characters.}

\item{toxvar}{Variable containing numeric toxicity codes.}

\item{relvar}{Variable containing numeric relation values 1-5.}

\item{gradevar}{Grade variable, numeric.}
}
\value{
This function creates a listing of the AE events for the Alliance DSMB report
}
\description{
This function is used to create commonly occurring AE listing output for Alliance DSMB reports.
}
\examples{

# With data loaded for A041501

crtlibn(d = 'A041501')

all.common <- aerpt.common(cytox, id = "dcntr_id", byvar = "arm", 
grade.cutoff = 3, rel.cutoff = 0, rel.include = 0)
all.common$grade3
all.common$grade45

}
\author{
Sawyer Jacobson
}
