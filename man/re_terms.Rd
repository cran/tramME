% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tramTMB.R
\name{re_terms}
\alias{re_terms}
\title{Create random effects data and initial paramaters}
\usage{
re_terms(ranef, data, negative)
}
\arguments{
\item{ranef}{a list of random effects formulas from \code{\link[lme4]{findbars}}}

\item{data}{data.frame containing the variables of the model}

\item{negative}{logical value that indicates whether the random effects have
a negative sign}
}
\value{
A list containing data and parameter values to be used in the TMB model.
}
\description{
Create random effects data and initial paramaters
}
