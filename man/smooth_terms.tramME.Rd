% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth.R
\name{smooth_terms.tramME}
\alias{smooth_terms.tramME}
\alias{smooth_terms}
\title{Extract and evaluate the smooth terms of a tramME model}
\usage{
\method{smooth_terms}{tramME}(object, k = 100, newdata = NULL, ...)
}
\arguments{
\item{object}{A \code{tramME} object.}

\item{k}{Integer, the number of points to be used to evaluate the smooth terms.
Ignored when \code{newdata} is supplied.}

\item{newdata}{A \code{data.frame} with new values for the smooth terms.
If \code{NULL}, the new data is set up based on the \code{model.frame} and
\code{k}. Smooths for which the supplied information in this input is incomplete
will be ignored.}

\item{...}{Optional arguments. \code{as.lm} is passed through this when it is necessary.}
}
\value{
A list of results from evaluating the smooth terms of the model.
}
\description{
Extract and evaluate the smooth terms of a tramME model
}
\examples{
data("mcycle", package = "MASS")
fit <- LmME(accel ~ s(times), data = mcycle)
plot(smooth_terms(fit))
}
