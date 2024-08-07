% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confband.R
\name{confband.tramME}
\alias{confband.tramME}
\title{Confidence intervals and bands from a \code{tramME} model}
\usage{
\method{confband}{tramME}(
  object,
  newdata,
  level = 0.95,
  type = c("trafo", "distribution", "survivor", "cumhazard"),
  adjust = FALSE,
  K = 40,
  cheat = K,
  q = NULL,
  baseline_only = FALSE,
  ...
)
}
\arguments{
\item{object}{The \code{tramME} object.}

\item{newdata}{A data frame of covariate values.}

\item{level}{Confidence level.}

\item{type}{The scale on which the condfidence bands are calculated.}

\item{adjust}{If \code{TRUE}, multiplicity-adjusted confidence bands are
calculated. (see Details)}

\item{K}{The number of grid points at which the outcome distribution is
evaluated.}

\item{cheat}{In the case of multiplicity adjustment (\code{adjust = TRUE}),
an option to decrease the number of grid points (\code{cheat < K}), for
faster calculations and increased numerical stability. (see Details)}

\item{q}{The quantiles at which the model is evaluated.}

\item{baseline_only}{If \code{TRUE}, only evaluate the baseline
transformation function and ignore the shift terms.}

\item{...}{Optional arguments passed to \code{\link[multcomp]{confint.glht}}.}
}
\value{
A matrix (in the case when \code{newdata} has a single row) or a
  list of matrices for each row of \code{newdata}.
}
\description{
Pointwise confidence intervals or multiplicity-adjusted confidence bands for
    transformation, distribution, survivor or cumulative hazard functions.
}
\details{
Similarly to \code{\link[mlt]{confband}}, this method evaluates the
  conditional distribution of the outcome on a selected scale given a number
  of grid-points and claculates the corresponding confidence intervals or
  bands (adjusting for multiplicity).

  The point estimates retured by this function could also be calculated with
  \code{\link{predict.tramME}} (when \code{newdata} does not contain
  response values and \code{K} is set to the number of grid points).
  While \code{\link{predict.tramME}} is designed to calculate a
  potentially large number of point estimates on a wider range of available
  scales, \code{confband} calculates the asymptotic intervals from the joint
  covariance matrix of the fixed and random effects. For technical reasons,
  a smaller set of \code{type} options are available, and the calculations
  are slower than with \code{\link{predict.tramME}}. The handling of random
  effects is currently stricter than in \code{predict.tramME}: No
  \code{ranef} option is available, and grouping factors for random effects
  supplied in \code{newdata} must have the same levels as the dataset used
  to fit the model.

  The multiplicity adjustment is done using
  \code{\link[multcomp]{confint.glht}}. The \code{cheat} argument reduces
  the dimensionality of the multivariate root-finding problem (see
  \code{\link[mvtnorm]{qmvt}}) for speed and (occasionally) numerical
  stability. The critical value for the confidence bands are obtained for
  \code{cheat < K} grid points, but the confidence bands are calculated for
  \code{K} grid points. As a result, the nominal level of the returned
  confidence band is not maintained, but the deviation is expected to be
  small if \code{cheat} is reasonably large. It is the user's responsibility
  to set this value, and by default \code{cheat = K}.
}
\section{Warning}{


  This method implements new functionality. Its user interface may
  be subject to change.
}

