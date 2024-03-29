% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth.R
\name{plot.smooth.tramME}
\alias{plot.smooth.tramME}
\title{Plot smooth terms of a tramME model.}
\usage{
\method{plot}{smooth.tramME}(
  x,
  which = seq_along(x),
  col = 1,
  fill = grey(0.5, 0.25),
  trafo = I,
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{A \code{smooth.tramME} object.}

\item{which}{Select terms to be printed by their indices}

\item{col}{Line color for the point estimates.}

\item{fill}{Fill color for the confidence intervals.}

\item{trafo}{Monotonic transformation to be applied on the smooth terms}

\item{add}{Add the plot to an existing figure.}

\item{...}{Optional parameters passed to the plotting functions.}
}
\description{
Plot smooth terms of a tramME model.
}
\examples{
data("mcycle", package = "MASS")
fit <- LmME(accel ~ s(times), data = mcycle)
plot(smooth_terms(fit, as.lm = TRUE))
}
