##' Mixed-effects version of \code{\link[tram]{Survreg}}
##' @inheritParams tram::Survreg
##' @inheritParams LmME
##' @param silent logical, make TMB functionality silent
##' @param resid logical, Should the score residuals also be calculated?
##' @param estinit logical, estimate a vector of initial values for the fixed effects parameters
##'   from a (fixed effects only) mlt model
##' @param initpar named list of initial parameter values, if \code{NULL}, it is ignored
##' @inheritParams mlt::mlt
##' @param nofit logical, if TRUE, creates the model object, but does not run the optimization
##' @param control list with controls for optimization
##' @return A \code{SurvregME} object.
##' @importFrom tram Survreg
##' @export
SurvregME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      dist = c("weibull", "logistic", "gaussian", "exponential",
                               "rayleigh", "loggaussian", "lognormal", "loglogistic"),
                      scale = 0,
                      silent = TRUE, resid = FALSE, do_update = FALSE,
                      estinit = TRUE, initpar = NULL,
                      fixed = NULL, nofit = FALSE,
                      control = optim_control(),
                      ...) {
  cl <- match.call()
  cl$call <- cl
  cl$dist <- match.arg(dist)
  cl$scale <- scale
  cl[[1L]] <- quote(tramME)
  cl$tram <- "Survreg"
  eval(cl, parent.frame())
}


##' Extract the coefficients of the fixed effects terms of an SurvregME model.
##' @param object An \code{SurvregME} object.
##' @param as.survreg If \code{TRUE}, return the transformed coefficients as in a
##'   \code{survival::survreg} object.
##' @inheritParams coef.LmME
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' library("survival")
##' fit <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats)
##' coef(fit, as.survreg = TRUE)
##' @importFrom stats coef
##' @export
coef.SurvregME <- function(object, as.survreg = FALSE, ...) {
  coef.LmME(object, as.lm = as.survreg, ...)
}
