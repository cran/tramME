##' Predict method for tramME objects
##'
##' Evaluates the _conditional_ distribution implied by a tramME model, given by a
##' set of covariates and random effects on a desired scale.
##
##' When \code{newdata} contains values of the response variable, prediction is only
##' done for those values. When no response values are supplied, prediction is done on
##' a grid of values.
##
##' Unfitted tramME models can also be used for prediction as long as the coefficent
##' parameter are set manually (with \code{coef<-}).
##'
##' When \code{ranef} is equal to "zero", a vector of zeros with the right size is
##' substituted.
##' @param object A tramME object
##' @param ranef Vector of random effects or the word "zero". See details.
##' @param ... Additional arguments, passed to \code{\link[mlt]{predict.mlt}}.
##' @inheritParams mlt::predict.ctm
##' @return A numeric matrix of the predicted values invisibly
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' predict(fit, type = "trafo") ## evaluate on the transformation function scale
##' nd <- sleepstudy
##' nd$Reaction <- NULL
##' pr <- predict(fit, newdata = nd, ranef = ranef(fit, raw = TRUE), type = "distribution",
##'               K = 100)
##' @importFrom stats predict
##' @export
predict.tramME <- function(object, newdata = NULL, ranef = NULL, ...) {
  if (!object$fitted && any(is.na(object$pars$coef)))
    stop("Coefficient values must be set before calling predict")
  if (is.null(ranef) && is.null(newdata)) {
    if (!object$fitted) ## NOTE: when unfitted, newdata must be supplied
      stop("newdata must be supplied for unfitted models")
    ranef <- ranef(object, raw = TRUE)
    newdata <- object$data$mf
  }
  stopifnot(!is.null(newdata) && !is.null(ranef))
  rsiz <- .re_size(object$model, newdata)
  if (identical(ranef, "zero"))
    ranef <- rep(0, sum(rsiz$bsize * rsiz$nlev))
  stopifnot(sum(rsiz$bsize * rsiz$nlev) == length(ranef)) ## check if the vector of random effects has enough distinct values
  re_struct <- .re_data(object$call$formula[-2L], data = newdata,
                        negative = FALSE) ## NOTE: set in .dummy_ctm
  ## TODO: .re_data drops unused levels of the grouping factor, which is what we want(?)
  ## Since it is an unexpected behavior, check if it works with more complex random effects
  ## structures (e.g. crossed)
  newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, ranef))
  ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
  fmlt <- .dummy_ctm(object$model, coef(object, with_baseline = TRUE))
  out <- predict(fmlt, newdata = newdata, ...)
  return(out)
}


##' Plotting method for tramME objects
##'
##' Plot the conditional distribution evaluated at a grid of possible response
##' values and a set of covariate and random effects values on a specified scale.
##'
##' When \code{ranef} is equal to "zero", a vector of zeros with the right size is
##' substituted.
##' @param x A tramME object
##' @param ranef Vector of random effects or the word "zero". See details.
##' @param ... Additional arguments, passed to \code{\link[mlt]{plot.mlt}}.
##' @inheritParams mlt::predict.ctm
##' @return A numeric matrix of the predicted values invisibly
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' plot(fit, K = 100, type = "density")
##' @importFrom graphics plot
##' @export
plot.tramME <- function(x, newdata = NULL, ranef = NULL, ...) {
  if (!x$fitted && any(is.na(x$pars$coef)))
    stop("Coefficient values must be set before calling predict")
  if (is.null(ranef) && is.null(newdata)) {
    if (!x$fitted) ## NOTE: when unfitted, newdata must be supplied
      stop("newdata must be supplied for unfitted models")
    ranef <- ranef(x, raw = TRUE)
    newdata <- x$data$mf
  }
  stopifnot(!is.null(newdata) && !is.null(ranef))
  rsiz <- .re_size(x$model, newdata)
  if (identical(ranef, "zero"))
    ranef <- rep(0, sum(rsiz$bsize * rsiz$nlev))
  stopifnot(sum(rsiz$bsize * rsiz$nlev) == length(ranef)) ## check if the vector of random effects has enough distinct values
  re_struct <- .re_data(x$call$formula[-2L], data = newdata,
                        negative = FALSE) ## NOTE: set in .dummy_ctm
  newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, ranef))
  ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
  fmlt <- .dummy_ctm(x$model, coef(x, with_baseline = TRUE))
  invisible(plot(fmlt, newdata = newdata, ...))
}
