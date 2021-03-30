##' Predict method for tramME objects
##'
##' Evaluates the _conditional_ distribution implied by a tramME model, given by a
##' set of covariates and random effects on a selected scale.
##'
##' When \code{newdata} contains values of the response variable, prediction is only
##' done for those values. In this case, if random effects vector (\code{ranef}) is not
##' supplied by the user, the function predicts the random effects from the model
##' using \code{newdata}.
##'
##' When no response values are supplied in \code{newdata}, the prediction is done
##' on a grid of values for each line of the dataset (see \code{\link[mlt]{predict.mlt}}
##' for information on how to control the setup of this grid).
##' In this case, the user has to specify the vector of random effects to avoid ambiguities.
##'
##' The linear predictor (\code{type = "lp"}) equals to the shift terms plus the random
##' effects terms _without the baseline transfromation function_.
##'
##' The linear predictor (\code{type = "lp"}) and the conditional quantile function
##' (\code{type = "quantile"}) are special in that they do not return results evaluated
##' on a grid, even when the response variable in \code{newdata} is missing. The probabilities
##' for the evaluation of the quantile function can be supplied with the \code{prob} argument
##' of \code{\link[mlt]{predict.mlt}}.
##'
##' In the case of \code{type = "quantile"}, when the some of the requested conditonal
##' quantiles fall outside of the support of the response distribution
##' (specified when the model was set up), the inversion of the CDF cannot be done exactly
##' and \code{tramME} returns censored values.
##'
##' When \code{ranef} is equal to "zero", a vector of zeros with the right size is
##' used.
##' @param object A \code{tramME} object.
##' @param ranef Random effects (either in named list format or a numeric vector)
##'   or the word "zero". See Details.
##' @param type The scale on which the predictions are evaluated:
##'   \itemize{
##'     \item lp: Linear predictor (Xb + Zg). For more information, see Details.
##'     \item trafo: The prediction evaluated on the scale of the
##'       transformation function.
##'     \item distribution: The prediction evaluated on the scale of the
##'       conditional CDF.
##'     \item survivor: The prediction evaluated on the scale of the
##'       (conditional) survivor function.
##'     \item density, logdensity: The prediction evaluated on the scale of
##'       the conditional (log-)PDF.
##'     \item hazard, loghazard, cumhazard: The prediction evaluated on the
##'       hazard/log-hazard/cumulative hazard scale.
##'     \item odds, logodds: The prediction evaluated on the (log-)odds scale.
##'     \item quantile: Return the quantiles of the conditional outcome distribution
##'       corresponding to \code{newdata}. For more information, see Details.
##'   }
##' @param ... Additional arguments, passed to \code{\link[mlt]{predict.mlt}}.
##' @inheritParams mlt::predict.ctm
##' @return A numeric vector/matrix of the predicted values (depending on the inputs)
##'   or a \code{response} object, when the some of the requested conditonal quantiles
##'   fall outside of the support of the response distribution specified when the model
##'   was set up (only can occur with \code{type = "quantile"}).
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' predict(fit, type = "trafo") ## evaluate on the transformation function scale
##' nd <- sleepstudy
##' nd$Reaction <- NULL
##' pr <- predict(fit, newdata = nd, ranef = ranef(fit), type = "distribution",
##'               K = 100)
##' @importFrom stats predict
##' @export
predict.tramME <- function(object, newdata = model.frame(object), ranef = NULL,
  type = c("lp", "trafo", "distribution", "survivor", "density",
           "logdensity", "hazard", "loghazard", "cumhazard",
           "odds", "logodds", "quantile"), ...) {
  type <- match.arg(type)

  ## NOTE: in case ranef is formatted, transfrom it to a vector
  if (is.list(ranef)) {
    ranef <- unname(unlist(lapply(ranef, function(x) c(t(x)))))
  }

  if (is.null(ranef) && !is.null(object$model$ranef)) {
    rn <- variable.names(object, "response")
    if (rn %in% names(newdata)) {
      ranef <- ranef(object, newdata = newdata, raw = TRUE)
    } else {
      ## FIXME: could try harder and assume that the ranef vector is the one
      ## fitted in the model if the levels of the grouping variable are the same
      stop("The random effects vector must be specified.")
    }
  }
  if (any(is.na(ranef))) {
    stop(paste("Could not calculate random effects vector.",
               "Please set up the model properly or supply",
               "random effects manually."))
  }

  ## -- set up and check random effects vector
  if (!is.null(object$model$ranef)) {
    rsiz <- .re_size(attr(object$param, "re")$blocksize, newdata)
    if (identical(ranef, "zero"))
      ranef <- rep(0, sum(rsiz$bsize * rsiz$nlev))
    ## NOTE: check if the vector of random effects has enough distinct values
    stopifnot(sum(rsiz$bsize * rsiz$nlev) == length(ranef))
  }

  ## -- linear predictor
  if (type == "lp") {
    X <- model.matrix(object, data = newdata, type = "fixef", with_baseline = FALSE)
    Zt <- model.matrix(object, data = newdata, type = "ranef")
    if (is.null(Zt)) {
      Zg <- 0
    } else {
      Zg <- as.numeric(Matrix::crossprod(Zt, ranef))
    }
    out <- as.numeric(X %*% coef(object, with_baseline = FALSE)) + Zg
    if (object$model$negative) out <- -out
    names(out) <- rownames(newdata)
    return(out)
  }

  ## -- other prediction types: wrapper around predict.mlt
  if (!is.null(object$model$ranef)) {
    ## FIXME: change this to model.matrix( , type = "ranef")
    re_struct <- re_terms(object$model$ranef, data = newdata,
                          negative = FALSE) ## NOTE: set in .cctm
    newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, ranef))
    ## ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
    fmlt <- .cctm(object$model$ctm, coef(object, with_baseline = TRUE, fixed = TRUE),
                  negative = object$model$negative)
  } else {
    fmlt <- object$model$ctm
    coef(fmlt) <- coef(object, with_baseline = TRUE, fixed = TRUE)
  }
  out <- predict(fmlt, newdata = newdata, type = type, ...)
  return(out)
}


##' Plotting method for tramME objects
##'
##' Plot the conditional distribution evaluated at a grid of possible response
##' values and a set of covariate and random effects values on a specified scale.
##'
##' When \code{ranef} is equal to "zero", a vector of zeros with the right size is
##' substituted.
##'
##' For more information on how to control the grid on which the functions are evaluated,
##' see the documentation of \code{\link[mlt]{predict.mlt}}.
##' @param x A \code{tramME} object.
##' @param ranef Random effects (either in named list format or a numeric vector)
##'   or the word "zero". See Details.
##' @param type The scale on which the predictions are evaluated:
##'   \itemize{
##'     \item trafo: The prediction evaluated on the scale of the
##'       transformation function.
##'     \item distribution: The prediction evaluated on the scale of the
##'       conditional CDF.
##'     \item survivor: The prediction evaluated on the scale of the
##'       (conditional) survivor function.
##'     \item density, logdensity: The prediction evaluated on the scale of
##'       the conditional (log-)PDF.
##'     \item hazard, loghazard, cumhazard: The prediction evaluated on the
##'       hazard/log-hazard/cumulative hazard scale.
##'     \item odds, logodds: The prediction evaluated on the (log-)odds scale.
##'     \item quantile: Return the quantiles of the conditional outcome distribution
##'       corresponding to \code{newdata}. For more information, see Details.
##'   }
##' @param ... Additional arguments, passed to \code{\link[mlt]{plot.mlt}}.
##' @inheritParams mlt::predict.ctm
##' @return A numeric matrix of the predicted values invisibly.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' plot(fit, K = 100, type = "density")
##' @importFrom graphics plot
##' @export
plot.tramME <- function(x, newdata = model.frame(x), ranef = NULL,
  type = c("trafo", "distribution", "survivor", "density",
           "logdensity", "hazard", "loghazard", "cumhazard",
           "odds", "logodds", "quantile"), ...) {
  type <- match.arg(type)

  ## NOTE: in case ranef is formatted, transfrom it to a vector
  if (is.list(ranef)) {
    ranef <- unname(unlist(lapply(ranef, function(x) c(t(x)))))
  }

  if (is.null(ranef) && !is.null(x$model$ranef)) {
    rn <- variable.names(x, "response")
    if (rn %in% names(newdata)) {
      ranef <- ranef(x, newdata = newdata, raw = TRUE)
    } else {
      ## FIXME: could try harder and assume that the ranef vector is the one
      ## fitted in the model if the levels of the grouping variable are the same
      stop("The random effects vector must be specified.")
    }
  }
  if (any(is.na(ranef))) {
    stop(paste("Could not calculate random effects vector.",
               "Please set up the model properly or supply",
               "random effects manually."))
  }

  if (!is.null(x$model$ranef)) {
    rsiz <- .re_size(attr(x$param, "re")$blocksize, newdata)
    if (identical(ranef, "zero"))
      ranef <- rep(0, sum(rsiz$bsize * rsiz$nlev))
    stopifnot(sum(rsiz$bsize * rsiz$nlev) == length(ranef)) ## check if the vector of random effects has enough distinct values

    ## FIXME: change this to model.matrix( , type = "ranef")
    re_struct <- re_terms(x$model$ranef, data = newdata,
                          negative = FALSE) ## NOTE: set in .dummy_ctm
    newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, ranef))
    ## ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
    fmlt <- .cctm(x$model$ctm, coef(x, with_baseline = TRUE, fixed = TRUE),
                  negative = x$model$negative)
  } else {
    fmlt <- x$model$ctm
    coef(fmlt) <- coef(x, with_baseline = TRUE, fixed = TRUE)
  }
  invisible(plot(fmlt, newdata = newdata, type = type, ...))
}
