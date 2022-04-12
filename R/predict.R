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
##' \code{ranef} can be different objects based on what we want to calculate and
##' what the other inputs are. If \code{ranef} is a \code{ranef.tramME}, we assume
##' that it contains the full set of random effects, but not the penalized coefficients
##' of the smooth terms. In this case \code{fix_smooth} must be \code{TRUE}. If
##' \code{ranef} is a named vector, we are fixing the supplied random effects (and
##' penalized coefficients) and predict the rest from \code{newdata} (\code{fix_smooth}
##' may also be used in this case). In this case, the random effects are identified
##' with the same naming convention as in \code{object$param$gamma}.
##'
##' If \code{ranef} is an unnamed vector, the function expects the
##' full set of necessary random effects (with or without penalized coefficients, depending
##' on \code{fix_smooth}). If \code{ranef = NULL} (the default), all random effects and
##' optionally penalized parameters (although this is not recommended) are predicted from
##' \code{newdata}. Finally, if \code{ranef} is equal to "zero", a vector of zeros with the
##' right size is used.
##' @param object A \code{tramME} object.
##' @param ranef Random effects it can be a \code{ranef.tramME} object, a named list,
##'   an unnamed list, \code{NULL} or the word "zero". See Details.
##' @param fix_smooth If \code{FALSE}, the random effects coefficients of the smooth
##'   terms are refitted to \code{newdata}. It's probably not what you want to do.
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
## TODO: later at least some of the predictions could run through predict.tramTMB
## to obtain the standard errors
predict.tramME <- function(object, newdata = model.frame(object),
  ranef = NULL, fix_smooth = TRUE,
  type = c("lp", "trafo", "distribution", "survivor", "density",
           "logdensity", "hazard", "loghazard", "cumhazard",
           "odds", "logodds", "quantile"), ...) {
  type <- match.arg(type)
  rfs <- .ranef_predict_setup(object, newdata, ranef, fix_smooth)
  Zt  <- rfs$Zt
  ranef <- rfs$ranef

  X <- model.matrix(object, data = newdata, type = "X",
                    keep_sign = FALSE, ignore_response = TRUE)$X
  if (length(ranef)) Zg <- as.numeric(Matrix::crossprod(Zt, ranef))
  else Zg <- 0
  cf <- object$param$beta
  bty <- attr(cf, "type")
  cfs <- cf[bty != "bl"]

  ## -- linear predictor
  if (identical(type, "lp")) {
    out <- as.numeric(X %*% cfs) + Zg
    names(out) <- rownames(newdata)
    return(out)
  }
  ## -- other prediction types: wrapper around predict.mlt
  if (length(idx <- which(bty[bty != "bl"] == "sm"))) {
    sXb <- as.numeric(X[, idx, drop = FALSE] %*% cfs[idx])
  } else sXb <- 0

  if (any((re_ <- sXb + Zg) != 0)) {
    newdata$re_ <- re_
    ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
    fmlt <- .cctm(object$model$ctm, cf[bty != "sm"],
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
##' substituted. For more details, see \code{\link[tramME]{predict.tramME}}.
##'
##' For more information on how to control the grid on which the functions are evaluated,
##' see the documentation of \code{\link[mlt]{predict.mlt}}.
##' @param x A \code{tramME} object.
##' @param ranef Random effects (either in named list format or a numeric vector)
##'   or the word "zero". See Details.
##' @param fix_smooth If \code{FALSE}, the random effects coefficients of the smooth
##'   terms are refitted to \code{newdata}. It's probably not what you want to do.
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
plot.tramME <- function(x, newdata = model.frame(x),
  ranef = NULL, fix_smooth = TRUE,
  type = c("trafo", "distribution", "survivor", "density",
           "logdensity", "hazard", "loghazard", "cumhazard",
           "odds", "logodds", "quantile"), ...) {
  type <- match.arg(type)

  rfs <- .ranef_predict_setup(x, newdata, ranef, fix_smooth)
  Zt  <- rfs$Zt
  ranef <- rfs$ranef

  X <- model.matrix(x, data = newdata, type = "X",
                    keep_sign = FALSE, ignore_response = TRUE)$X
  if (length(ranef)) Zg <- as.numeric(Matrix::crossprod(Zt, ranef))
  else Zg <- 0
  cf <- x$param$beta
  bty <- attr(cf, "type")
  cfs <- cf[bty != "bl"]

  if (length(idx <- which(bty[bty != "bl"] == "sm"))) {
    sXb <- as.numeric(X[, idx, drop = FALSE] %*% cfs[idx])
  } else sXb <- 0

  if (any((re_ <- sXb + Zg) != 0)) {
    newdata$re_ <- re_
    ## NOTE: create a dummy mlt model where REs enter as fixed parameter FEs
    fmlt <- .cctm(x$model$ctm, cf[bty != "sm"],
                  negative = x$model$negative)
  } else {
    fmlt <- x$model$ctm
    coef(fmlt) <- coef(x, with_baseline = TRUE, fixed = TRUE)
  }
  invisible(plot(fmlt, newdata = newdata, type = type, ...))
}


## Helper function to handle the different ways ranef can be supplied to
## \code{predict.tramME} and \code{plot.tramME}.
## It sets up the random effects vector, either by reading out saved results,
## or by predicting elements, using \code{ranef.tramME}.
## @param object A \code{tramME} object.
## @param newdata A \code{data.frame} with the new data points. May or may not contain
##   values for the outcome.
## @param ranef Random effects it can be a \code{ranef.tramME} object, a named list,
##   an unnamed list, \code{NULL} or the word "zero". See Details.
## @param fix_smooth If \code{FALSE}, the random effects coefficients of the smooth
##    terms are refitted to \code{newdata}. It's probably not what you want to do.
.ranef_predict_setup <- function(object, newdata, ranef, fix_smooth) {
  Zt <- model.matrix(object, data = newdata, type = "Zt", keep_sign = FALSE,
                     drop_unused_groups = TRUE, ignore_response = TRUE)$Zt
  ty <- attr(Zt, "type")
  pn <- attr(Zt, "parnames")
  reo <- rep(NA, length(pn))
  names(reo) <- pn
  ## 1) As a ranef.tramME object: it contains _all_ random effects, but not the
  ## penalized parameters of the smooth terms (fix_smooth should be set)
  ## Convert it to a named vector.
  if (inherits(ranef, "ranef.tramME")) {
    ## do name checks
    ren <- pn[ty == "re"]
    gn <- unique(sapply(strsplit(ren, "|", fixed = TRUE), `[`, 1))
    stopifnot(identical(gn, names(ranef)))
    sn <- unique(sapply(strsplit(ren, ":", fixed = TRUE), `[`, 2))
    stopifnot(identical(sn, as.vector(sapply(ranef, rownames))))
    vn <- unique(sapply(strsplit(ren, "|", fixed = TRUE), function(n) {
      strsplit(n, ":", fixed = TRUE)[[2]][1]
    }))
    stopifnot(identical(vn, as.vector(sapply(ranef, colnames))))
    reo[ren] <- unname(unlist(lapply(ranef, function(x) c(t(x)))))
  }
  ## 2) "zero"
  if (identical(ranef, "zero")) {
    reo[] <- 0
  }
  gty <- attr(object$param$gamma, "type")
  if (fix_smooth) {
    reo[pn[ty == "sm"]] <- unname(object$param$gamma[gty == "sm"])
  }
  if (is.numeric(ranef)) {
    if (!length(ren <- names(ranef))) {
      ## 3) Unnamed vector: Should contain all necessary values (and nothing more)
      ## There is nothing to predict in this case.
      stopifnot(length(ranef) == sum(idx <- is.na(reo)))
      reo[idx] <- ranef
    } else {
      ## 4) Named vector: fix the random effects & penalized coefs _by name_.
      ## Predict the remainder. Ambiguities are not tolerated.
      stopifnot(all(ren %in% pn))
      reo[ren] <- ranef
    }
  }
  if (any(is.na(reo))) {
    gg <- reo[!is.na(reo)]
    if (!length(gg)) gg <- NULL
    reo[] <- ranef(object, newdata = newdata, param = list(gamma = reo[!is.na(reo)]),
                       raw = TRUE, fix_smooth = fix_smooth)
  }
  if (any(is.na(reo))) {
    stop(paste("Could not calculate random effects vector.",
               "Please set up the model properly or supply",
               "random effects manually."))
  }
  return(list(Zt = Zt, ranef = reo))
}
