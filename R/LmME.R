##' Mixed-effects Additive Normal Linear Regression Model
##'
##' Estimates the normal linear model parameterized as a linear transformation
##' model.
##'
##' @inheritParams tramME
##' @details
##'
##' The additive mixed-effects normal linear model is a special case of the
##'   mixed-effects additive transformation model family, with the
##'   transformation function restricted to be linear and the inverse link set
##'   to the standard Gaussian CDF (see Hothorn et al., 2018). This function
##'   estimates this model with the transformation model parameterization, and
##'   offers features that are typically not available in other mixed-effects
##'   additive implementations, such as stratum-specific variances, and censored
##'   and/or truncated observations.
##'
##' The model extends \code{\link[tram:Lm]{tram::Lm}} with random effects and
##'   (optionally penalized) additive terms. For details on mixed-effect
##'   transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with negative
##'   parameters (i.e. \code{negative = TRUE} in \code{\link[tram]{tram}}).
##'
##' The results can be transformed back to the usual linear mixed/additive model
##'   parametrization with specific methods provided by \code{tramME}. The
##'   differences between the two parametrizations are discussed in Tamasi and
##'   Hothorn (2021).
##'
##' @inherit tramME references
##' @return A \code{LmME} model object.
##' @examples
##' library("survival")
##' data("sleepstudy", package = "lme4")
##' ## Create a version of the response with 200 ms detection limit and 50 ms
##' ## step sizes
##' ub <- ceiling(sleepstudy$Reaction / 50) * 50
##' lb <- floor(sleepstudy$Reaction / 50) * 50
##' lb[ub == 200] <- 0
##' sleepstudy$Reaction_ic <- Surv(lb, ub, type = "interval2")
##' m <- LmME(Reaction_ic ~ Days + (Days | Subject), data = sleepstudy)
##' summary(m)
##' coef(m, as.lm = TRUE)
##' @importFrom tram Lm
##' @export
LmME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                 silent = TRUE, resid = FALSE, do_update = FALSE,
                 estinit = TRUE, initpar = NULL,
                 fixed = NULL, nofit = FALSE,
                 control = optim_control(),
                 ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Lm"
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("LmME", class(out))
  return(out)
}


##' Extract the coefficients of an \code{LmME} model
##'
##' Extracts the fixed effects coefficents (default behavior), the baseline
##' parameters or all (baseline, fixed and random) coefficients of the model.
##'
##' See also the documentation of \code{\link{coef.tramME}}.
##' @param object An \code{LmME} object.
##' @param as.lm If \code{TRUE}, return the transformed coefficients as in a
##'   \code{lmerMod} object.
##' @param ... Optional arguments passed to \code{\link{coef.tramME}}.
##' @inheritParams coef.tramME
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' coef(fit, as.lm = TRUE)
##' @importFrom stats coef
##' @export
## FIXME: remove FE params for smooth terms?
coef.LmME <- function(object, as.lm = FALSE, fixed = TRUE, ...) {
  class(object) <- class(object)[-1L]
  if (!as.lm)
    return(coef(object, fixed = fixed, ...))

  if (!is.null(object$model$ctm$bases$interacting))
    stop("Cannot compute scaled coefficients with strata.")

  par <- coef(object, with_baseline = TRUE, fixed = fixed, ...)
  rn <- variable.names(object, "response")

  scidx <- grep(rn, names(par), fixed = TRUE)
  icidx <- which(names(par) == "(Intercept)")
  sig <- 1 / par[scidx]
  par <- c(-par[icidx], par[-c(icidx, scidx)]) * sig
  return(par)
}


##' Extract the SD of the error term of an LmME model.
##' @param object An \code{LmME} object.
##' @param ... Optional argument (for consistency with generic).
##' @return A numeric value of the transformed sigma parameter.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' sigma(fit)
##' @importFrom stats sigma
##' @export
sigma.LmME <- function(object, ...) {
  if (!is.null(object$model$ctm$bases$interacting))
    stop("Cannot compute residual standard error with strata.")
  class(object) <- class(object)[-1L]
  par <- coef(object, with_baseline = TRUE, fixed = TRUE)
  rn <- variable.names(object, "response")
  scidx <- grep(rn, names(par), fixed = TRUE)
  sig <- 1 / par[scidx]
  return(unname(sig))
}


##' Get the variance-covariance matrix of the parameters of an LmME model
##'
##' @param object A fitted \code{LmME} object.
##' @param as.lm If \code{TRUE}, return the covariance matrix of the same
##'   parameterization as used by \code{\link[lme4]{lmer}}.
##' @inheritParams confint.LmME
##' @return A numeric covariance matrix.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' vcov(fit) ## transformation model parameterization
##' vcov(fit, as.lm = TRUE) ## LMM parameterization
##' vcov(fit, as.lm = TRUE, pargroup = "fixef") ## cov of fixed effects
##' @importFrom stats vcov
##' @export
## FIXME: standradize sdreport call throught the package
vcov.LmME <- function(object, as.lm = FALSE, parm = NULL,
                      pargroup = c("all", "fixef", "ranef"), ...) {
  pargroup <- pargroup[1L]
  class(object) <- class(object)[-1L]
  if (!as.lm) {
    return(vcov(object, parm, pargroup, ...))
  }
  data <- object$tmb_obj$env$data
  data$as_lm <- as.numeric(as.lm)
  newobj <- update(object$tmb_obj, data = data,
                   parameters = .get_par(object$tmb_obj, full = TRUE))
  sdr <- TMB::sdreport(newobj, getReportCovariance = TRUE)
  ## FIXME: robustify cov matrix calculations
  vc <- sdr$cov
  pr <- names(sdr$value)
  pr <- switch(pargroup, all = pr %in% c("b", "sigma", "th"),
               fixef = pr == "b", ranef = pr == "th",
               stop("Unknown parameter group for parameterization with as.lm = TRUE."))
  nm <- c("(Intercept)",
          names(coef(object, with_baseline = TRUE, fixed = FALSE))[-(1:2)],
          "(Sigma)",
          names(varcov(object, as.theta = TRUE, full = TRUE, fixed = FALSE)))
  if (length(parm)) {
    pr <- pr & (nm %in% parm)
  }
  vc <- vc[pr, pr, drop = FALSE]
  colnames(vc) <- rownames(vc) <- nm[pr]
  return(vc)
}

##' Variances and correlation matrices of random effects of an LmME object
##'
##' The returned parameters are the transformed versions of the original parameters that
##' correspond to the normal linear mixed model parameterization.
##'
##' The function only returns the correlation matrices that belong to actual random effects
##' (defined for groups in the data) and ignores the random effects parameters of the smooth
##' shift terms. To extract these, the user should use \code{varcov} with \code{full = TRUE}.
##' @param x An \code{LmME} object.
##' @param sigma Standard deviation of the error term in the LMM parameterization (should
##'   not be set manually, only for consistency with the generic method)
##' @param as.lm If \code{TRUE}, return the variances and correlations that correspond to
##'   a normal linear mixed model (i.e. \code{lmerMod}).
##' @param ... Optional arguments (for consistency with generic)
##' @return A list of vectors with variances and correlation matrices corresponding to the
##'   various grouping variables.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' VarCorr(fit) ## tranformation model parameterization
##' VarCorr(fit, as.lm = TRUE) ## LMM parameterization
##' @importFrom nlme VarCorr
##' @importFrom stats sigma
##' @export VarCorr
##' @export
VarCorr.LmME <- function(x, sigma = 1, as.lm = FALSE, ...) {
  if (!as.lm) {
    class(x) <- class(x)[-1L]
    return(VarCorr(x))
  }
  if (missing(sigma))
    sigma <- sigma(x)
  class(x) <- class(x)[-1L]
  vc <- VarCorr(x)
  vc <- lapply(vc, function(xx) {
    xx$var <- xx$var * sigma^2
    xx
  })
  class(vc) <- c("VarCorr.tramME", class(vc))
  return(vc)
}


##' Extract the variance-covariance matrix of the random effects of an LmME model
##' @param object A \code{LmME} object.
##' @param as.lm If \code{TRUE}, the returned values correspond to the LMM
##'   parameterization.
##' @param as.theta Logical value, if \code{TRUE}, the values are returned
##'   in their reparameterized form.
##' @param full Logical value; if \code{TRUE}, return all random effects elements,
##'   if \code{FALSE}, do not return the random effects parameters of the smooth
##'   terms.
##' @param ... Optional arguments (unused).
##' @return A list of the covariance matrices or a vector of theta values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' varcov(fit, as.lm = TRUE)
##' varcov(fit, as.theta = TRUE, as.lm = TRUE)
##' @export
varcov.LmME <- function(object, as.lm = FALSE, as.theta = FALSE, full = FALSE, ...) {
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    return(varcov(object, as.theta = as.theta, full = full, ...))
  }
  sig <- sigma(object)
  vc <- varcov.tramME(object, as.theta = FALSE, full = full, ...)
  vc <- lapply(vc, function(x) x * sig^2)
  if (as.theta) {
    th <- varcov(object, as.theta = TRUE, full = full) ## NOTE: for the names
    bls <- attr(object$param, "re")$blocksize
    if (full)
      bls <- c(bls, rep(1, length(attr(object$param, "sm")$re_dims)))
    th[] <- .vc2th(vc, bls)
    return(th)
  }
  return(vc)
}


##' Confidence intervals for LmME model parameters
##'
##' Confidence intervals for model parameters on their original scale,
##' optionally consistent with the linear mixed-model specification.
##' When \code{as.lm = TRUE}, only Wald CIs are available.
##' @param object An \code{LmME} object.
##' @param parm Names of the parameters to extract.
##' @param as.lm Logical. If \code{TRUE}, return results consistent with the normal linear
##'   mixed model parameterization.
##' @param pargroup The name of the parameter group to extract. With \code{as.lm = FALSE},
##'   the available options are described in \code{confint.tramME}. When \code{as.lm = TRUE},
##'   the following options are available:
##'   \itemize{
##'     \item all: Fixed effects and variance components parameters.
##'     \item fixef: Fixed effects parameters (including FE parameters of the smooth terms).
##'     \item ranef: Variance components parameters (including the smoothing parameters of
##'       the random effects).
##'   }
##' @param ... Optional parameters passed to \code{confint.tramME}
##' @inheritParams confint.tramME
##' @return A matrix with lower and upper bounds.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' confint(fit) ## transformation model parameterization
##' confint(fit, as.lm = TRUE) ## LMM parameterization
##' confint(fit, as.lm = TRUE, pargroup = "fixef", estimate = TRUE)
##' confint(fit, as.lm = TRUE, parm = "(Sigma)") ## error SD
##' @importFrom stats confint qnorm
##' @export
confint.LmME <- function(object, parm = NULL, level = 0.95,
                         as.lm = FALSE,
                         pargroup = c("all", "fixef", "ranef"),
                         type = c("Wald", "wald", "profile"),
                         estimate = FALSE, ...) {
  type <- tolower(match.arg(type))
  pargroup <- pargroup[1L]
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    fc <- match.call()
    fc$object <- object
    fc[[1L]] <- quote(confint)
    return(eval(fc))
  }

  b <- coef(object, fixed = FALSE, as.lm = TRUE)
  sig <- sigma(object)
  th <- varcov(object, fixed = FALSE, as.theta = TRUE, full = TRUE)
  pr <- c(b, sig, th)

  vc <- vcov(object, as.lm = TRUE)
  se <- sqrt(diag(vc))
  nm <- colnames(vc)

  g <- switch(pargroup, all = rep(TRUE, length(nm)), fixef = nm %in% names(b),
              ranef = nm %in% names(th), stop("Unknown parameter group."))
  if (length(parm)) {
    g <- g & (nm %in% parm)
  }

  if (type == "wald") {
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    ci <- pr + se %o% fac
  } else {
    stop("Only Wald confidence intervals are available with as.lm = TRUE")
  }

  if (estimate) {
    ci <- cbind(ci, pr)
  }
  rownames(ci) <- nm
  colnames(ci) <- if (estimate) c("lwr", "upr", "est") else c("lwr", "upr")
  ci <- ci[g, , drop = FALSE]
  return(ci)

}

##' Extract the conditional modes of random effects of an LmME model
##'
##' The \code{condVar} option is not implemented for \code{ranef.LmME}.
##' Setting \code{raw=TURE} will return the raw random effects estimates from
##' the transformation model parameterization.
##' @param object A fitted LmME object.
##' @param as.lm If \code{TRUE}, return the transformed conditional modes as in a
##'   normal linear mixed effects model.
##' @param ... Optional parameters passed to \code{ranef.tramME}.
##' @return A numeric vector or a \code{ranef.tramME} object depending on the inputs.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' ranef(fit, raw = TRUE) ## transformation model parameterization!
##' ranef(fit, as.lm = TRUE)
##' @importFrom nlme ranef
##' @importFrom stats sigma
##' @export
## FIXME: with condVar?
## FIXME: change this when ranef.tramME changes
ranef.LmME <- function(object, as.lm = FALSE, ...) {
  if (!as.lm || isTRUE(list(...)$raw)) {
    class(object) <- class(object)[-1L]
    return(ranef(object, ...))
  }
  if (!is.null(list(...)$param))
    stop("Setting the param argument is not supported with as.lm = TRUE.")
  if (isTRUE(list(...)$condVar))
    warning("condVar option is not available with as.lm = TRUE")
  sig <- sigma(object)
  class(object) <- class(object)[-1L]
  re <- ranef(object, condVar = FALSE, ...)
  out <- lapply(re, function(x) x * sig)
  return(out)
}

##' Residuals of a LmME model
##'
##' Calculates the score residuals of an intercept term fixed at 0.
##' In the case of an LmME model, this is equal to the residual of an LMM.
##' @param object An \code{LmME} object.
##' @param as.lm If \code{TRUE}, return the residuals as in a normal linear
##'   mixed effects model.
##' @inheritParams residuals.tramME
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' resid(fit)
##' @importFrom stats residuals
##' @export
residuals.LmME <- function(object, as.lm = FALSE, ...) {
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    return(residuals(object, ...))
  }
  if (!is.null(list(...)$param))
    stop("Setting the param argument is not supported with as.lm = TRUE.")
  sig <- sigma(object)
  class(object) <- class(object)[-1L]
  res <- residuals(object, ...)
  return(res * sig)
}
