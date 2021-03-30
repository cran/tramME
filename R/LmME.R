##' ME version of tram::Lm
##' @inheritParams tram::Lm
##' @param silent Logical. Make \pkg{TMB} functionality silent.
##' @param resid Logical. If \code{TRUE}, the score residuals are also calculated.
##'   This comes with some performance cost.
##' @param do_update Logical. If \code{TRUE}, the model is set up so that the weights and the
##'   offsets are updateable. This comes with some performance cost.
##' @param estinit logical, estimate a vector of initial values for the fixed effects parameters
##'   from a (fixed effects only) mlt model
##' @param initpar named list of initial parameter values, if \code{NULL}, it is ignored
##' @inheritParams mlt::mlt
##' @param nofit logical, if TRUE, creates the model object, but does not run the optimization
##' @param control list with controls for optimization
##' @return A LmME object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Lm
##' @export
LmME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                 silent = TRUE, resid = FALSE, do_update = FALSE,
                 estinit = TRUE, initpar = NULL,
                 fixed = NULL, nofit = FALSE,
                 control = optim_control(),
                 ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Lm"
  mod <- eval(fc, parent.frame())

  ## -- sanitize initial parameter settings
  if (is.null(mod$ranef) || nofit || !is.null(initpar)) {
    estinit <- FALSE
  }

  ## -- create model frame
  if (missing(data) || !inherits(data, "tramME_data")) {
    fc <- cl
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"),
               names(fc), 0L)
    fc <- fc[c(1L, m)]
    fc$formula <- .combine_formulas(mod$formula)
    fc[[1L]] <- quote(tram::tram_data)
    out <- eval(fc, parent.frame())
    dat <- out$mf
    class(dat) <- c("tramME_data", class(dat))
  } else {
    dat <- data
  }

  cf <- coef(mod$ctm)

  ## -- create terms required by tramTMB
  mmlt <- mlt::mlt(mod$ctm, data = dat, offset = model.offset(dat),
              weights = model.weights(dat), ## TODO: offset and weights might not be needed
              fixed = fixed, dofit = estinit)
  fe <- fe_terms(mmlt)
  re <- re_terms(mod$ranef, dat, mod$negative)
  inp <- tramTMB_inputs(mod, fe, re, dat, param = initpar)

  cf[names(fixed)] <- fixed

  mp <- list()
  if (!is.null(fixed)) {
    idx <- which(names(cf) %in% names(fixed))
    bb <- rep(NA, length(cf))
    bb[-idx] <- seq_along(bb[-idx])
    mp <- list(beta = as.factor(bb))
  }

  ## -- create the tramTMB object
  obj <- tramTMB(inp$data, inp$parameters, inp$constraint, inp$negative,
                 map = mp, resid = resid, do_update = do_update, silent = silent)

  ## -- model fitting
  if (!nofit) {

    if (is.null(initpar) && !estinit) {
      par <- .optim_start(obj, resp = dat[[1]])
    } else par <- NULL

    opt <- optim_tramTMB(obj, par = par,
                         method = control$method, control = control$control,
                         trace = control$trace, ntry = control$ntry,
                         scale = control$scale)
    parm <- .get_par(obj)
  } else {
    opt <- NULL
    parm <- list(beta = cf, theta = rep(NA, length(.get_par(obj)$theta)))
  }

  param <- .gen_param(parm, fe = list(names = names(cf)),
                      re = list(names = re$names, blocksize = re$blocksize,
                                levels = re$levels, termsize = re$termsize),
                      varnames = names(dat))
  structure(list(call = cl, model = mod, data = dat, tmb_obj = obj, opt = opt,
                 param = param),
            class = c("LmME", "tramME"))
}


##' Extract the coefficients of the fixed effects terms of an LmME model.
##' @param object An \code{LmME} object.
##' @param as.lm If \code{TRUE}, return the transformed coefficients as in a
##'   \code{lmerMod} object.
##' @param ... Optional arguments passed to \code{coef.tramME}.
##' @inheritParams coef.tramME
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' coef(fit, as.lm = TRUE)
##' @importFrom stats coef
##' @export
coef.LmME <- function(object, as.lm = FALSE, fixed = TRUE, ...) {
  class(object) <- class(object)[-1L]
  if (!as.lm)
    return(coef(object, fixed = fixed, ...))

  if (!is.null(object$model$ctm$bases$interacting))
    stop("Cannot compute scaled coefficients with strata.")

  par <- coef(object, with_baseline = TRUE, fixed = fixed)
  rn <- variable.names(object, "response")

  scidx <- grep(rn, names(par), fixed = TRUE)
  icidx <- grep("(Intercept)", names(par), fixed = TRUE)
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
  class(object) <- class(object)[-1L]
  sig <- 1 / par[scidx]
  return(unname(sig))
}


##' Get the variance-covariance matrix of the parameters of an LmME model
##'
##' \code{pargroup = "baseline"} with the option \code{as.survreg = TRUE} is
##' not available for \code{LmME} objects.
##' @param object A fitted \code{LmME} object.
##' @param as.lm If \code{TRUE}, return the covariance matrix of the same
##'   parametrization as used by \code{\link[lme4]{lmer}}.
##' @inheritParams vcov.tramME
##' @return A numeric covariance matrix.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' vcov(fit) ## transformation model parametrization
##' vcov(fit, as.lm = TRUE) ## LMM parametrization
##' ## cov of coefficient AND other terms with 'Days' in names
##' vcov(fit, as.lm = TRUE, parm = "Days", pmatch = TRUE)
##' vcov(fit, as.lm = TRUE, parm = "^Days", pmatch = TRUE) ## var of coefficient only
##' vcov(fit, as.lm = TRUE, pargroup = "fixef") ## cov of fixed effects
##' @importFrom stats vcov
##' @export
vcov.LmME <- function(object, as.lm = FALSE, parm = NULL,
                      pargroup = c("all", "fixef", "shift", "baseline", "ranef"),
                      pmatch = FALSE, ...) {
  pargroup <- match.arg(pargroup)
  class(object) <- class(object)[-1L]
  if (!as.lm) {
    return(vcov(object, parm, pargroup, pmatch, ...))
  }
  if (!is.null(object$model$ctm$bases$interacting))
    stop("Covariance matrix for scaled parameters with strata.")
  stopifnot(pargroup != "baseline")
  vc <- vcov(object, pargroup = "all", ...)

  b <- .get_cf(object)[.idx(object, fixed = FALSE, pargroup = "fixef")]
  bs <- .get_cf(object)[.idx(object, fixed = FALSE, pargroup = "shift")]
  th <- .get_vc(object, as.theta = TRUE)

  pr <- c(b, th)
  vn <- variable.names(object, "response")
  sidx <- grep(vn, names(pr), fixed = TRUE)
  iidx <- grep("^\\(Intercept\\)$", names(pr))

  bls <- attr(object$param, "re")$blocksize
  np <- length(.idx(object, fixed = FALSE, pargroup = "ranef"))

  ic <- pr[iidx]
  sig <- 1 / pr[sidx]

  ## Construct the Jacobian block-by-block
  bl1 <- Matrix::Matrix(0, nrow = length(b), ncol = length(pr))
  bl1[iidx, iidx] <- -sig
  bl1[iidx, sidx] <- ic * sig^2
  bl1[-iidx ,sidx] <- -bs * sig^2
  bl1[-c(iidx, length(b)), -c(iidx, sidx)] <-
    Matrix::Matrix(diag(sig, nrow = length(bs), ncol = length(pr)-2))
  bl1[length(b), sidx] <- -sig^2
  if (np > 0) {
    bl2 <- Matrix::Matrix(0, nrow = np, ncol = length(b))
    bl2[, sidx] <- sapply(bls, function(x){
      v <- rep(0, x * (x+1) / 2)
      v[1:x] <- -sig
      v
    })
    bl3 <- Matrix::Matrix(diag(np))
    ja <- rbind(bl1, cbind(bl2, bl3))
  } else {
    ja <- bl1
  }

  ## Delta method
  out <- ja %*% vc %*% Matrix::t(ja)
  idx <- .idx(object, pargroup = pargroup, which = parm, pmatch = pmatch, altpar = "lm")
  out <- as.matrix(out[idx, idx])
  colnames(out) <- names(idx)
  rownames(out) <- names(idx)
  return(out)
}


##' Variances and correlation matrices of random effects of an LmME object
##'
##' The returned parameters are the transformed versions of the original parameters that
##' correspond to the normal linear mixed model parametrization.
##' @param x An \code{LmME} object.
##' @param sigma Standard deviation of the error term in the LMM parametrization (should
##'   not be set manually, only for consistency with the generic method)
##' @param as.lm If \code{TRUE}, return the variances and correlations that correspond to
##'   a normal linear mixed model (i.e. \code{lmerMod}).
##' @param ... Optional arguments (for consistency with generic)
##' @return A list of vectors with variances and correlation matrices corresponding to the
##'   various grouping variables.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' VarCorr(fit) ## tranformation model parametrization
##' VarCorr(fit, as.lm = TRUE) ## LMM parametrization
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
##'   parametrization.
##' @param as.theta Logical value, if \code{TRUE}, the values are returned
##'   in their reparameterized form.
##' @param ... Optional arguments (unused).
##' @return A list of the covariance matrices or a vector of theta values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' varcov(fit, as.lm = TRUE)
##' varcov(fit, as.theta = TRUE, as.lm = TRUE)
##' @export
varcov.LmME <- function(object, as.lm = FALSE, as.theta = FALSE, ...) {
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    return(varcov(object, as.theta = as.theta, ...))
  }
  sig <- sigma(object)
  vc <- varcov.tramME(object, as.theta = FALSE, ...)
  vc <- lapply(vc, function(x) x * sig^2)
  if (as.theta) {
    th <- varcov(object, as.theta = TRUE) ## NOTE: for the names
    th[] <- .vc2th(vc, attr(object$param, "re")$blocksize)
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
##' @param as.lm Logical. If \code{TRUE}, return results consistent with the normal linear
##'   mixed model parametrization.
##' @param ... Optional parameters passed to \code{confint.tramME}
##' @inheritParams confint.tramME
##' @return A matrix with lower and upper bounds.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' confint(fit) ## transformation model parametrization
##' confint(fit, as.lm = TRUE) ## LMM parametrization
##' confint(fit, as.lm = TRUE, pargroup = "fixef", estimate = TRUE)
##' confint(fit, as.lm = TRUE, parm = "(Sigma)") ## error SD
##' @importFrom stats confint qnorm
##' @export
confint.LmME <- function(object, parm = NULL, level = 0.95,
                         as.lm = FALSE,
                         pargroup = c("all", "fixef", "shift", "baseline", "ranef"),
                         type = c("Wald", "wald", "profile"),
                         estimate = FALSE,
                         pmatch = FALSE, ...) {
  type <- tolower(match.arg(type))
  pargroup <- match.arg(pargroup)
  if (!as.lm) {
    fc <- match.call()
    fc[[1L]] <- quote(confint.tramME)
    return(eval(fc))
  }

  b <- coef(object, fixed = FALSE, with_baseline = TRUE, as.lm = TRUE)
  th <- varcov(object, as.theta = TRUE, as.lm = TRUE)
  pr <- c(b, th)

  idx <- .idx(object, pargroup = pargroup, which = parm, pmatch = pmatch, altpar = "lm")
  par <- pr[idx]

  if (any(is.na(pr)) || length(par) == 0) {
    nc <- if (estimate) 3 else 2
    ci <- matrix(NA, nrow = length(par), ncol = nc)
    rownames(ci) <- names(par)
    colnames(ci) <- if (nc == 3) c("lwr", "upr", "est") else c("lwr", "upr")
    return(ci)
  }

  if (type == "wald") {
    ses <- sqrt(diag(vcov(object, as.lm = TRUE)))[idx]
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    ci <- par + ses %o% fac
  } else {
    stop("Only Wald confidence intervals are available with as.lm = TRUE")
  }
  colnames(ci) <- c("lwr", "upr")
  rownames(ci) <- names(idx)
  if (estimate) {
    ci <- cbind(ci, par)
    colnames(ci) <- c("lwr", "upr", "est")
  }
  return(ci)
}


##' Extract the conditional modes of random effects of an LmME model
##'
##' The \code{condVar} option is not implemented for \code{ranef.LmME}.
##' Setting \code{raw=TURE} will return the raw random effects estimates from
##' the transformation model parametrization.
##' @param object A fitted LmME object.
##' @param as.lm If \code{TRUE}, return the transformed conditional modes as in a
##'   normal linear mixed effects model.
##' @param ... Optional parameters passed to \code{ranef.tramME}.
##' @return A numeric vector or a \code{ranef.tramME} object depending on the inputs.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' ranef(fit, raw = TRUE) ## transformation model parametrization!
##' ranef(fit, as.lm = TRUE)
##' @importFrom nlme ranef
##' @importFrom stats sigma
##' @export
## FIXME: with condVar?
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
