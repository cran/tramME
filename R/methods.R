##' Set coefficients of a tramME model.
##'
##' Sets the whole vector of coefficients of a tramME model. The parameters of
##' the baseline transformation function should respect the restrictions of
##' the parameter space. This is checked before setting the new parameter values.
##' When called on a fitted tram object, the function sets it to unfitted and removes
##' all parts that come from the estimation.
##' @param object A tramME object (fitted or unfitted).
##' @param value Numeric vector of new coefficient values.
##' @return An unfitted tramME object with the new coefficient values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' coef(mod) <- c(-1, 0.5, 1)
##' @importFrom mlt "coef<-"
##' @export
"coef<-.tramME" <- function(object, value) {
  stopifnot(.check_coef(value, object$model))
  if (object$fitted) {
    warning("The model object is fitted. Setting it to unfitted.")
    object$fitted <- FALSE
    object$tmb_obj <- object$data <- object$opt <- object$tmb_sdr <- NULL
  }
  object$pars$coef[] <- value
  return(object)
}


##' Generic method for \code{"varcov<-"}
##' @param object A model object
##' @param value The new value of the covariance matrix
##' @return An object with the same class as \code{object}, with updated
##'   variance-covariance matrix.
##' @export
"varcov<-" <- function(object, value)
  UseMethod("varcov<-")


##' Set the values of the random effects covariance matrices of a tramME model.
##'
##' Sets the list containing the covariance matrices of a tramME model. The matrices have
##' to be positive definite. Just as in \code{"coef<-"}, when the function is called
##' on a fitted object, it will be set to unfitted.
##'
##' The supplied list does not have to be named, and the names will be ignored.
##' When multiple grouping factors are present, the function assumes the same order as in the
##' object to be modified. Hence, it might be a good idea to call \code{varcov} first, and
##' modify this list to make sure that the input has the right structure.
##' @param object A tramME object (fitted or unfitted).
##' @param value A list of positive definite covariance matrices.
##' @return An unfitted tramME object with the new coefficient values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' vc <- varcov(mod)
##' vc[[1]] <- matrix(c(1, 0, 0, 2), ncol = 2)
##' varcov(mod) <- vc
##' @export
"varcov<-.tramME" <- function(object, value) {
  vc <- object$pars$varcov
  stopifnot(.check_varcov(value, vc))
  if (object$fitted) {
    warning("The model object is fitted. Setting it to unfitted.")
    object$fitted <- FALSE
    object$tmb_obj <- object$data <- object$opt <- object$tmb_sdr <- NULL
  }
  vc <- mapply(function(m1, m2) { m1[] <- m2; m1 }, m1 = vc, m2 = value,
               SIMPLIFY = FALSE)
  object$pars$varcov <- vc
  return(object)
}


##' Generic method for \code{varcov}
##' @param object A model object
##' @param ... Optional parameters
##' @return A variance-covariance matrix.
##' @export
varcov <- function(object, ...)
  UseMethod("varcov")


##' Extract the variance-covariance matrix of the random effects
##'
##' Returns the covariance matrix of the random effects as saved in the tramME object.
##' The returned values correspond to the transformation model parametrization.
##' @param object A tramME object (fitted or unfitted).
##' @param ... Optional arguments (unused)
##' @return A list of the covariance matrices.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' varcov(fit)
##' @export
varcov.tramME <- function(object, ...) {
  object$pars$varcov
}


##' Extract the coefficients of the fixed effects terms.
##' @param object A tramME object (fitted or unfitted)
##' @param with_baseline If TRUE, include the baseline parameters, too.
##' @param ... Optional parameters (ignored).
##' @return Numeric vector of parameter values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' coef(mod, with_baseline = TRUE)
##' @importFrom stats coef
##' @export
coef.tramME <- function(object, with_baseline = FALSE, ...) {
  cf <- object$pars$coef
  if (with_baseline) {
    return(cf[.paridx(object$model, pargroup = "fixef")])
  } else {
    return(cf[.paridx(object$model, pargroup = "shift")])
  }
}


##' Extract the coefficients of the fixed effects terms of an LmME model.
##' @param object An LmME object (fitted or unfitted).
##' @param as.lm If TRUE, return the transformed coefficients as in a
##'   \code{lmerMod} object.
##' @param ... optional parameters passed to \code{coef.tramME}
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' coef(fit, as.lm = TRUE)
##' @importFrom stats coef
##' @export
coef.LmME <- function(object, as.lm = FALSE, ...) {
  class(object) <- class(object)[-1L]
  if (!as.lm) return(coef(object, ...))
  if (!object$fitted) { ## Manually transform
    par <- object$pars$coef
    sig <- 1 / par[2L]
    par <- c(-par[1L], par[-(1:2)]) * sig
  } else { ## get from the tmb report
    par <- object$tmb_sdr$value
    par <- par[1:(object$model$fixef$npar["all"]-1)]
  }
  names(par) <- c("(Intercept)", object$model$fixef$names)
  return(par)
}


##' Extract the SD of the error term of an LmME model.
##' @param object An LmME object (fitted or unfitted).
##' @param ... Optional argument (for consistency with generic)
##' @return A numeric value of the transformed sigma parameter.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' sigma(fit)
##' @importFrom stats sigma
##' @export
sigma.LmME <- function(object, ...) {
  if (!object$fitted)
    return(unname(1 / object$pars$coef[2L]))
  sigma <- object$tmb_sdr$value
  sigma <- unname(sigma[object$model$fixef$npar["all"]])
  return(sigma)
}


##' Get the log-likelihood of the model
##' @param object A fitted tramME model
##' @param ... Optional argument (for consistency with generic)
##' @return A numeric value of the log-likelihood at its optimum.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' logLik(fit)
##' @importFrom stats logLik
##' @export
logLik.tramME <- function(object, ...) {
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  ll <- -object$opt$value
  df <- object$model$fixef$npar["all"] + object$model$ranef$npar
  nobs <- sum(object$data$weights)
  out <- structure(ll, df = df, nobs = nobs)
  class(out) <- "logLik"
  return(out)
}


##' Comparison of nested tramME models.
##'
##' Calculates information criteria and LR ratio test for nested tramME models.
##' The calculation of the degrees of freedom is problematic, because the
##' parameter space is restricted.
##'
##' Currently only supports the comparison of two models. Additional arguments
##' will be ignored.
##'
##' The nestedness of the models is not checked.
##' @param object A fitted tramME model.
##' @param object2 A fitted tramME model.
##' @param ... Optional arguments, for compatibility with the generic. (Ignored)
##' @return A data.frame with the calculated statistics.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod1 <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' mod2 <- LmME(Reaction ~ Days + (Days || Subject), data = sleepstudy)
##' anova(mod1, mod2)
##' @importFrom stats anova pchisq
##' @export
anova.tramME <- function(object, object2, ...) {
  stopifnot(inherits(object2, "tramME"))
  ll_  <- lapply(list(object, object2), logLik)
  stopifnot(attr(ll_[[1]], "nobs") == attr(ll_[[2]], "nobs"))
  out <- data.frame(Df = sapply(ll_, attr, "df"), logLik = unlist(ll_),
                    AIC = sapply(list(object, object2), "AIC"),
                    BIC = sapply(list(object, object2), "BIC"))
  rownames(out) <- c("Model 1", "Model 2")
  ord <- order(out$Df)
  out <- out[ord, ]
  out$Chisq <- 2 * pmax(0, c(NA, diff(out$logLik)))
  out$`Chisq df` <- c(NA, diff(out$Df))
  out$`Pr(>Chisq)` <- pchisq(out$Chisq, df = out$`Chisq df`,lower.tail = FALSE)
  class(out) <- c("anova.tramME", class(out))
  attr(out, "title") <- "Model comparison"
  attr(out, "ordering") <- ord ## ordered by complexity
  attr(out, "models") <-
    sapply(list(object, object2), function(x) x$call$formula)## in original order
  return(out)
}


##' Printing \code{anova.tramME} table
##' @param x A \code{anova.tramME} object.
##' @param ... Optional arguments passed to \code{\link[stats]{printCoefmat}}
##' @return Invisibly retrurns the \code{anova.tramME} object.
##' @inheritParams stats::printCoefmat
##' @importFrom stats printCoefmat
##' @export
print.anova.tramME <- function(x, digits = max(getOption("digits") - 2L, 3L),
                              signif.stars = getOption("show.signif.stars"), ...) {
  cat(attr(x, "title"), "\n\n", sep = "")
  cat(paste0("\t Model ", 1:nrow(x), ": ", attr(x, "models"), "\n"), "\n", sep = "")
  printCoefmat(x, digits = digits, signif.stars = signif.stars,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = NULL,
               zap.ind = 1L:ncol(x), tst.ind = 5L, na.print = "", ...)
  return(invisible(x))
}


##' Extract the conditional modes and conditional variances of random effects
##'
##' @param object A fitted tramME object.
##' @param condVar If TRUE, include the conditional variances as attributes.
##' @param raw Return the unformatted RE estimates as fitted by the model.
##' @param ... Optional arguments (for consistency with generic)
##' @return Depending on the value of raw, either a numeric vector or a
##'   \code{ranef.tramME} object which contains the conditional mode and variance
##'   estimates by grouping factors.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy, order = 5)
##' ranef(fit, raw = TRUE)
##' ranef(fit)
##' @importFrom nlme ranef
##' @aliases ranef
##' @export ranef
##' @export
ranef.tramME <- function(object, condVar = FALSE, raw = FALSE, ...) {
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  par <- unname(object$tmb_sdr$par.random) ## NOTE: alternatively, using tmb_object$env$parList(...)
  if (raw) return(par)
  out <- .re_format(object$model$ranef, par)
  if (condVar) {
    cv <- object$tmb_sdr$diag.cov.random ## NOTE: alternatively from summary(object$tmb_sdr, "random")
    cv <- .re_format(object$model$ranef, cv)
    out <- mapply(FUN = function(cm, cv) {
      attr(cm, "condVar") <- cv
      cm
    }, cm = out, cv = cv, SIMPLIFY = FALSE)
  }
  class(out) <- c("ranef.tramME", class(out))
  return(out)
}


##' Extract the conditional modes of random effects of an LmME model
##'
##' The \code{condVar} option is not implemented for \code{ranef.LmME}.
##' Setting \code{raw=TURE} will return the raw random effects estimates from
##' the transformation model parametrization.
##' @param object A fitted LmME object.
##' @param as.lm If TRUE, return the transformed conditional modes as in a
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
## TODO: with condVar?
ranef.LmME <- function(object, as.lm = FALSE, ...) {
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    return(ranef(object, ...))
  }
  if (isTRUE(list(...)$condVar))
    warning("condVar option is not available with as.lm = TRUE")
  sig <- sigma(object)
  class(object) <- class(object)[-1L]
  re <- ranef(object, condVar = FALSE)
  out <- lapply(re, function(x) x * sig)
  return(out)
}


##' Calculate the variance-covariance matrix of the parameters
##'
##' Extracts the covariance matrix of the selected parameters. The returned values
##' are on the same scale as the estimated parameter values, i.e. the standard
##' deviations of the random effect terms are on log scale.
##'
##' The argument \code{parm} defines the indices or the names of the parameters
##' of interest within the selected \code{pargroup}. When \code{pmatch = TRUE},
##' partial matching of parameter names is allowed.
##' @param object A fitted tramME object.
##' @param parm The indeces or names of the parameters of interest. See in details.
##' @param pargroup fixef: fixed-effects, shift: shift parameters, all: fixed
##'   effects and variance component parameters, baseline: parameters of the
##'   baseline transformation function, ranef: variance components parameters.
##' @param pmatch Logical. If \code{TRUE}, partial name matching is allowed.
##' @param ... Optional arguments
##' @return A numeric covariance matrix.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy, order = 10)
##' vcov(fit)
##' vcov(fit, pargroup = "ranef")
##' vcov(fit, pargroup = "baseline")
##' vcov(fit, parm = "Reaction") ## same as previous
##' @importFrom stats vcov
##' @export
vcov.tramME <- function(object, parm = NULL,
                       pargroup = c("all", "fixef", "shift", "baseline", "ranef"),
                       pmatch = FALSE, ...) {
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  pargroup <- match.arg(pargroup)
  out <- object$tmb_sdr$cov.fixed
  idx <- .paridx(object$model, parm, pargroup, pmatch)
  out <- as.matrix(out[idx, idx])
  colnames(out) <- names(idx)
  rownames(out) <- names(idx)
  return(out)
}


##' Get the variance-covariance matrix of the parameters of an LmME model
##'
##' \code{pargroup = "baseline"} is not available for \code{LmME} objects.
##' @param object A fitted LmME object.
##' @param as.lm If TRUE, return the covariance matrix of the transformed
##'   parameters as in a \code{lmerMod} object.
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
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  pargroup <- match.arg(pargroup)
  if (!as.lm) {
    class(object) <- class(object)[-1L]
    return(vcov(object, parm, pargroup, pmatch))
  }
  out <- object$tmb_sdr$cov
  idx <- .paridx(object$model, parm, pargroup, pmatch, altpar = "lm")
  out <- as.matrix(out[idx, idx])
  colnames(out) <- names(idx)
  rownames(out) <- names(idx)
  return(out)
}


##' Variances and correlation matrices of random effects
##'
##' This function calculates the variances and correlations from \code{varcov.tramME}.
##' @param x A tramME object
##' @param ... optional arguments (for consistency with the generic method)
##' @return A list of vectors with variances and correlation matrices corresponding to the
##'   various grouping variables.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' VarCorr(fit)
##' @importFrom nlme VarCorr
##' @aliases VarCorr
##' @export VarCorr
##' @export
VarCorr.tramME <- function(x, ...) {
  vc <- varcov(x)
  nl <- sapply(x$model$ranef$levels, length)
  out <- mapply(function(xx, n) {
    nm <- rownames(xx)
    v <- diag(xx)
    h <- diag(1/sqrt(v), ncol = length(v), nrow = length(v))
    c <- h %*% xx %*% h
    names(v) <- rownames(c) <- colnames(c) <- nm
    out <- list(var = v, corr = c)
    attr(out, "nlevels") <- nl[n]
    out
  }, xx = vc, n = names(vc), SIMPLIFY = FALSE)
  names(out) <- names(vc)
  class(out) <- c("VarCorr.tramME", class(out))
  return(out)
}


##' Print method for the variance-correlation parameters of a tramME object
##' @param x A VarCorr.tramME object
##' @param sd Logical. Print standard deviations instead of variances.
##' @param digits Number of digits
##' @param ... optional arguments
##' @return Invisibly returns the input VarCorr.tramME object.
##' @export
print.VarCorr.tramME <- function(x, sd = TRUE,
                                digits = max(getOption("digits") - 2L, 3L),
                                ...) {
  for (i in seq_along(x)) {
    cat("\nGrouping factor: ", names(x)[i], " (", attr(x[[i]], "nlevels"),
        " levels)", sep = "")
    if (sd) {
      cat("\nStandard deviation:\n")
      vs <- sqrt(x[[i]]$var)
    } else {
      cat("\nVariance:\n")
      vs <- x[[i]]$var
    }
    print(signif(vs, digits))
    cr <- x[[i]]$corr
    if (nrow(cr) > 1) {
      cat("\nCorrelations:\n")
      pcr <- format(cr, digits = digits, justify = "right")
      pcr[upper.tri(pcr, diag = TRUE)] <- ""
      pcr <- pcr[-1, -ncol(pcr), drop = FALSE]
      print(noquote(pcr, right = TRUE))
    }
  }
  cat("\n")
  invisible(x)
}


##' Variances and correlation matrices of random effects of an LmME object
##'
##' The returned parameters are the transformed versions of the original parameters,
##' and correspond to the normal linear mixed model parametrization.
##' @param x An LmME object.
##' @param sigma Standard deviation of the error term in the LMM parametrization (should
##'   not be set manually, only for consistency with the generic method)
##' @param as.lm If TRUE, return the variances and correlations that correspond to
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


##' Confidence intervals for tramME model parameters
##'
##' Confidence intervals for model parameters on their original scale.
##' Either Wald CI or profile CI by root finding. Multicore computations
##' are supported in the case of profile confidence intervals, but snow
##' support is yet to be implemented.
##' @param object A fitted tramME object.
##' @param level Confidence level.
##' @param type Type of the CI: either Wald or profile.
##' @param estimate Logical, add the point estimates in a thrid column
##' @param parallel Method for parallel computation
##' @param ncpus Number of cores to use for parallel computation
##' @param ... Optional parameters
##' @inheritParams vcov.tramME
##' @return A matrix with lower and upper bounds.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' confint(fit)
##' confint(fit, pargroup = "shift", estimate = TRUE)
##' exp(confint(fit, 1:2, pargroup = "ranef")) ## CIs for the SDs of the REs
##' @importFrom stats confint qnorm qchisq
##' @export
confint.tramME <- function(object, parm = NULL,
                          level = 0.95,
                          pargroup = c("all", "fixef", "shift", "baseline", "ranef"),
                          type = c("Wald", "wald", "profile"),
                          estimate = FALSE,
                          pmatch = FALSE,
                          parallel = c("no", "multicore", "snow"),
                          ncpus = getOption("profile.ncpus", 1L), ...) {
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  type <- tolower(match.arg(type))
  pargroup <- match.arg(pargroup)
  plist <- .parallel_default(parallel, ncpus)

  ## --- Indices, point estimates
  idx <- .paridx(object$model, parm, pargroup, pmatch)
  par <- object$opt$par[idx] ## NOTE: alternative?

  if (type == "wald") { ## --- Wald CI
    ses <- sqrt(diag(vcov(object)))[idx]
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    fac <- qnorm(a)
    ci <- par + ses %o% fac

  } else if (type == "profile") { ## --- Profile CI

    fun <- function(x) {
      TMB::tmbroot(object$tmb_obj, x, target = 0.5 * qchisq(level, df = 1), ...)
    }

    if (plist$do_parallel) {
      if (plist$parallel == "multicore") {
        ci <- parallel::mclapply(idx, fun, mc.cores = ncpus)
      } else if (plist$parallel == "snow") {
        ## TODO: add snow support
        stop("No snow support yet")
      }
    } else {
      ci <- lapply(idx, fun)
    }
    ci <- do.call("rbind", ci)

  } else { ## ---
    stop("Only Wald or profile types are allowed.")
  }

  colnames(ci) <- c("lwr", "upr")
  rownames(ci) <- names(idx)
  if (estimate) {
    ci <- cbind(ci, par)
    colnames(ci) <- c("lwr", "upr", "est")
  }
  return(ci)
}


##' Confidence intervals for LmME model parameters
##'
##' Confidence intervals for model parameters on their original scale,
##' optionally consistent with the linear mixed-model specification.
##' When \code{as.lm = TRUE}, only Wald CIs are available.
##' @param object A fitted LmME object.
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
  if (!object$fitted)
    stop("Not a fitted tramME model!")
  type <- tolower(match.arg(type))
  pargroup <- match.arg(pargroup)
  if (!as.lm) {
    class(object) <- class(object)[-1L] ## TODO: with match.call to make it less error-prone
    return(confint(object, parm, level, pargroup, type, estimate, pmatch, ...))
  }
  idx <- .paridx(object$model, parm, pargroup, pmatch, altpar = "lm")
  par <- object$tmb_sdr$value[idx]
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


##' Return variable names.
##'
##' Returns the variable names corresponding the selected group.
##' The returned names are derived names as tramME uses them. For example, when the
##' response is a Surv object, \code{variable.names} returns the name of that object, and
##' the names of the variables used to create it.
##' @param object a tramME object (fitted or unfitted)
##' @param which
##'   all: all non-eliminated variable names,
##'   response: response variable,
##'   grouping: grouping factors for random effects,
##'   shifting: shifting variables,
##'   interacting: interacting variables.
##' @param ... optional parameters
##' @return A vector of variable names.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' variable.names(mod)
##' variable.names(mod, "response")
##' @importFrom stats variable.names
##' @export
variable.names.tramME <- function(object,
    which = c("all", "response", "grouping", "shifting", "interacting"), ...) {
  which <- match.arg(which)
  rn <- object$model$response$name
  fen <- object$model$fixef$names
  ren <- names(object$model$ranef$names)
  out <- switch(which,
                all = c(rn, fen, ren),
                response = rn,
                grouping = ren,
                shifting = unname(variable.names(object$model$fixef$bases$shifting)),
                interacting = unname(variable.names(object$model$fixef$bases$interacting)),
                stop("Unknown variable group."))
  return(out)
}


##' Print tramME model
##' @param x A fitted or unfitted tramME model
##' @param digits Number of significant digits
##' @param ... Optional arguments (for consistency with the generic)
##' @return The original tramME obejct invisibly
##' @export
print.tramME <- function(x, digits = max(getOption("digits") - 2L, 3L), ...) {
  mnm <- .model_name(x$model)
  cat("\n", mnm, "\n", sep = "")
  cat("\n\tFormula: ")
  print(x$call$formula)
  if (x$fitted) {
    cat("\n\tFitted to dataset ")
    print(x$call$data)
  } else {
    cat("\n\tNot fitted\n")
  }
  fe <- coef(x)
  if (x$fitted && length(fe) > 0) {
    cat("\n\tEstimated fixed effects:\n")
    cat("\t========================\n\n")
    print(signif(fe, digits))
  } else if (!x$fitted && all(!is.na(fe)) && length(fe) > 0) {
    cat("\n\tFixed effects values set to:\n")
    cat("\t============================\n\n")
    print(signif(fe, digits))
  }
  vc <- VarCorr(x)
  if (x$fitted) {
    cat("\n\tEstimated random effects parameters:\n")
    cat("\t====================================\n")
    print(vc, digits  = digits)
  } else if (all(!is.na(unlist(vc)))) {
    cat("\n\tRandom effects parameters set to:\n")
    cat("\t=================================\n")
    print(vc, digits  = digits)
  }
  if (x$fitted) {
    ll <- logLik(x)
    cat("\n\tLog-likelihood: ", round(ll, digits),
        " (df = ", attr(ll, "df"), ")", sep ="")
  }
  cat("\n\n")
  invisible(x)
}


##' Summary method for tramME model
##'
##' @param object A tramME object
##' @param ... Optional arguments (for consistency with the generic)
##' @return A summary.tramME object.
##' @importFrom stats pnorm
##' @export
summary.tramME <- function(object, ...) {
  np <- object$model$fixef$npar["shift"]
  if (object$fitted) {
    ll <- logLik(object)
    se <- sqrt(diag(vcov(object, pargroup = "shift")))
  } else {
    ll <- NA
    se <- rep(NA, np)
  }
  b <- coef(object)
  zval <- b / se
  coef <- cbind(Estimate = b, `Std. Error` = se,
    `z value` = zval,
    `Pr(>|z|)` = 2 * pnorm(abs(zval), lower.tail = FALSE))
  rownames(coef) <- names(b)
  structure(
    list(name    = .model_name(object$model),
         formula = object$call$formula,
         wtd     = if (any(object$data$weights != 1)) TRUE else FALSE,
         fitted  = object$fitted,
         data    = object$call$data,
         conv    = object$opt$convergence == 0,
         nwarn   = length(object$opt$warnings),
         coef    = coef,
         varcorr = VarCorr(object),
         ll      = ll),
    class = "summary.tramME")
}


##' Print method for tramME model summary
##'
##' @param x a \code{summary.tramME} object
##' @param ... Optional arguments passed to \code{\link[stats]{printCoefmat}}
##' @return The input summary.tramME object, invisibly.
##' @inheritParams stats::printCoefmat
##' @export
print.summary.tramME <- function(x, digits = max(getOption("digits") - 2L, 3L),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("\n", x$name, "\n", sep = "")
  cat("\n\tFormula: ")
  print(x$formula)
  if (x$fitted) {
    wmsg <- if (x$wtd) " (weighted estimation)" else ""
    cat("\n\tFitted to dataset ", paste0("\033[0;32m", x$data, "\033[0m"),
        wmsg, "\n", sep = "")
    if (!x$conv)
      cat("\t\033[0;31mOptimizer did not achieve convergence!\033[0m\n") ## TODO: add later optimizer name
    if (x$nwarn > 0)
      cat("\tThere were", x$nwarn, "warning messages captured during optimization.",
          "\n", sep = " ")
  } else {
    cat("\n\t\033[0;31mNot fitted\033[0m\n")
  }
  cat("\n\tFixed effects:\n")
  cat("\t==============\n\n")
  if (nrow(x$coef) == 0) {
    cat("No estimated shift coefficients.\n")
  } else {
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1L:2L,
               tst.ind = 3L, na.print = "NA", ...)
  }
  cat("\n\tRandom effects:\n")
  cat("\t===============\n")
  print(x$varcorr, digits  = digits)
  if (x$fitted) {
    cat("\n\tLog-likelihood: ", round(x$ll, digits),
        " (df = ", attr(x$ll, "df"), ")", sep ="")
  } else {
    cat("\n\tLog-likelihood: ", NA, sep ="")
  }
  cat("\n\n")
  invisible(x)
}
## TODO: test colored output with knitr
