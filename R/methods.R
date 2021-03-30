##' Set coefficients of a tramME model.
##'
##' Sets the whole vector of coefficients of a tramME model. The parameters of
##' the baseline transformation function should respect the restrictions of
##' the parameter space. This is checked before setting the new parameter values
##' provided that the parameters for the variance components has already been set.
##' If the model contains fixed coefficient parameters, the input should also respect
##' that.
##' When called on a fitted tram object, the function sets it to unfitted and removes
##' all parts that come from the estimation.
##' @param object A \code{tramME} object.
##' @param value Numeric vector of new coefficient values.
##' @return A \code{tramME} object with the new coefficient values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' coef(mod) <- c(-1, 0.5, 1)
##' @importFrom mlt "coef<-"
##' @export
"coef<-.tramME" <- function(object, value) {
  object <- duplicate(object) ## NOTE: force copy
  object$param <- .set_cf(object, value)
  if (!is.null(object$opt)) {
    warning(paste("The model object has already been fitted.",
                  "Removing optimization results."))
    object$opt <- NULL
  }
  return(object)
}


##' Generic method for \code{"varcov<-"}
##' @param object A model object.
##' @param value The new value of the covariance matrix.
##' @param ... Optional inputs.
##' @return An object with the same class as \code{object}, with updated
##'   variance-covariance matrix of random effects.
##' @export
"varcov<-" <- function(object, ..., value)
  UseMethod("varcov<-")


##' Set the values of the random effects covariance matrices of a tramME model.
##'
##' Sets the list containing the covariance matrices of a tramME model. The matrices have
##' to be positive definite. Just as in \code{"coef<-"}, when the function is called
##' on a fitted object, the function will remove the infromation about the optimization.
##'
##' The supplied list has to be named with the same names as implied by the model.
##' Hence, it might be a good idea to call \code{varcov} first, and
##' modify this list to make sure that the input has the right structure.
##'
##' The new values can also be supplied in a form that corresponds to the reparametrization
##' used by the \code{tramTMB} model (see the option \code{as.theta = TRUE}).
##' @param object A \code{tramME} object.
##' @param value A list of positive definite covariance matrices.
##' @param as.theta Logical value, if \code{TRUE}, indicating that the new values are supplied
##'   in their reparameterized form.
##' @param ... Optional arguments (ignored).
##' @return A  \code{tramME} object with the new coefficient values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' vc <- varcov(mod)
##' vc[[1]] <- matrix(c(1, 0, 0, 2), ncol = 2)
##' varcov(mod) <- vc
##' @export
"varcov<-.tramME" <- function(object, as.theta = FALSE, ..., value) {
  object <- duplicate(object) ## NOTE: force copy
  object$param <- .set_vc(object, val = value, as.theta = as.theta)
  if (!is.null(object$opt)) {
    warning(paste("The model object has already been fitted.",
                  "Removing optimization results."))
    object$opt <- NULL
  }
  return(object)
}


##' Generic method for \code{varcov}
##' @param object A model object.
##' @param ... Optional inputs.
##' @return A variance-covariance matrix.
##' @export
varcov <- function(object, ...)
  UseMethod("varcov")


##' Extract the variance-covariance matrix of the random effects
##'
##' Returns the covariance matrix of the random effects as saved in the \code{tramME}
##' object.
##' The returned values correspond to the transformation model parametrization.
##' @param object A \code{tramME} object.
##' @param as.theta Logical value, if \code{TRUE}, the values are returned
##'   in their reparameterized form.
##' @param ... Optional arguments (unused).
##' @return A list of the covariance matrices or a vector of theta values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' varcov(fit)
##' varcov(fit, as.theta = TRUE)
##' @export
varcov.tramME <- function(object, as.theta = FALSE, ...) {
  .get_vc(object, as.theta = as.theta)
}


##' Extract the coefficients of the fixed effects terms.
##' @param object A \code{tramME} object.
##' @param with_baseline If \code{TRUE}, also include the baseline parameters.
##' @param fixed If \code{TRUE}, also include the fixed parameters.
##' @param ... Optional parameters (ignored).
##' @return Numeric vector of parameter values.
##' @examples
##' library("survival")
##' mod <- SurvregME(Surv(time, status) ~ rx + (1 | litter/rx), data = rats,
##'                  dist = "exponential", nofit = TRUE)
##' coef(mod, with_baseline = TRUE)
##' coef(mod, with_baseline = TRUE, fixed = FALSE)
##' @importFrom stats coef
##' @export
coef.tramME <- function(object, with_baseline = FALSE, fixed = TRUE, ...) {
  cf <- .get_cf(object)
  if (with_baseline)
    pargroup <- "fixef"
  else pargroup <- "shift"
  cf[.idx(object, fixed = fixed, pargroup = pargroup)]
}


##' Get the log-likelihood of the model
##' @param object A \code{tramME} object.
##' @param param An optional vector of parameter values in the structure
##'   (beta, theta).
##' @param newdata An optional data.frame to calculate the out-of-sample
##'   log-likelihood.
##' @param ... Optional argument (for consistency with generic).
##' @return A numeric value of the log-likelihood.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' logLik(fit)
##' @importFrom stats logLik update
##' @export
## FIXME: weights & offset arguments. Maybe taking do_update of the object into account.
logLik.tramME <- function(object,
                          param = c(coef(object, with_baseline = TRUE, fixed = FALSE),
                                  varcov(object, as.theta = TRUE)),
                          newdata = NULL, ...) {
  np <- length(coef(object, with_baseline = TRUE, fixed = FALSE)) +
    length(varcov(object, as.theta = TRUE))
  stopifnot(length(param) == np)
  if (!is.null(newdata)) {
    object <- update(object, ctm = object$model$ctm, data = newdata, nofit = TRUE)
  }
  if (any(is.na(param))) {
    ll <- NA
  } else {
    if (!.check_par(object$tmb_obj, param))
      stop("The supplied parameters do not satisfy the parameter constraints.")
    ll <- -object$tmb_obj$fn(param)
  }
  df <- length(param)
  nobs <- sum(object$tmb_obj$env$data$weights)
  structure(ll, df = df, nobs = nobs, class = "logLik")
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
##' @param object A \code{tramME} object.
##' @param object2 A \code{tramME} object.
##' @param ... Optional arguments, for compatibility with the generic. (Ignored)
##' @return A \code{data.frame} with the calculated statistics.
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
  out <- data.frame(npar = sapply(ll_, attr, "df"), logLik = unlist(ll_),
                    AIC = sapply(list(object, object2), "AIC"),
                    BIC = sapply(list(object, object2), "BIC"))
  rownames(out) <- c("Model 1", "Model 2")
  ord <- order(out$npar)
  out <- out[ord, ]
  out$Chisq <- 2 * pmax(0, c(NA, diff(out$logLik)))
  out$`Chisq df` <- c(NA, diff(out$npar))
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
##' @param ... Optional arguments passed to \code{vcov.tramTMB}
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
  pargroup <- match.arg(pargroup)

  b <- .get_cf(object)[.idx(object, fixed = FALSE, pargroup = "fixef")]
  th <- .get_vc(object, as.theta = TRUE)
  pr <- c(b, th)
  if (any(is.na(pr))) {
    out <- matrix(NA, nrow = length(pr), ncol = length(pr))
  } else {
    if ("method" %in% names(list(...))) {
      out <- vcov(object$tmb_obj, par = pr, ...) ## if method is specified try that
    } else {
      if (is.null(object$tmb_obj$env$random) && !object$tmb_obj$env$resid) {
        method <- "analytical" ## when anylitical makes sense do that
      } else method <- "optimHess" ## default
      out <- try(vcov(object$tmb_obj, par = pr, method = method, ...), silent = TRUE)
      if (inherits(out, "try-error")) {
        ## NOTE: numDeriv is often more stable numerically
        out <- vcov(object$tmb_obj, par = pr, method = "numDeriv", ...)
      }
    }
  }
  idx <- .idx(object, pargroup = pargroup, which = parm, pmatch = pmatch)
  out <- as.matrix(out[idx, idx])
  colnames(out) <- names(idx)
  rownames(out) <- names(idx)
  return(out)
}


##' Variances and correlation matrices of random effects
##'
##' This function calculates the variances and correlations from \code{varcov.tramME}.
##' @param x A \code{tramME} object
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
  if (!is.list(vc)) {
    out <- list()
  } else {
    lv <- attr(x$param, "re")$levels
    nl <- sapply(lv, length)
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
  }
  class(out) <- c("VarCorr.tramME", class(out))
  return(out)
}


##' Print method for the variance-correlation parameters of a tramME object
##' @param x A \code{VarCorr.tramME} object.
##' @param sd Logical. Print standard deviations instead of variances.
##' @param digits Number of digits
##' @param ... optional arguments
##' @return Invisibly returns the input VarCorr.tramME object.
##' @export
print.VarCorr.tramME <- function(x, sd = TRUE,
                                digits = max(getOption("digits") - 2L, 3L),
                                ...) {
  if (length(x) == 0) {
    cat("\nNo random effects.\n")
  } else {
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
  }
  cat("\n")
  invisible(x)
}


##' Return variable names.
##'
##' Returns the variable names corresponding the selected group.
##' The returned names are the names as they are used by tramME. For example,
##' when the response is a \code{Surv} object, \code{variable.names} returns
##' the name of that object, and not the names of the variables used to create it.
##' @param object a tramME object (fitted or unfitted)
##' @param which \enumerate{
##'   \item all: all variables,
##'   \item response: response variable,
##'   \item grouping: grouping factors for random effects,
##'   \item shifting: shifting variables,
##'   \item interacting: interacting variables.
##'   }
##' @param ... optional parameters
##' @return A vector of variable names.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' variable.names(mod)
##' variable.names(mod, "response")
##' @importFrom stats variable.names
##' @export
## NOTE: Should REs w/o corresponding FE be addedd to shifting?
variable.names.tramME <- function(object,
    which = c("all", "response", "grouping", "shifting", "interacting"), ...) {
  which <- match.arg(which)
  unname(switch(which,
    grouping = {
      unique(unlist(sapply(object$model$ranef, function(x) all.vars(x[[3]]))))
    },
    all = c(variable.names(object$model$ctm, which = "all", ...),
            variable.names(object, which = "grouping", ...)),
    variable.names(object$model$ctm, which = which, ...)))
}


##' Confidence intervals for tramME model parameters
##'
##' Confidence intervals for model parameters on their original scale.
##' Either Wald CI or profile CI by root finding. Multicore computations
##' are supported in the case of profile confidence intervals, but snow
##' support is yet to be implemented.
##' @param object A \code{tramME} object.
##' @param level Confidence level.
##' @param type Type of the CI: either Wald or profile.
##' @param estimate Logical, add the point estimates in a thrid column.
##' @param parallel Method for parallel computation.
##' @param ncpus Number of cores to use for parallel computation.
##' @param ... Optional parameters.
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

  type <- tolower(match.arg(type))
  pargroup <- match.arg(pargroup)
  plist <- .parallel_default(parallel, ncpus)

  ## --- Indices, point estimates
  b <- .get_cf(object)[.idx(object, fixed = FALSE, pargroup = "fixef")]
  th <- .get_vc(object, as.theta = TRUE)
  pr <- c(b, th)

  idx <- .idx(object, pargroup = pargroup, which = parm, pmatch = pmatch)
  par <- pr[idx]

  if (any(is.na(pr)) || length(par) == 0) {
    nc <- if (estimate) 3 else 2
    ci <- matrix(NA, nrow = length(par), ncol = nc)
    rownames(ci) <- names(par)
    colnames(ci) <- if (nc == 3) c("lwr", "upr", "est") else c("lwr", "upr")
    return(ci)
  }

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
        ## FIXME: add snow support
        stop("No snow support yet")
      }
    } else {
      ci <- lapply(idx, fun)
    }
    ci <- do.call("rbind", ci)
  }

  colnames(ci) <- c("lwr", "upr")
  rownames(ci) <- names(idx)
  if (estimate) {
    ci <- cbind(ci, par)
    colnames(ci) <- c("lwr", "upr", "est")
  }
  return(ci)
}


##' Extract the conditional modes and conditional variances of random effects
##'
##' @param object A \code{tramME} object.
##' @param param An optional vector of parameter values in the structure
##'   (beta, theta).
##' @param newdata An optional \code{data.frame} of new observations for which the
##'   new random effects values are predicted.
##' @param condVar If \code{TRUE}, include the conditional variances as attributes.
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
##' @importFrom stats update
##' @importFrom nlme ranef
##' @aliases ranef
##' @export ranef
##' @export
ranef.tramME <- function(object,
                         param = c(coef(object, with_baseline = TRUE, fixed = FALSE),
                                   varcov(object, as.theta = TRUE)),
                         newdata = NULL,
                         condVar = FALSE, raw = FALSE, ...) {
  np <- length(coef(object, with_baseline = TRUE, fixed = FALSE)) +
    length(varcov(object, as.theta = TRUE))
  stopifnot(length(param) == np)

  if (!is.null(newdata)) {
    object <- update(object, ctm = object$model$ctm, data = newdata, nofit = TRUE)
  }

  if (any(is.na(param))) {
    re <- .get_par(object$tmb_obj)$gamma
    re <- rep(NA, length(re))
  } else {
    re <- .get_par(object$tmb_obj, param)$gamma
  }

  if (raw) {
    return(re)
  }
  if (length(re) == 0) {
    return(NULL)
  }

  re_ <- attr(object$param, "re")
  out <- .re_format(re, re_$termsize, re_$names, re_$blocksize, re_$levels)

  if (condVar && !any(is.na(param))) {
    cv <- TMB::sdreport(object$tmb_obj, par.fixed = param)$diag.cov.random
    cv <- .re_format(cv, re_$termsize, re_$names, re_$blocksize, re_$levels)
    out <- mapply(FUN = function(cm, cv) {
      attr(cm, "condVar") <- cv
      cm
    }, cm = out, cv = cv, SIMPLIFY = FALSE)
  }
  class(out) <- c("ranef.tramME", class(out))
  return(out)
}


##' Residuals of a tramME model
##'
##' Calculates the score residuals of an intercept term fixed at 0.
##' @param object A \code{tramME} object.
##' @param param An optional vector of parameter values in the structure
##'   (beta, theta).
##' @param newdata An optional data.frame.
##' @param ... Optional arguments (currently ignored).
##' @examples
##' library("survival")
##' fit <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats)
##' resid(fit)
##' @importFrom stats residuals update
##' @export
## FIXME: weights & offset arguments. Maybe taking do_update of the object into account.
residuals.tramME <- function(object,
                             param = c(coef(object, with_baseline = TRUE, fixed = FALSE),
                                       varcov(object, as.theta = TRUE)),
                             newdata = NULL, ...) {
  np <- length(coef(object, with_baseline = TRUE, fixed = FALSE)) +
    length(varcov(object, as.theta = TRUE))
  stopifnot(length(param) == np)

  if (!object$tmb_obj$env$resid && !is.null(newdata)) {
    object <- update(object, ctm = object$model$ctm,
                     data = newdata, resid = TRUE, nofit = TRUE)
  } else if (!object$tmb_obj$env$resid) {
    object <- update(object, resid = TRUE, nofit = TRUE)
  } else if (!is.null(newdata)) {
    object <- update(object, ctm = object$model$ctm, data = newdata, nofit = TRUE)
  }

  if (any(is.na(param))) {
    r <- rep(NA, length(object$tmb_obj$env$parameters$alpha0))
  } else {
    if (!.check_par(object$tmb_obj, param))
      stop("The supplied parameters do not satisfy the parameter constraints.")
    r <- object$tmb_obj$resid(param)
  }
  names(r) <- rownames(object$data)
  return(r)
}


##' Print tramME model
##' @param x A \code{tramME} object.
##' @param digits Number of significant digits
##' @param ... Optional arguments (for consistency with the generic)
##' @return The original \code{tramME} object invisibly
##' @export
print.tramME <- function(x, digits = max(getOption("digits") - 2L, 3L), ...) {
  mnm <- .model_name(x)
  cat("\n", mnm, "\n", sep = "")
  cat("\n\tFormula: ")
  print(x$call$formula)
  fitted <-!is.null(x$opt)
  if (fitted) {
    cat("\n\tFitted to dataset ")
    print(x$call$data)
  } else {
    cat("\n\tNot fitted\n")
  }
  fe <- coef(x, fixed = TRUE)
  fe2 <- coef(x, fixed = FALSE)
  fix <- setdiff(names(fe), names(fe2))
  names(fe)[match(fix, names(fe), nomatch = 0L)] <- paste(fix, "(fixed)")
  if (fitted && length(fe) > 0) {
    cat("\n\tFixed effects parameters:\n")
    cat("\t=========================\n\n")
    print(signif(fe, digits))
  } else if (!fitted && all(!is.na(fe)) && length(fe) > 0) {
    cat("\n\tFixed effects parameters set to:\n")
    cat("\t================================\n\n")
    print(signif(fe, digits))
  }
  vc <- VarCorr(x)
  if (length(vc) > 0) {
    if (fitted) {
      cat("\n\tRandom effects parameters:\n")
      cat("\t===========================\n")
      print(vc, digits  = digits)
    } else if (all(!is.na(unlist(vc)))) {
      cat("\n\tRandom effects parameters set to:\n")
      cat("\t=================================\n")
      print(vc, digits  = digits)
    }
  }
  ll <- logLik(x)
  if (!is.na(ll)) {
    cat("\n\tLog-likelihood: ", round(ll, digits),
        " (df = ", attr(ll, "df"), ")", sep ="")
  }
  cat("\n\n")
  invisible(x)
}


##' Summary method for tramME model
##'
##' @param object A \code{tramME} object
##' @param ... Optional arguments (for consistency with the generic)
##' @return A summary.tramME object.
##' @importFrom stats pnorm
##' @export
summary.tramME <- function(object, ...) {
  ll <- logLik(object)
  b <- coef(object, with_baseline = FALSE, fixed = FALSE)
  b2 <- coef(object, with_baseline = FALSE, fixed = TRUE)
  if (length(b) > 0) {
    se <- sqrt(diag(vcov(object, pargroup = "shift")))
  } else se <- numeric(0)
  zval <- b / se
  coef <- cbind(Estimate = b, `Std. Error` = se,
    `z value` = zval,
    `Pr(>|z|)` = 2 * pnorm(abs(zval), lower.tail = FALSE))
  rownames(coef) <- names(b)
  structure(
    list(name    = .model_name(object),
         formula = object$call$formula,
         wtd     = !is.null(model.weights(object$data)),
         fitted  = !is.null(object$opt),
         data    = object$call$data,
         conv    = object$opt$convergence == 0,
         nwarn   = length(object$opt$warnings),
         coef    = coef,
         fixed   = b2[!(names(b2) %in% names(b))],
         varcorr = VarCorr(object),
         ll      = ll),
    class = "summary.tramME")
}


##' Print method for tramME model summary
##'
##' @param x A \code{summary.tramME} object.
##' @param fancy Logical, if \code{TRUE}, use color in outputs.
##' @param ... Optional arguments passed to \code{\link[stats]{printCoefmat}}
##' @return The input summary.tramME object, invisibly.
##' @inheritParams stats::printCoefmat
##' @export
print.summary.tramME <- function(x,
  fancy = !isTRUE(getOption("knitr.in.progress")) && interactive(),
  digits = max(getOption("digits") - 2L, 3L),
  signif.stars = getOption("show.signif.stars"),
  ...) {
  cat("\n", x$name, "\n", sep = "")
  cat("\n\tFormula: ")
  print(x$formula)
  wmsg <- if (x$wtd) " (weighted estimation)" else ""
  if (fancy) {
    if (x$fitted) {
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
  } else {
    if (x$fitted) {
      cat("\n\tFitted to dataset", x$data, wmsg, "\n", sep = " ")
      if (!x$conv)
        cat("\tOptimizer did not achieve convergence!\n") ## TODO: add later optimizer name
      if (x$nwarn > 0)
        cat("\tThere were", x$nwarn, "warning messages captured during optimization.",
            "\n", sep = " ")
    } else {
      cat("\n\tNot fitted\n")
    }
  }
  cat("\n\tFixed effects parameters:\n")
  cat("\t=========================\n\n")
  if (nrow(x$coef) == 0) {
    cat("No estimated shift coefficients.\n")
  } else {
  printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1L:2L,
               tst.ind = 3L, na.print = "NA", ...)
  }
  if (length(x$fixed) > 0) {
    cat("\n\tFixed coefficients:\n")
    cat("\t===================\n\n")
    print(signif(x$fixed, digits))
  }
  cat("\n\tRandom effects:\n")
  cat("\t===============\n")
  print(x$varcorr, digits  = digits)
  cat("\n\tLog-likelihood: ", round(x$ll, digits),
      " (df = ", attr(x$ll, "df"), ")", sep ="")
  cat("\n\n")
  invisible(x)
}

##' Extract model frame from a tramME model
##' @param formula A \code{tramME} object
##' @param ... Optional arguments (currently ignored)
##' @importFrom stats model.frame
##' @export
model.frame.tramME <- function(formula, ...) {
  formula$data
}

##' Generic method for \code{"offset"}
##' @param object An object.
offset <- function(object)
  UseMethod("offset")

##' Default method for \code{"offset"}
##'
##' Overloads the original \code{\link[stats]{offset}} function.
##' @inheritParams stats::offset
offset.default <- function(object)
  stats::offset(object)

##' Get the offset vector of a tramME object.
##' @param object A \code{tramME} object.
offset.tramME <- function(object)
  model.offset(model.frame(object))


##' Get the observation weight vector of a tramME object.
##' @param object A \code{tramME} object.
##' @param ... Optional arguments (ignored).
##' @importFrom stats weights
weights.tramME <- function(object, ...) {
  model.weights(model.frame(object))
}


##' Generic method for \code{"offset<-"}
##' @param object A model object.
##' @param value The new vector of the offsets.
##' @return An object with the same class as \code{object}, with updated
##'   offset vector.
"offset<-" <- function(object, value)
  UseMethod("offset<-")


##' Set the values of the offsets of a tramME model.
##'
##' This method updates the internal \code{tramTMB} object, the \code{model.frame}
##' of the \code{tramME} object and the function call to propagate the change.
##' @note It works only when the \code{tramME} model is defined with \code{do_update = TRUE}.
##' @param object A \code{tramME} object defined with \code{do_update = TRUE}.
##' @param value A vector of new offset values.
##' @return A  \code{tramME} object with the new offset values.
## @examples
## data("sleepstudy", package = "lme4")
## mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE,
##             do_update = TRUE)
## offset(mod)
## offset(mod) <- rep(1, nrow(sleepstudy))
## offset(mod)
## FIXME: this solution allows the tramME model to diverge from the tramTMB object,
## which is very dangerous
## ff <- function(m, a) {offset(m) <- rep(0, nrow(model.frame(m))) + a; print(offset(m)); logLik(m)}
## ff(mm, 1)
## offset(mm)
## mm$tmb_obj$env$data$offset
"offset<-.tramME" <- function(object, value) {
  if (!object$tmb_obj$env$do_update)
    stop("The model is not defined with the option do_update = TRUE. Try updating first.")
  if (!is.null(object$opt)) {
    warning(paste("The model object has already been fitted.",
                  "Removing optimization results."))
    object$opt <- NULL
  }
  stopifnot(nrow(object$data) == length(value))
  value <- as.numeric(value)
  stopifnot(all(!is.na(value)))
  object$data[["(offset)"]] <- value ## 1: update model.frame
  object$tmb_obj$env$data$offset <- value ## 2: update tramTMB
  cl <- match.call()
  oc <- as.list(object$call)
  oc$offset <- cl$value
  object$call <- as.call(oc) ## 3: update call
  return(object)
}


##' Generic method for \code{"weights<-"}
##' @param object A model object.
##' @param value The new vector of the weights.
##' @return An object with the same class as \code{object}, with updated
##'   weight vector.
"weights<-" <- function(object, value)
  UseMethod("weights<-")


##' Set the values of the observation weights of a tramME model.
##'
##' This method updates the internal \code{tramTMB} object, the \code{model.frame}
##' of the \code{tramME} object and the function call to propagate the change.
##' @note It works only when the \code{tramME} model is defined with \code{do_update = TRUE}.
##' @param object A \code{tramME} object defined with \code{do_update = TRUE}.
##' @param value A vector of new weight values.
##' @return A  \code{tramME} object with the new weight values.
## @examples
## library("survival")
## data("eortc", package = "coxme")
## mod <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = eortc, nofit = TRUE,
##                do_update = TRUE)
## weights(mod)
## weights(mod) <- sample(1:3, nrow(eortc), replace = TRUE)
## weights(mod)
"weights<-.tramME" <- function(object, value) {
  if (!object$tmb_obj$env$do_update)
    stop("The model is not defined with the option do_update = TRUE. Try updating first.")
  if (!is.null(object$opt)) {
    warning(paste("The model object has already been fitted.",
                  "Removing optimization results."))
    object$opt <- NULL
  }
  stopifnot(nrow(object$data) == length(value))
  value <- as.numeric(value)
  stopifnot(all(!is.na(value)))
  object$data[["(weights)"]] <- value ## 1: update model.frame
  object$tmb_obj$env$data$weights <- value ## 2: update tramTMB
  cl <- match.call()
  oc <- as.list(object$call)
  oc$weights <- cl$value
  object$call <- as.call(oc) ## 3: update call
  return(object)
}

##' Model matrix for tramME mdoels
##'
##' Creates the model matrix of fixed and random effects corresponding a \code{tramME}
##' model from a \code{data.frame} of response and covariate values.
##' @param object A \code{tramME} object.
##' @param data A \code{data.frame} containing the variable values.
##' @param type Either \code{"fixef"} or \code{"ranef"}.
##' @param with_baseline Logical; indicating whether the returned fixed effects model
##'   matrix should contain the columns corresponding to the baseline transfromation.
##'   (ignored when \code{type = "ranef"})
##' @param ... Additional arguments.
##' @note The model matrix of the random effects is a sparse matrix and it is transposed
##'   to be directly used with Matrix::crossprod which is faster than transposing and
##'   multiplying.
##' @importFrom stats model.matrix
##' @export
model.matrix.tramME <- function(object, data = model.frame(object),
                                type = c("fixef", "ranef"),
                                with_baseline = TRUE,
                                ...) {
  type <- match.arg(type)
  switch(type,
    fixef = {
      if (with_baseline) {
        model.matrix(object$model$ctm, data = data, ...)
      } else {
        if (is.null(object$model$ctm$bases$shifting)) {
          NULL
        } else {
          model.matrix(object$model$ctm$bases$shifting, data = data, ...)
        }
      }
    },
    ranef = {
      if (is.null(object$model$ranef)) {
        NULL
      } else {
        re_terms(object$model$ranef, data = data,
                 negative = object$model$negative)$Zt
      }
    }
  )
}

##' Fit the model.
##' @param object An object.
##' @param ... Optional parameters.
##' @export
fitmod <- function(object, ...) {
  UseMethod("fitmod")
}

##' Call the optimizer on a tramME object
##' @param object A \code{tramME} object.
##' @inheritParams LmME
##' @export
fitmod.tramME <- function(object, initpar = NULL, control = optim_control(), ...) {
  ## NOTE: force copy tramTMB object, to avoid accidentally creating tramMEs that share the
  ## tramTMB
  obj <- duplicate(object$tmb_obj)
  opt <- optim_tramTMB(obj, par = initpar, method = control$method,
                       control = control$control,
                       trace = control$trace, ntry = control$ntry,
                       scale = control$scale)
  parm <- .get_par(obj)
  att <- attributes(object$param)
  param <- .gen_param(parm, fe = att$fe,
                      re = att$re,
                      varnames = att$varnames)
  object$tmb_obj <- obj
  object$param <- param
  object$opt <- opt
  return(object)
}

##' Duplicate a tramME object
##'
##' In general, this is not necessary for the usual usage of tramME.
##' It is only written to avoid errors stemming from the fact that
##' some parts of the tramME object are modified in place.
##' @param object A \code{tramME} object.
##' @param ... Optional arguments (currently ignored).
duplicate.tramME <- function(object, ...) {
  newobj <- object
  par <- list(beta = coef(object, with_baseline = TRUE),
              theta = varcov(object, as.theta = TRUE))
  newobj$tmb_obj <- duplicate(object$tmb_obj)
  att <- attributes(object$param)
  param <- .gen_param(par, fe = att$fe,
                      re = att$re,
                      varnames = att$varnames)
  newobj$param <- param
  return(newobj)
}
