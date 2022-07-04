##' Set coefficients of a tramME model.
##'
##' Sets the whole vector of fixed-effects coefficients of a tramME model.
##' The parameters of the baseline transformation function should respect the
##' restrictions of the parameter space. This is checked before setting the new
##' parameter values provided that the parameters for the variance components has
##' already been set.
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
  pr <- .get_par(object$tmb_obj, fixed = FALSE)
  nb0 <- length(pr$beta0)
  nb <- length(pr$beta)
  stopifnot(length(value) == (nb0 + nb))
  if (all(!is.na(value)) && all(!is.na(object$param$theta))) { ## NOTE: check constraints and update tmb
    if (!.check_par(object$tmb_obj, c(value, object$param$theta))) {
      stop(paste("The assigned parameter values do not satisfy the constraints",
                 "implied by the model.\n\t",
                 "Please check BOTH coef and varcov."))
    }
  }
  pr <- .get_par(object$tmb_obj, c(value, object$param$theta)) ## NOTE: to handle fixed values
  object$param$beta[] <- c(pr$beta0, pr$beta)
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
##'
##' All random effects variance parameters must be supplied. When there are penalized smooth
##' terms in the model variance parameters corresponding to these should also be part of the
##' input list.
##' @param object A \code{tramME} object.
##' @param value A list of positive definite covariance matrices.
##' @param as.theta Logical value, if \code{TRUE}, indicating that the new values are supplied
##'   in their reparameterized form.
##' @param ... Optional arguments (ignored).
##' @return A new \code{tramME} object with the new coefficient values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' vc <- varcov(mod)
##' vc[[1]] <- matrix(c(1, 0, 0, 2), ncol = 2)
##' varcov(mod) <- vc
##' @export
## TODO: It does not support fixed parameter values for theta/varcov atm
"varcov<-.tramME" <- function(object, as.theta = FALSE, ..., value) {
  object <- duplicate(object) ## NOTE: force copy
  ## object$param <- .set_vc(object, val = value, as.theta = as.theta)
  att <- attributes(object$param)
  bls <- c(att$re$blocksize, rep(1, length(att$sm$re_dims)))
  if (!as.theta) {
    stopifnot(identical(lapply(value, dim), lapply(object$param$varcov, dim)))
    stopifnot(all(sapply(value, is.pd)))
    th_ <- .vc2th(value, bls)
  } else {
    stopifnot(length(value) == length(.get_par(object$tmb_obj, fixed = FALSE)$theta))
    th_ <- value
  }
  if (all(!is.na(object$param$beta)) && all(!is.null(th_))) { ## NOTE: check constraints and update tmb
    if (!is.null(object$tmb_obj$env$map$beta))
      b_ <- object$param$beta[!is.na(object$tmb_obj$env$map$beta)]
    else b_ <- object$param$beta
    if (!.check_par(object$tmb_obj, c(b_, th_))) {
      stop(paste("The assigned parameter values do not satisfy the constraints",
                 "implied by the model.\n\t",
                 "Please check BOTH coef and varcov."))
    }
  }
  vc_ <- object$param$varcov
  if (as.theta) {
    val_ <- .th2vc(value, bls)
  } else {
    val_ <- value
  }
  ## FIXME: .upd_param could be used here
  object$param$theta[] <- th_
  att <- attributes(object$param$varcov)
  object$param$varcov <- mapply(function(old, new) {old[] <- new[]; old}, ## NOTE: to keep the names
                                vc_, val_, SIMPLIFY = FALSE)
  attributes(object$param$varcov) <- att
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
##' @param full Logical value; if \code{TRUE}, return all random effects elements,
##'   if \code{FALSE}, do not return the random effects parameters of the smooth
##'   terms.
##' @param ... Optional arguments (unused).
##' @return A list of the covariance matrices or a vector of theta values.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' varcov(fit)
##' varcov(fit, as.theta = TRUE)
##' @export
## FIXME: fixed parameters?
varcov.tramME <- function(object, as.theta = FALSE, full = FALSE, ...) {
  ## ty <- attr(object$param$theta, "type")
  if (!as.theta) {
    bls <- c(attr(object$param, "re")$blocksize,
             rep(1, length(attr(object$param, "sm")$re_dims)))
    ## ty <- ty[cumsum(bls * (bls + 1) / 2)]
    ty <- attr(object$param$varcov, "type")
    out <- object$param$varcov
  } else {
    ty <- attr(object$param$theta, "type")
    out <- object$param$theta
  }
  if (full) {
    g <- rep(TRUE, length(ty))
  } else {
    g <- ty == "re"
  }
  return(out[g])
}


##' Extract the coefficients of the fixed effects terms.
##' @param object A \code{tramME} object.
##' @param with_baseline If \code{TRUE}, also include the baseline parameters and the
##'   fixed effects parameters from the smooth terms.
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
## FIXME: should probably rename with_baseline to full but keep with_baseline for backward compatibility
coef.tramME <- function(object, with_baseline = FALSE, fixed = TRUE, ...) {
  cf <- object$param$beta
  g <- rep(TRUE, length(cf))
  if (!with_baseline)
    g <- g & (attr(cf, "type") == "sh")
  if (!fixed)
    g <- g & !attr(cf, "fixed")
  cf[g]
}


##' Get the log-likelihood of the tramME model
##'
##' Evaluates the log-likelihood function. New parameter values and data can
##' optionally be supplied. In the latter case, the function returns the
##' out-of-sample log-likelihood.
##'
##' @details
##'
##' By default, \code{param} is set to the estimated (or previously set)
##'   parameters. If the parameter vectors in the model are incomplete (contain
##'   \code{NA} elemets), the returned log-likelihood will also be \code{NA},
##'   unless the user provides new values.
##'
##' Setting \code{type = "fix_smooth"} fixes the random effects terms that
##'   correspond to penalized smooths at their estimated values, so that they
##'   are not refitted when \code{newdata} is supplied. This is consistent with
##'   treating these parameter regularized fixed terms, i.e. as 'new-style'
##'   random effects described by Hodges (2014, Chapter 13).
##'
##' The \code{"fix_smooth"} and \code{"penalized"} options for \code{type} are
##'   just for convenience.  The same functionality can be achieved by setting
##'   \code{param$gamma} to the desired values.  \code{"penalized"} respects the
##'   values of \code{param$gamma} if both are supplied, while
##'   \code{"fix_smooth"} overwrites them with the fitted values if there are
##'   ambiguities.
##'
##' @section Type of the log-likelihood:
##'
##' By default, \code{logLik} calculates the _integrated_ (or marginal)
##'   log-likelihood by integrating over the random effects. By fixing the
##'   random effects, the value of the log-likelihood changes, because TMB won't
##'   integrate over these random effects.  This will result in the _penalized_
##'   log-likelihood (conditional log-likelihood + penalty for smooth terms and
##'   random effects, see example).
##'
##' By setting \code{type = "penalized"}, the function will 'fix' all random
##'   effects and penalized parameters of the smooth terms at their predicted
##'   levels, and calcualte the penalized log-likelihood. In this sense, setting
##'   \code{type = "fix_smooth"} will result in a hybrid log-likelihood value,
##'   where the 'true' random effects (c.f.  Hodges 2014, Ch. 13) are integrated
##'   out, while it includes the penalty values for the penalized parameters of
##'   the smooths terms.
##'
##' In general, it is not clear which type of log-likelihood we should calculate
##'   when we want to evaluate models based on their out-of-sample
##'   log-likelihood values.  The context and the model setup are key in these
##'   cases. Please make sure you know what you want to calculate to avoid
##'   misunderstandings.
##'
##' @references
##'
##' Hodges, James S. (2014). Richly Parameterized Linear Models: Additive, Time
##'   Series, and Spatial Models Using Random Effects. Chapman & Hall/CRC Texts
##'   in Statistical Science Series.
##'
##' @param object A \code{tramME} object.
##' @param param An optional named list of parameter values (beta and theta).
##'   See details.  Optionally, gamma elements can also be added, which leads to
##'   'fixing' those random effects terms at the supplied values.
##' @param newdata An optional data.frame to calculate the out-of-sample
##'   log-likelihood.
##' @param type The type of the likelihood to be calculated:
##'   \itemize{
##'     \item integrated (default when \code{newdata = NULL}): The marginal
##'           log-likelihood, calculated by integrating out the random effects.
##'     \item fix_smooth (default when \code{newdata} is supplied): Treating the
##'           penalized parameters of the smooth terms as fixed at their
##'           posterior mode predictions and only integrating out the 'true'
##'           random effects. (Consistent with the functionality of
##'           \code{\link[tramME]{ranef.tramME}} and
##'           \code{\link[tramME]{residuals.tramME}} when
##'           \code{fix_smooth = TRUE}.)
##'     \item penalized: Treat all parameters as fixed, return the penalized
##'           log-likelihood (conditional log-likelihood + penalty for smooth
##'           terms and random effects). This is equivalent to fixing all random
##'           effect values.
##'   }
##'   See details.
##' @param ... Optional argument (for consistency with generic).
##' @return A numeric value of the log-likelihood.
##'
##' @examples
##'
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' logLik(fit)
##'
##' data("mcycle", package = "MASS")
##' fit <- LmME(accel ~ s(times), data = mcycle)
##' logLik(fit) < logLik(fit, type = "penalized")
##'
##' @importFrom stats logLik update
##' @export
## TODO: weights & offset arguments. Maybe taking do_update of the object into account.
## TODO: option to calculate the 'unpenalized' logLik
## TODO: what should we use for OOS loglik (w/ REs + smooths) integrated? penalized?
## unpenalized? partially integrated/penalized?
logLik.tramME <- function(object, param = NULL, newdata = NULL,
                          type = c("integrated", "fix_smooth", "penalized"),
                          ...) {
  type <- match.arg(type, several.ok = TRUE)
  if (!is.null(newdata) && length(type) > 1) type <- "fix_smooth"
  else type <- type[1L]
  rf <- switch(type, fix_smooth = "smooth", penalized = "all", 0)
  param <- .default_param(object, param, rf)

  if (!is.null(newdata) || !is.null(param$gamma)) {
    if (!is.null(param$gamma))
      object$call$fixed[names(param$gamma)] <- param$gamma
    ex <- quote(update(object, ctm = object$model$ctm,
                       smooth = .tramME_smooth(object),
                       nofit = TRUE))
    if (!is.null(newdata))
      ex$data <- quote(newdata)
    object <- eval(ex)
  }
  param <- c(param$beta, param$theta)
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
## TODO: more than two models
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
##' @param pargroup The name of the parameter group to return:
##'   \itemize{
##'     \item all: All parameters.
##'     \item fixef: Fixed effects parameters.
##'     \item shift: Shift parameters.
##'     \item baseline: Parameters of the baseline transformation function.
##'     \item ranef: Variance components parameters.
##'     \item smooth: Paramaters that belong to the smooth shift terms
##'       (both FE and smoothing parameters).
##'   }
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
## FIXME: vcov of penalized parameters of smooth terms (random effects)
vcov.tramME <- function(object, parm = NULL,
                        pargroup = c("all", "fixef", "shift", "baseline", "ranef", "smooth"),
                        pmatch = FALSE, ...) {
  pargroup <- match.arg(pargroup)
  b <- coef(object, with_baseline = TRUE, fixed = FALSE)
  th <- varcov(object, as.theta = TRUE, full = TRUE)
  pr <- c(b, th)
  if (any(is.na(pr))) {
    out <- matrix(NA, nrow = length(pr), ncol = length(pr))
  } else {
    if ("method" %in% names(list(...))) {
      out <- vcov(object$tmb_obj, par = pr, ...) ## if method is specified try that
    } else {
      if (is.null(object$tmb_obj$env$random) && !object$tmb_obj$env$resid) {
        method <- "analytical" ## when analytical makes sense do that
      } else method <- "optimHess" ## default
      out <- try(vcov(object$tmb_obj, par = pr, method = method, ...), silent = TRUE)
      if (inherits(out, "try-error")) {
        ## NOTE: numDeriv is often more stable numerically
        out <- vcov(object$tmb_obj, par = pr, method = "numDeriv", ...)
      }
    }
  }
  g <- .choose_parm(object, parm = parm, pargroup = pargroup, pmatch = pmatch,
                    fixed = FALSE)
  out <- out[g, g, drop = FALSE]
  colnames(out) <- rownames(out) <- names(g[g])
  return(out)
}

## FIXME: move to util
.choose_parm <- function(object, parm, pargroup, pmatch, fixed = FALSE) {
  b  <- coef(object, with_baseline = TRUE, fixed = TRUE)
  th <- varcov(object, as.theta = TRUE, full = TRUE)
  nm <- c(names(b), names(th))
  fx <- c(attr(object$param$beta, "fixed"), attr(object$param$theta, "fixed"))
  ty <- c(attr(object$param$beta, "type"), attr(object$param$theta, "type"))
  if (!fixed) {
    ty <- ty[!fx]
    nm <- nm[!fx]
  }
  g  <- switch(pargroup, all = rep(TRUE, length(ty)), fixef = ty %in% c("bl", "sh"),
               shift = ty == "sh", baseline = ty == "bl", ranef = ty == "re",
               smooth = ty == "sm")
  if (is.character(parm)) {
    if (pmatch) {
      g2 <- lapply(parm, grepl, nm)
      if (length(g2) > 1) {
        g2 <- do.call("|", g2)
      } else {
        g2 <- g2[[1]]
      }
      g <- g & g2
    } else {
      g <- g & (nm %in% parm)
    }
  }
  if (is.numeric(parm)) {
    g2 <- seq_along(g)
    g2[g] <- seq_along(g2[g])
    g2[!g] <- 0
    g <- g2 %in% parm
  }
  names(g) <- nm
  g
}

##' Variances and correlation matrices of random effects
##'
##' This function calculates the variances and correlations from \code{varcov.tramME}.
##'
##' The function only returns the correlation matrices that belong to actual random effects
##' (defined for groups in the data) and ignores the random effects parameters of the smooth
##' shift terms. To extract these, the user should use \code{varcov} with \code{full = TRUE}.
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
  vc <- varcov(x, full = FALSE)
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
##' Returns the variable names corresponding to different variable groups in a tramME
##' model.
##'
##' @details
##' The returned names are the names as they are used by tramME. For example,
##' when the response is a \code{Surv} object, \code{variable.names} returns
##' the name of that object, and not the names of the variables used to create it.
##'
##' @param object a tramME object (fitted or unfitted)
##' @param which \enumerate{
##'   \item all: all variables,
##'   \item response: response variable,
##'   \item grouping: grouping factors for random effects,
##'   \item shifting: shifting variables,
##'   \item interacting: interacting variables,
##'   \item smooth: variables in smooth terms,
##'   \item ranef: all random effects variables (covariates with random slopes and grouping
##'   factors).
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
## FIXME: Should REs w/o corresponding FE be addedd to shifting?
variable.names.tramME <- function(object,
  which = c("all", "response", "grouping", "shifting", "interacting", "smooth",
            "ranef"),
  ...) {
  uuapply <- function(x, fun)  unique(unlist(lapply(x, fun)))
  chv2av <- function(x) uuapply(x, function(y) all.vars(str2lang(y)))
  which <- match.arg(which)
  unname(switch(which,
    grouping = {
      uuapply(object$model$ranef, function(x) all.vars(x[[3]]))
    },
    ranef = {
      uuapply(object$model$ranef, function(x) all.vars(x))
    },
    smooth = {
      uuapply(object$model$smooth, function(x) {
        v <- x$term
        if (x$by != "NA") v <- c(v, x$by)
        chv2av(v)
      })
    },
    all = unique(c(variable.names(object$model$ctm, which = "all", ...),
                   variable.names(object, which = "smooth", ...),
                   variable.names(object, which = "ranef", ...))),
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
## TODO: add bootstrap, calculate CIs for random effects (e.g. parameters of smooth terms)
confint.tramME <- function(object, parm = NULL,
                           level = 0.95,
                           pargroup = c("all", "fixef", "shift", "baseline", "ranef", "smooth"),
                           type = c("Wald", "wald", "profile"),
                           estimate = FALSE,
                           pmatch = FALSE,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("profile.ncpus", 1L), ...) {

  type <- tolower(match.arg(type))
  pargroup <- match.arg(pargroup)
  plist <- .parallel_default(parallel, ncpus)

  ## --- Point estimates
  b <- coef(object, with_baseline = TRUE, fixed = FALSE)
  th <- varcov(object, as.theta = TRUE, full = TRUE)
  pr <- c(b, th)

  ## --- Select parameters
  g <- .choose_parm(object, parm = parm, pargroup = pargroup, pmatch = pmatch,
                    fixed = FALSE)
  idx <- which(g)
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
        ## TODO: add snow support
        stop("No snow support yet")
      }
    } else {
      ci <- lapply(idx, fun)
    }
    ci <- do.call("rbind", ci)
  }

  colnames(ci) <- c("lwr", "upr")
  rownames(ci) <- names(g[g])
  if (estimate) {
    ci <- cbind(ci, par)
    colnames(ci) <- c("lwr", "upr", "est")
  }
  return(ci)
}


##' Point estimates and conditional variances of random effects.
##'
##' Extract the conditional modes and conditional variances of random effects in
##' a formatted or unformatted way.
##'
##' @details
##'
##' \code{raw = TRUE} returns the whole vector of random effects (i.e. with
##'   parameters of smooth shift terms), while \code{raw = FALSE} only returns
##'   the formatted list of actual random effects (i.e. for grouped
##'   observations) values. For the conceptual differences between the two types
##'   of random effects, see Hodges (2014, Chapter 13).
##'
##' The conditional variances of the fixed random effects are set to \code{NA}.
##'
##' @section Warning:
##'
##' The function has several optional arguments that allow great flexibilty
##'   beyond its most basic usage. The user should be careful with setting
##'   these, because some combinations might not return sensical results.  Only
##'   limited sanity checks are performed.
##'
##' @inherit logLik.tramME references
##' @inheritParams logLik.tramME
##' @param newdata An optional \code{data.frame} of new observations for which the
##'   new random effects values are predicted.
##' @param fix_smooth Logical; it is set to \code{TRUE} by default, if
##'   \code{newdata} is supplied.  The random effects parameters corresponding
##'   the smooth terms are fixed and not fitted (posterior mode) to
##'   \code{newdata} instead they are treated just like fixed effects
##'   parameters. See details.
##' @param condVar If \code{TRUE}, include the conditional variances as attributes.
##'   Only works with \code{raw = FALSE}.
##' @param raw Return the unformatted RE estimates as fitted by the model.
##' @param ... Optional arguments (for consistency with generic)
##' @return Depending on the value of \code{raw}, either a numeric vector or a
##'   \code{ranef.tramME} object which contains the conditional mode and variance
##'   estimates by grouping factors.
##' @examples
##'
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy, order = 5)
##' ranef(fit, raw = TRUE)
##' ranef(fit)
##'
##' @importFrom stats update
##' @importFrom nlme ranef
##' @aliases ranef
##' @export ranef
##' @export
ranef.tramME <- function(object, param = NULL, newdata = NULL,
                         fix_smooth = !is.null(newdata), condVar = FALSE,
                         raw = FALSE, ...) {
  if (is.null(object$model$smooth) && is.null(object$model$ranef))
      return(NULL)
  param <- .default_param(object, param, if (fix_smooth) "smooth" else 0)
  ## NOTE: if fix_smooth and no other REs, we know the size of the
  ## output and can return it w/o any other steps
  if (fix_smooth && !is.null(object$model$smooth) &&
      is.null(object$model$ranef)) {
    re <- as.vector(object$param$gamma)
    if (raw) return(re)
    return(structure(list(), class = c("ranef.tramME", "list")))
  }
  if (!is.null(newdata) || !is.null(param$gamma)) {
    if (!is.null(param$gamma))
      object$call$fixed[names(param$gamma)] <- param$gamma
    ex <- quote(update(object, ctm = object$model$ctm,
                       smooth = .tramME_smooth(object),
                       nofit = TRUE))
    if (!is.null(newdata))
      ex$data <- quote(newdata)
    object <- eval(ex)
  }
  pg <- param$gamma
  param <- c(param$beta, param$theta)
  if (any(is.na(param))) {
    re <- object$param$gamma
    re <- rep(NA, length(re))
    re[names(pg)] <- pg
  } else {
    re <- .get_par(object$tmb_obj, param)$gamma
  }
  attributes(re) <- NULL
  if (raw) {
    return(re)
  }
  if (length(re) == 0) {
    return(NULL)
  }
  ## NOTE: remove REs that belong to penalized smooth terms
  nsm <- attr(object$param$gamma, "type") != "sm"
  re <- re[nsm]
  re_ <- attr(object$param, "re")
  out <- .re_format(re, re_$termsize, re_$names, re_$blocksize, re_$levels)
  if (condVar && !any(is.na(param))) {
    cv <- TMB::sdreport(object$tmb_obj, par.fixed = param)$diag.cov.random
    ## NOTE: remove SEs that belong to penalized smooth terms
    fg <- attr(object$param$gamma, "fixed")
    cv <- cv[nsm[!fg]] ## drop
    fg <- fg[nsm]
    cvf <- rep(NA, length(fg))
    cvf[!fg] <- cv
    cvf <- .re_format(cvf, re_$termsize, re_$names, re_$blocksize, re_$levels)
    out <- mapply(FUN = function(cm, cv) {
      attr(cm, "condVar") <- cv
      cm
    }, cm = out, cv = cvf, SIMPLIFY = FALSE)
  }
  class(out) <- c("ranef.tramME", class(out))
  return(out)
}


##' Residuals of a tramME model
##'
##' Calculates the score residuals of an intercept term fixed at 0.
##'
##' @param newdata An optional \code{data.frame} of observations for which we
##'   want to calculate the residuals.
##' @inheritParams ranef.tramME
##' @examples
##'
##' library("survival")
##' fit <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats)
##' resid(fit)
##'
##' @importFrom stats residuals update
##' @export
## TODO: more detailed descriptions
## TODO: alternative: quantile residuals?
## TODO: weights & offset arguments. Maybe taking do_update of the object into account.
residuals.tramME <- function(object,
                             param = NULL,
                             newdata = NULL,
                             fix_smooth = !is.null(newdata),
                             ...) {
  param <- .default_param(object, param, if (fix_smooth) "smooth" else 0)
  if (!object$tmb_obj$env$resid || !is.null(newdata) ||
      !is.null(param$gamma)) {
    if (!is.null(param$gamma))
      object$call$fixed[names(param$gamma)] <- param$gamma
    ex <- quote(update(object, ctm = object$model$ctm,
                       smooth = .tramME_smooth(object),
                       nofit = TRUE, resid = TRUE))
    if (!is.null(newdata)) ex$data <- quote(newdata)
    object <- eval(ex)
  }
  param <- c(param$beta, param$theta)
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

## @param re_fix What parts of the random effects vector should be fixed?
## smooth: penalized smoozhing parameters; all: all random effecst
## It handles the overlaps between \code{param$gamma} and \code{re_fix}
## differently. \code{re_fix = "smooth"} overwrites while \code{re_fix = "all"}
## respects the manually supplied values.
.default_param <- function(object, param, re_fix) {
  if (is.null(param)) param <- list()
  stopifnot(is.list(param))
  ## -- defaults for beta and theta
  if (is.null(param$beta)) param$beta <- coef(object, with_baseline = TRUE, fixed = FALSE)
  if (is.null(param$theta)) param$theta <- varcov(object, as.theta = TRUE, full = TRUE)
  stopifnot(length(param$beta) == length(coef(object, with_baseline = TRUE, fixed = FALSE)))
  stopifnot(length(param$theta) == length(varcov(object, as.theta = TRUE, full = TRUE)))
  ## -- fixing random effects
  g <- object$param$gamma
  if (length(g)) stopifnot(all(names(param$gamma) %in% names(g)))
  if (re_fix == "smooth" && !is.null(object$model$smooth)) {
    g <- g[attr(g, "type") == "sm"]
    param$gamma[names(g)] <- g
  } else if (re_fix == "all" && length(g)) {
    nms <- setdiff(names(g), names(param$gamma))
    param$gamma[nms] <- g[nms]
  }
  return(param)
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
  formula <- x$model$formula
  if (inherits(formula, "fake_formula"))
    formula <- eval(x$call$formula, envir = environment(formula))
  print(formula)
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
  if (!is.null(x$model$smooth)) {
    sm <- edf_smooth(x)
    cat("\n\tSmooth shift terms (edf):\n")
    cat("\t=========================\n\n")
    print(signif(sm, digits))
  }
  vc <- VarCorr(x)
  if (length(vc) > 0) {
    if (fitted) {
      cat("\n\tRandom effects parameters:\n")
      cat("\t==========================\n")
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
        " (npar = ", attr(ll, "df"), ")", sep ="")
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
  smooth <- data.frame(edf = edf_smooth(object))
  formula <- object$model$formula
  if (inherits(formula, "fake_formula"))
    formula <- eval(object$call$formula, envir = environment(formula))
  structure(
    list(name    = .model_name(object),
         formula = formula,
         wtd     = !is.null(model.weights(object$data)),
         fitted  = !is.null(object$opt),
         data    = object$call$data,
         conv    = object$opt$convergence == 0,
         nwarn   = length(object$opt$warnings),
         coef    = coef,
         fixed   = b2[!(names(b2) %in% names(b))],
         varcorr = VarCorr(object),
         smooth  = smooth,
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
  if (length(x$fixed)) {
    cat("\n\tFixed coefficients:\n")
    cat("\t===================\n\n")
    print(signif(x$fixed, digits))
  }
  if (length(x$smooth)) {
    cat("\n\tSmooth shift terms:\n")
    cat("\t===================\n\n")
    printCoefmat(x$smooth, digits = digits,
                 signif.stars = signif.stars, ...)
  }
  if (length(x$varcorr)) {
    cat("\n\tRandom effects:\n")
    cat("\t===============\n")
    print(x$varcorr, digits  = digits)
  }
  cat("\n\tLog-likelihood: ", round(x$ll, digits),
      " (npar = ", attr(x$ll, "df"), ")", sep ="")
  cat("\n\n")
  invisible(x)
}


## Generic method for \code{"offset"}
## @param object An object.
offset <- function(object)
  UseMethod("offset")

## Default method for \code{"offset"}
##
## Overloads the original \code{\link[stats]{offset}} function.
## @inheritParams stats::offset
offset.default <- function(object)
  stats::offset(object)

## Get the offset vector of a tramME object.
## @param object A \code{tramME} object.
offset.tramME <- function(object)
  model.offset(model.frame(object))


## Get the observation weight vector of a tramME object.
## @param object A \code{tramME} object.
## @param ... Optional arguments (ignored).
##' @importFrom stats weights
weights.tramME <- function(object, ...) {
  model.weights(model.frame(object))
}


## Generic method for \code{"offset<-"}
## @param object A model object.
## @param value The new vector of the offsets.
## @return An object with the same class as \code{object}, with updated
##   offset vector.
"offset<-" <- function(object, value)
  UseMethod("offset<-")


## Set the values of the offsets of a tramME model.
##
## This method updates the internal \code{tramTMB} object, the \code{model.frame}
## of the \code{tramME} object and the function call to propagate the change.
## @note It works only when the \code{tramME} model is defined with \code{do_update = TRUE}.
## @param object A \code{tramME} object defined with \code{do_update = TRUE}.
## @param value A vector of new offset values.
## @return A  \code{tramME} object with the new offset values.
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


## Generic method for \code{"weights<-"}
## @param object A model object.
## @param value The new vector of the weights.
## @return An object with the same class as \code{object}, with updated
##   weight vector.
"weights<-" <- function(object, value)
  UseMethod("weights<-")


## Set the values of the observation weights of a tramME model.
##
## This method updates the internal \code{tramTMB} object, the \code{model.frame}
## of the \code{tramME} object and the function call to propagate the change.
## @note It works only when the \code{tramME} model is defined with \code{do_update = TRUE}.
## @param object A \code{tramME} object defined with \code{do_update = TRUE}.
## @param value A vector of new weight values.
## @return A  \code{tramME} object with the new weight values.
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

## FIXME: remove
## Fit the model.
## @param object An object.
## @param ... Optional parameters.
## @export
fitmod <- function(object, ...) {
  UseMethod("fitmod")
}

## Call the optimizer on a tramME object
## @param object A \code{tramME} object.
## @inheritParams LmME
## @export
fitmod.tramME <- function(object, initpar = NULL, control = optim_control(), ...) {
  ## NOTE: force copy tramTMB object, to avoid accidentally creating tramMEs that share the
  ## tramTMB
  obj <- duplicate(object$tmb_obj)
  opt <- optim_tramTMB(obj, par = initpar, method = control$method,
                       control = control$control,
                       trace = control$trace, ntry = control$ntry,
                       scale = control$scale)
  ## parm <- .get_par(obj)
  ## att <- attributes(object$param)
  ## param <- .gen_param(parm, fe = att$fe,
  ##                     re = att$re,
  ##                     varnames = att$varnames)
  param <- .upd_param(object$param, obj)
  object$tmb_obj <- obj
  object$param <- param
  object$opt <- opt
  return(object)
}

## Duplicate a tramME object
##
## In general, this is not necessary for the usual usage of tramME.
## It is only written to avoid errors stemming from the fact that
## some parts of the tramME object are modified in place.
## @param object A \code{tramME} object.
## @param ... Optional arguments (currently ignored).
## @return A copy of the original tramME object
duplicate.tramME <- function(object, ...) {
  newobj <- object
  newobj$tmb_obj <- duplicate(object$tmb_obj)
  return(newobj)
}
