##' Mixed-effects Additive Parametric Cox Regression Model
##'
##' Estimates a mixed-effects additive transformation model with flexible smooth
##' parameterization for the baseline transformation (log-cumulative baseline
##' hazard) and the inverse link set to the CDF of the standard minimum extreme
##' value distribution (see Hothorn et al., 2018).
##'
##' @inheritParams tramME
##' @details
##'
##' The model extends \code{\link[tram:Coxph]{tram::Coxph}} with random effects and
##'   (optionally penalized) additive terms. For details on mixed-effect
##'   transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with positive
##'   parameters (i.e. \code{negative = FALSE} in \code{\link[tram]{tram}}).
##'
##' @inherit tramME references
##' @return A \code{CoxphME} model object.
##' @examples
##' library("survival")
##' rats$litter <- factor(rats$litter)
##' m <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats,
##'              log_first = TRUE)
##' summary(m)
##' @importFrom tram Coxph
##' @export
CoxphME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                    silent = TRUE, resid = FALSE, do_update = FALSE,
                    estinit = TRUE, initpar = NULL,
                    fixed = NULL, nofit = FALSE,
                    control = optim_control(),
                    ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Coxph"
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("CoxphME", class(out))
  return(out)
}


##' Mixed-effects Additive Continuous Outcome Logistic Regression Model
##'
##' Estimates a mixed-effects additive transformation model with flexible smooth
##' parameterization for the baseline transformation and the inverse link set to
##' the CDF of the standard logistic distribution (see Hothorn et al., 2018).
##'
##' @inheritParams tramME
##' @details
##'
##' The model extends \code{\link[tram:Colr]{tram::Colr}} with random effects and
##'   (optionally penalized) additive terms. For details on mixed-effect
##'   transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with positive
##'   parameters (i.e. \code{negative = FALSE} in \code{\link[tram]{tram}}).
##'
##' @inherit tramME references
##' @return A \code{ColrME} model object.
##' @examples
##' data("neck_pain", package = "ordinalCont")
##' m <- ColrME(vas ~ time * laser + (1 | id), data = neck_pain,
##'             bounds = c(0, 1), support = c(0, 1), order = 6)
##' summary(m)
##' @importFrom tram Colr
##' @export
ColrME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                   silent = TRUE, resid = FALSE, do_update = FALSE,
                   estinit = TRUE, initpar = NULL,
                   fixed = NULL, nofit = FALSE,
                   control = optim_control(),
                   ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Colr"
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("ColrME", class(out))
  return(out)
}


##' Non-normal (Box-Cox-type) Linear Mixed-effects Additive Regression Model
##'
##' Estimates a mixed-effects additive transformation model with flexible smooth
##' parameterization for the baseline transformation and the inverse link set to
##' the CDF of the standard Gaussian distribution (see Hothorn et al., 2018).
##'
##' @inheritParams tramME
##' @details
##'
##' The model extends \code{\link[tram:BoxCox]{tram::BoxCox}} with random effects and
##'   (optionally penalized) additive terms. For details on mixed-effect
##'   transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with negative
##'   parameters (i.e. \code{negative = TRUE} in \code{\link[tram]{tram}}).
##'
##' @inherit tramME references
##' @return A \code{BoxCoxME} model object.
##' @examples
##' data("sleepstudy", package = "lme4")
##' m <- BoxCoxME(Reaction ~ s(Days) + (Days | Subject), data = sleepstudy)
##' summary(m)
##' @importFrom tram BoxCox
##' @export
BoxCoxME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                     silent = TRUE, resid = FALSE, do_update = FALSE,
                     estinit = TRUE, initpar = NULL,
                     fixed = NULL, nofit = FALSE,
                     control = optim_control(),
                     ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "BoxCox"
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("BoxCoxME", class(out))
  return(out)
}


##' Mixed-effects Additive Lehmann-alternative Linear Regression Model
##'
##' Estimates a mixed-effects additive transformation model with flexible smooth
##' parameterization for the baseline transformation and the inverse link set to
##' the CDF of the standard maximum extreme value distribution (see Hothorn et
##' al., 2018).
##'
##' @inheritParams tramME
##' @details
##'
##' The model extends \code{\link[tram:Lehmann]{tram::Lehmann}} with random
##'   effects and (optionally penalized) additive terms. For details on
##'   mixed-effect transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with negative
##'   parameters (i.e. \code{negative = TRUE} in \code{\link[tram]{tram}}).
##'
##' @inherit tramME references
##' @return A \code{LehmannME} model object.
##' @examples
##' data("sleepstudy", package = "lme4")
##' m <- LehmannME(Reaction ~ s(Days) + (Days | Subject), data = sleepstudy)
##' summary(m)
##' @importFrom tram Lehmann
##' @export
LehmannME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      silent = TRUE, resid = FALSE, do_update = FALSE,
                      estinit = TRUE, initpar = NULL,
                      fixed = NULL, nofit = FALSE,
                      control = optim_control(),
                      ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Lehmann"
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("LehmannME", class(out))
  return(out)
}


##' Mixed-effects Additive Transformation Models for Ordered Categorical
##' Responses
##'
##' Estimates mixed-effects additive transformation models for ordered
##' categorical responses with various link functions.
##'
##' @inheritParams tramME
##' @inheritParams tram::Polr
##' @details
##'
##' The transformation function is parameterized as a step function on a scale
##'   defined by the link function (see Hothorn et al., 2018).
##'
##' The model extends \code{\link[tram:Polr]{tram::Polr}} with random effects and
##'   (optionally penalized) additive terms. For details on mixed-effect
##'   transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with negative
##'   parameters (i.e. \code{negative = TRUE} in \code{\link[tram]{tram}}).
##'
##' @inherit tramME references
##' @return A \code{PolrME} model object.
##' @examples
##' data("soup", package = "ordinal")
##' m <- PolrME(SURENESS | SOUPFREQ ~ PROD + (1 | RESP/PROD),
##'             data = soup, method = "probit")
##' summary(m)
##' @importFrom tram Polr
##' @export
PolrME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                   method = c("logistic", "probit", "loglog", "cloglog"),
                   silent = TRUE, resid = FALSE, do_update = FALSE,
                   estinit = TRUE, initpar = NULL,
                   fixed = NULL, nofit = FALSE,
                   control = optim_control(),
                   ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Polr"
  args$method <- match.arg(method)
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("PolrME", class(out))
  return(out)
}

## A helper function to force the evaluation of name types It helps to avoid
## mixups stemming from the unfortunate but necessary mixing of standard and
## non-standard evaluation chains when setting up tramME models.
get_names <- function(args, env) {
  ## nms <- names(which(sapply(args, is.name)))
  ## args[nms] <- mget(nms, env, inherits = TRUE)
  if (length(idx <- which(sapply(args, is.name)))) {
    nms <- args[idx]
    nms <- sapply(nms, deparse1)
    args[names(nms)] <- mget(nms, env, inherits = TRUE)
  }
  return(args)
}

##' Mixed-effects Additive transformation models
##'
##' A general function to define and fit \code{tramME} models.
##'
##' @details
##'
##' The specific model functions (\code{\link[tramME]{LmME}},
##' \code{\link[tramME]{BoxCoxME}}, \code{\link[tramME]{ColrME}}, etc.) are
##' wrappers around this function.
##'
##' For a general description of the transformation model family, see Hothorn et
##'   al. (2018), for details on the mixed-effects extension, see Tamasi and
##'   Hothorn (2021).
##'
##' @section Warning:
##'
##' Typically, the \code{tramME} function shouldn't be called directly; it is
##'   only exported to allow the advanced users to define their \code{tramME}
##'   models in a more flexible way from their basic building blocks.
##'
##' @references
##'
##' Hothorn, Torsten, Lisa Möst, and Peter Bühlmann. "Most Likely
##'   Transformations."  Scandinavian Journal of Statistics 45, no. 1 (March
##'   2018): 110–34.  <doi:10.1111/sjos.12291>
##'
##' Tamasi, Balint, and Torsten Hothorn. "tramME: Mixed-Effects Transformation
##'   Models Using Template Model Builder." The R Journal 13, no. 2 (2021):
##'   398–418. <doi:10.32614/RJ-2021-075>
##'
##' @inheritParams tram::tram
##' @param formula A formula describing the model. Smooth additive terms are
##'   defined the way as in \code{mgcv}, and random effects consistently with
##'   the notation used in \code{lme4}.
##' @param tram Parameter vector for the \code{tram} model type.
##' @param call The original function call (to be passed from the wrapper).
##' @param ctm A model object of the \code{ctm} class that descibes the
##'   fixed-effects part of the \code{tramME} model.
##' @param smooth A \code{tramME_smooth} object that describes the smooth
##'   additive elements of the \code{tramME} model.
##' @param negative Logical; if \code{TRUE}, the model is parameterized with
##'   negative coefficinets for the elements of the linear predictor.
##' @param silent Logical. Make \pkg{TMB} functionality silent.
##' @param resid Logical. If \code{TRUE}, the score residuals are also calculated.
##'   This comes with some performance cost.
##' @param do_update Logical. If \code{TRUE}, the model is set up so that the weights and the
##'   offsets are updateable. This comes with some performance cost.
##' @param estinit Logical. Estimate a vector of initial values for the fixed effects parameters
##'   from a (fixed effects only) mlt model
##' @param initpar Named list of initial parameter values, if \code{NULL}, it is ignored
##' @inheritParams mlt::mlt
##' @param nofit logical, if TRUE, creates the model object, but does not run the optimization
##' @param control list with controls for optimization
##' @param ... Optional arguments to \code{\link[tram]{tram}}
##' @importFrom stats na.omit model.offset model.weights
##' @export
tramME <- function(formula, data, subset, weights, offset, na.action,
                   tram = NULL, call = NULL,
                   ctm = NULL, smooth = NULL, negative = NULL,
                   silent = TRUE, resid = FALSE, do_update = FALSE,
                   estinit = TRUE, initpar = NULL,
                   fixed = NULL, nofit = FALSE,
                   control = optim_control(), ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- NULL
  args <- get_names(args, environment(formula)) ## XXX
  mod <- do.call("tramME_model", args = args)
  if (is.null(call <- substitute(call))) call <- cl

  ## -- Create model.frame
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"),
             names(call), 0L)
  fc <- call[c(1L, m)]
  ## -- NOTE: fake a tramME object and use model.frame.tramME
  fc$formula <- structure(list(model = mod), class = "tramME")
  ## --
  fc[[1L]] <- quote(model.frame)
  dat <- eval(fc, parent.frame())

  ## -- sanitize initial parameter settings
  if (is.null(mod$ranef) || nofit || !is.null(initpar)) {
    estinit <- FALSE
  }

  ## Additional parameter constraints for certain model types
  ## (from tram::Survreg & tram::Aareg)
  if (isTRUE(sub("ME$", "", tram) == "Survreg")) {
    cf <- coef(mod$ctm)
    dist <- list(...)$dist
    scale <- list(...)$scale
    scalecf <- grep(names(dat)[1], names(cf), fixed = TRUE)
    if (dist == "exponential")
      scale <- 1
    if (dist == "rayleigh")
      scale <- 0.5
    if (scale > 0) {
      fix <- rep(1 / scale, length(scalecf))
      names(fix) <- names(cf)[scalecf]
      fixed <- c(fixed, fix)
    }
  }

  if (isTRUE(sub("ME$", "", tram) == "Aareg")) {
    cf <- coef(mod$ctm)
    nm <- names(cf)
    nm <- nm[grep("Bs1", nm)]
    fix <- numeric(length(nm))
    names(fix) <- nm
    fixed <- c(fixed, fix)
  }

  ## -- create terms required by tramTMB
  ## NOTE: fixed can contain elements that don't enter the FE only
  ## mlt model (e.g. 'fixed' random effects)
  ## to avoid errors, remove these from fixed temporarily
  fixed2 <- fixed[names(fixed) %in% names(coef(mod$ctm))]
  mmlt <- mlt::mlt(mod$ctm, data = dat, offset = model.offset(dat),
                   weights = model.weights(dat),
                   fixed = fixed2, dofit = estinit)
  fe <- fe_terms(mmlt)
  re <- re_terms(mod$ranef, dat, mod$negative)
  sm <- sm_terms(mod$smooth, dat, mod$negative)

  param <- .param(fe, re, sm, fixed)

  inp <- tramTMB_inputs(mod, fe, re, sm, data = dat, param = param,
                        initpar = initpar)

  ## -- create the tramTMB object
  obj <- tramTMB(inp$data, inp$parameters, inp$constraint, inp$negative,
                 map = inp$map, resid = resid, do_update = do_update,
                 silent = silent)

  ## -- model fitting
  if (!nofit) {
    if (is.null(initpar) && !estinit) {
      par <- .optim_start(obj, resp = dat[[1]])
    } else
      par <- NULL

    opt <- optim_tramTMB(obj, par = par,
                         method = control$method, control = control$control,
                         trace = control$trace, ntry = control$ntry,
                         scale = control$scale, ok_warnings = control$ok_warnings)
    param <- .upd_param(param, obj)
  } else {
    opt <- NULL
  }
  structure(list(call = call, model = mod, data = dat, tmb_obj = obj, opt = opt,
                 param = param), class = "tramME")
}

