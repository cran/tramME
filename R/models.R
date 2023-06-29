##' Mixed-effects version of \code{\link[tram]{Coxph}}
##' @inheritParams LmME
##' @inheritParams tram::Coxph
##' @return A CoxphME object.
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


##' Mixed-effects version of \code{\link[tram]{Colr}}
##' @inheritParams LmME
##' @inheritParams tram::Colr
##' @return A ColrME object.
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


##' Mixed-effects version of \code{\link[tram]{BoxCox}}
##' @inheritParams LmME
##' @inheritParams tram::BoxCox
##' @return A BoxCoxME object.
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


##' Mixed-effects version of \code{\link[tram]{Lehmann}}
##' @inheritParams LmME
##' @inheritParams tram::Lehmann
##' @return A LehmannME object.
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


##' Mixed-effects version of \code{\link[tram]{Polr}}
##' @inheritParams LmME
##' @inheritParams tram::Polr
##' @return A PolrME object.
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
  nms <- names(which(sapply(args, is.name)))
  args[nms] <- mget(nms, env)
  args
}

##' General function to define and fit \code{tramME} models
##'
##' @details
##'
##' The specific model functions (\code{\link[tramME]{LmME}},
##' \code{\link[tramME]{BoxCoxME}}, \code{\link[tramME]{ColrME}}, etc.) are
##' wrappers around this function.
##'
##' @section Warning:
##'
##' Typically, the \code{tramME} function shouldn't be called directly; it is
##'   only exported to allow the advanced users to define their \code{tramME}
##'   models in a more flexible way from their basic building blocks.
##'
##' @inheritParams LmME
##' @param tram Parameter vector for the \code{tram} model type.
##' @param call The original function call (to be passed from the wrapper).
##' @param ctm A model object of the \code{ctm} class that descibes the
##'   fixed-effects part of the \code{tramME} model.
##' @param smooth A \code{tramME_smooth} object that describes the smooth
##'   additive elements of the \code{tramME} model.
##' @param negative Logical; if \code{TRUE}, the model is parameterized with
##'   negative coefficinets for the elements of the linear predictor.
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
  args <- get_names(args, environment()) ## XXX
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
                         scale = control$scale)
    param <- .upd_param(param, obj)
  } else {
    opt <- NULL
  }
  structure(list(call = call, model = mod, data = dat, tmb_obj = obj, opt = opt,
                 param = param), class = "tramME")
}

