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
  cl <- match.call()
  ## cl$formula <- formula ##
  cl$call <- cl
  cl[[1L]] <- quote(tramME)
  cl$tram <- "Coxph"
  eval(cl, parent.frame())
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
  cl <- match.call()
  cl$call <- cl
  cl[[1L]] <- quote(tramME)
  cl$tram <- "Colr"
  eval(cl, parent.frame())
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
  cl <- match.call()
  cl$call <- cl
  cl[[1L]] <- quote(tramME)
  cl$tram <- "BoxCox"
  eval(cl, parent.frame())
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
  cl <- match.call()
  cl$call <- cl
  cl[[1L]] <- quote(tramME)
  cl$tram <- "Lehmann"
  eval(cl, parent.frame())
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
  cl <- match.call()
  cl$call <- cl
  clmethod <- match.arg(method)
  cl[[1L]] <- quote(tramME)
  cl$tram <- "Polr"
  eval(cl, parent.frame())
}


##' General function to define and fit tramME models
##'
##' The specific model types (\code{\link[tramME]{LmME}},
##' \code{\link[tramME]{BoxCoxME}}, \code{\link[tramME]{ColrME}}, etc.) are
##' wrappers around this function.
##'
##' @section Warning:
##'
##' You should not call directly this function. Only exported for technical reasons.
##'
##' @inheritParams LmME
##' @param tram Parameter vector for the \code{tram} model type.
##' @param call The original function call (to be passed from the wrapper).
##' @importFrom stats na.omit model.offset model.weights
##' @export
tramME <- function(formula, tram, call,
                   data, subset, weights, offset, na.action,
                   silent = TRUE, resid = FALSE, do_update = FALSE,
                   estinit = TRUE, initpar = NULL,
                   fixed = NULL, nofit = FALSE,
                   control = optim_control(), ...) {
  fc <- match.call()
  fc$call <- NULL
  fc[[1L]] <- quote(tramME_model)
  mod <- eval(fc, parent.frame())

  call <- substitute(call)

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
  if (sub("ME$", "", tram) == "Survreg") {
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

  if (sub("ME$", "", tram) == "Aareg") {
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
  ## param <- .gen_param(obj, fe, re, sm, dat, nofit)
  structure(list(call = call, model = mod, data = dat, tmb_obj = obj, opt = opt,
                 param = param),
            class = c(paste0(sub("ME$", "", tram), "ME"), "tramME"))

}

