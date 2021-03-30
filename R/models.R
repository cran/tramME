##' Mixed-effects version of \code{\link[tram]{Coxph}}
##' @inheritParams LmME
##' @inheritParams tram::Coxph
##' @return A CoxphME object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Coxph
##' @export
CoxphME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                    silent = TRUE, resid = FALSE, do_update = FALSE,
                    estinit = TRUE, initpar = NULL,
                    fixed = NULL, nofit = FALSE,
                    control = optim_control(),
                    ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Coxph"
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
            class = c("CoxphME", "tramME"))
}


##' Mixed-effects version of \code{\link[tram]{Colr}}
##' @inheritParams LmME
##' @inheritParams tram::Colr
##' @return A ColrME object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Colr
##' @export
ColrME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                   silent = TRUE, resid = FALSE, do_update = FALSE,
                   estinit = TRUE, initpar = NULL,
                   fixed = NULL, nofit = FALSE,
                   control = optim_control(),
                   ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Colr"
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
            class = c("ColrME", "tramME"))
}


##' Mixed-effects version of \code{\link[tram]{BoxCox}}
##' @inheritParams LmME
##' @inheritParams tram::BoxCox
##' @return A BoxCoxME object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram BoxCox
##' @export
BoxCoxME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                     silent = TRUE, resid = FALSE, do_update = FALSE,
                     estinit = TRUE, initpar = NULL,
                     fixed = NULL, nofit = FALSE,
                     control = optim_control(),
                     ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "BoxCox"
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

    opt <- optim_tramTMB(obj, method = control$method, control = control$control,
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
            class = c("BoxCoxME", "tramME"))
}


##' Mixed-effects version of \code{\link[tram]{Lehmann}}
##' @inheritParams LmME
##' @inheritParams tram::Lehmann
##' @return A LehmannME object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Lehmann
##' @export
LehmannME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      silent = TRUE, resid = FALSE, do_update = FALSE,
                      estinit = TRUE, initpar = NULL,
                      fixed = NULL, nofit = FALSE,
                      control = optim_control(),
                      ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Lehmann"
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

    opt <- optim_tramTMB(obj, method = control$method, control = control$control,
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
            class = c("LehmannME", "tramME"))
}


##' Mixed-effects version of \code{\link[tram]{Polr}}
##' @inheritParams LmME
##' @inheritParams tram::Polr
##' @return A PolrME object.
##' @importFrom stats na.omit model.offset model.weights
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

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Polr"
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
    } else
      par <- NULL

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
            class = c("PolrME", "tramME"))
}
