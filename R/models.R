##' General ME tram model
##'
##' @param silent Logical, if TRUE, prints all tracing information.
##' @param nofit Logical, if TRUE, creates the model objects, but does not run the optimization.
##' @param optim_control List of optional arguments for the optimizer.
##' @inheritParams tram::tram
##' @importFrom Matrix Matrix
##' @importMethodsFrom Matrix t
##' @importFrom mlt mlt
##' @importFrom tram Lm Coxph BoxCox Colr Polr Survreg Lehmann
##' @useDynLib tramME, .registration = TRUE
.tramME <- function(formula, data, na.action, silent, nofit, optim_control) {
  fcall <- match.call(definition = sys.function(-1), call = sys.call(-1),
                      expand.dots = TRUE) ## matching a higher level call

  ## --- Check RE terms and return tram model if FE only
  fc <- fcall
  m <- !(names(fc) %in% c("silent", "nofit", "optim_control"))
  fc <- fc[m]
  fname <- sub("ME", "", fc[[1L]])
  fc[[1L]] <- str2lang(paste0("tram::", fname))
  fc$formula <- .nobars(fc$formula)
  if (nofit)
    fc$model_only <- TRUE
  if (identical(fc$formula, fcall$formula)) {
    warnm <- paste0("The model does not contain random effects. Falling back on tram::",
                    fname)
    warning(warnm)
    mod_fe <- eval(fc, parent.frame(n = 2L)) ## NOTE: 2 levels up
    return(mod_fe)
  }

  ## --- Get data for the extended model
  ## NOTE: this will handle subsets and na.actions properly
  fc <- fcall
  m <- match(c("formula", "data", "subset", "offset", "weights", "na.action"), names(fc),
             nomatch = 0L)
  fc <- fc[c(1L, m)]
  fc[[1L]] <- quote(tram_data)
  fc$formula <- .subbars(fcall$formula)
  trdat2 <- eval(fc, parent.frame(n = 2L))

  ## --- Define the corresponding ctm (NOTE: FE only!)
  fc <- fcall
  fc$formula <- .nobars(fcall$formula)
  fc[[1L]] <- str2lang(paste0("tram::", fname))
  fc$model_only <- TRUE
  ctrm1 <- eval(fc, parent.frame(n = 2L))

  ## --- Create a non-fitted mlt model from the FE only ctm AND the data of the
  ##     extended model
  mltm <- mlt(ctrm1, data = trdat2$mf, weights = trdat2$weights,
              offset = trdat2$offset, dofit = FALSE)
  mltdat <- .mlt_data(mltm)
  if (is.null(trdat2$weights)) {
    mltdat$weights <- rep(1, nrow(trdat2$mf))
  } else {
    mltdat$weights <- trdat2$weights
  }
  ## NOTE: hard-coded bec some edge-cases would break it
  ## changed from:
  ## mltdat$negative <- get("negative", environment(ctrm1$bases$shifting))
  mltdat$negative <- switch(fname, Lm = TRUE, Coxph = FALSE, Survreg = TRUE, Colr = FALSE,
                            Polr = TRUE, BoxCox = TRUE, Lehmann = TRUE)

  ## --- RE specification
  redat <- .re_data(fcall$formula[-2L], trdat2$mf, negative = mltdat$negative)

  ## --- Constraints
  constr <- list(ui = Matrix::bdiag(mltdat$constr$ui, redat$ui),
                 ci = c(mltdat$constr$ci, redat$ci))
  ui <- as.matrix(constr$ui)
  ci <- constr$ci

  ## --- Model structure to return
  vnm <- variable.names(ctrm1)
  model <- list(
    name = fname,
    response = list(name = ctrm1$response, basis = ctrm1$bases$response),
    ranef = list(termsize = redat$termsize, blocksize = redat$blocksize,
                 npar = redat$npar, names = redat$names, levels = redat$levels),
    fixef = list(names = vnm[vnm != ctrm1$response],
                 bases = list(shifting = ctrm1$bases$shifting,
                              interacting = ctrm1$bases$interacting),
                 npar = c(all = length(mltdat$pargroup),
                          baseline = sum(mltdat$pargroup == "baseline"),
                          shift = sum(mltdat$pargroup == "shift")),
                 parnames = names(coef(ctrm1))),
    distr = mltdat$errdistr,
    negative = mltdat$negative,
    constraint = constr
  )

  if (!nofit) {
    ## --- Model type TODO: bundle these steps into a dedicated .tramME_data function
    mt <- which(fname ==
      c("Lm", "BoxCox", "Colr", "Survreg", "Coxph", "Polr", "Lehmann")) - 1L

    ## --- Error distribution
    ed <- which(mltdat$errdistr ==
      c("normal", "logistic", "minimum extreme value", "maximum extreme value")) - 1L

    ## --- Initial parameter values TODO: to a dedicated function later
    estim <- mlt(ctrm1, data = trdat2$mf, weights = trdat2$weights,
              offset = trdat2$offset, dofit = TRUE)
    params <- list(beta = estim$coef, gamma = rep(0, nrow(redat$Zt)),
                   theta = rep(0, redat$npar))

    ## --- Set up the data input NOTE: add truncated observations later
    Z <- t(redat$Zt)
    datalist <- list(
      modtype = mt, errdist = ed,
      MMl = mltdat$censl$ay, MMr = mltdat$censr$ay,
      MMil = mltdat$censi$ayl, MMir = mltdat$censi$ayr,
      MMe = mltdat$exact$ay, MMeprime = mltdat$exact$aypr,
      Zl = Z[mltdat$censl$which, , drop = FALSE],
      Zr = Z[mltdat$censr$which, , drop = FALSE],
      Zi = Z[mltdat$censi$which, , drop = FALSE],
      Ze = Z[mltdat$exact$which, , drop = FALSE],
      re_termsize = redat$termsize, re_blocksize = redat$blocksize,
      offsetl = mltdat$offset[mltdat$censl$which], offsetr = mltdat$offset[mltdat$censr$which],
      offseti = mltdat$offset[mltdat$censi$which], offsete = mltdat$offset[mltdat$exact$which],
      weightsl = mltdat$weights[mltdat$censl$which], weightsr = mltdat$weights[mltdat$censr$which],
      weightsi = mltdat$weights[mltdat$censi$which], weightse = mltdat$weights[mltdat$exact$which])

    ## --- Model
    obj <- TMB::MakeADFun(data = datalist, parameters = params, random = "gamma",
                          DLL = "tramME", silent = silent)
    ## ===================================================================
    ## FIXME: clean up
    ## TODO: list of custom optimiziers and starting values
    ## --- Optimize
    control.outer <- optim_control$outer
    control.outer$method <- "nlminb"
    control.outer$kkt2.check <- FALSE
    control.outer$trace <- !silent
    control.optim <- optim_control$optim
    if (FALSE) {
    ##obj$env$validpar <- function(par) all(ui %*% par - ci > 0) ## NOTE: is this necessary
      opt_time <- system.time(
        opt <- alabama::auglag(par = obj$par, fn = obj$fn, gr = obj$gr,
          hin = function(par) ui %*% par - ci,
          hin.jac = function(par) ui,
          control.outer = list(method = "nlminb", kkt2.check = FALSE, trace = !silent)))
    }
    if (TRUE) {
      #obj$env$validpar <- function(par) all(ui %*% par - ci > 0) ## NOTE: is this necessary
      warn <- NULL
      opt_time <- system.time(
        opt <- withCallingHandlers(
          alabama::auglag(par = obj$par, fn = obj$fn, gr = obj$gr,
            hin = function(par) ui %*% par - ci,
            hin.jac = function(par) ui,
            control.outer = control.outer, control.optim = control.optim),
          warning = function(w) {
            warn <<- append(warn, conditionMessage(w))
            invokeRestart("muffleWarning")
          }))
    }
    if (FALSE) { ## NOTE: alternative solution, does not work
      obj$env$validpar <- function(par) all(ui %*% par - ci > 0)
      opt_time <- system.time(
        opt <- optim(par = obj$par, fn = obj$fn, gr = obj$gr, method = "BFGS",
          control = list(trace = !silent)))
    }
    ## ===================================================================
    ## --- Collect results
    param_all <- obj$env$parList(opt$par, obj$env$last.par.best)
    param_fix <- opt$par
    sdr <- TMB::sdreport(obj, par.fixed = param_fix)
  }

  ## --- Output
  out <- list(call = fcall, model = model, fitted = !nofit)
  if (!nofit) {
    out$tmb_obj <- obj
    out$data <- list(mf = trdat2$mf, offset = mltdat$offset, weights = mltdat$weights)
    out$opt <- opt
    out$opt$opt_time <- opt_time
    if (!sdr$pdHess)
      warn <- append("Non-positive definite Hessian matrix!", warn)
    out$opt$warnings <- warn
    out$opt$max_gr_comp <- max(abs(sdr$gradient.fixed))
    out$tmb_sdr <- sdr
    ## --- ML estimates of parameters to return
    cf <- param_all$beta
    names(cf) <- names(ctrm1$coef)
    rep <- obj$report(unlist(param_all))
    cv <- mapply(FUN = function(sd, cr, nm) {
      ss <- diag(sd, nrow = length(sd), ncol = length(sd))
      m <- crossprod(ss, cr) %*% ss
      colnames(m) <- rownames(m) <- nm
      m
    }, sd = rep$sd_rep, cr = rep$corr_rep, nm = redat$names,
    SIMPLIFY = FALSE)
    names(cv) <- names(redat$names)
    out$pars <- list(coef = cf, varcov = cv)
  } else {
    ## --- Blank structure of parameters
    cv <- lapply(redat$names, function(nm) {
      m <- matrix(NA, nrow = length(nm), ncol = length(nm))
      colnames(m) <- rownames(m) <- nm
      m
    })
    names(cv) <- names(redat$names)
    out$pars <- list(coef = ctrm1$coef,
                     varcov = cv)
  }
  class(out) <- "tramME"
  return(out)
}

##' ME version of tram::Colr
##' @inheritParams tram::Colr
##' @inheritParams .tramME
##' @return A ColrME object.
##' @importFrom stats na.omit
##' @export
ColrME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                   silent = TRUE, nofit = FALSE,
                   optim_control = list(outer = list(), optim = list()),
                   ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("ColrME", class(out))
  return(out)
}

##' ME version of tram::Lm
##' @inheritParams tram::Lm
##' @inheritParams .tramME
##' @return A LmME object.
##' @importFrom stats na.omit
##' @export
LmME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                 silent = TRUE, nofit = FALSE,
                 optim_control = list(outer = list(), optim = list()),
                 ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("LmME", class(out))
  return(out)
}

##' ME version of tram::Polr
##' @inheritParams tram::Polr
##' @inheritParams .tramME
##' @return A PolrME object.
##' @importFrom stats na.omit
##' @export
PolrME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                   method = c("logistic", "probit", "loglog", "cloglog"),
                   silent = TRUE, nofit = FALSE,
                   optim_control = list(outer = list(), optim = list()),
                   ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("PolrME", class(out))
  return(out)
}

##' ME version of tram::Coxph
##' @inheritParams tram::Coxph
##' @inheritParams .tramME
##' @return A CoxphME object.
##' @importFrom stats na.omit
##' @export
CoxphME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                    silent = TRUE, nofit = FALSE,
                    optim_control = list(outer = list(), optim = list()),
                    ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("CoxphME", class(out))
  return(out)
}

##' ME version of tram::BoxCox
##' @inheritParams tram::BoxCox
##' @inheritParams .tramME
##' @return A BoxCoxME object.
##' @importFrom stats na.omit
##' @export
BoxCoxME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                     silent = TRUE, nofit = FALSE,
                     optim_control = list(outer = list(), optim = list()),
                     ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("BoxCoxME", class(out))
  return(out)
}

## FIXME: Exponential and Rayleigh (scale > 0) are not compatible with the current
## implementation.
##' ME version of tram::Survreg
##' @inheritParams tram::Survreg
##' @inheritParams .tramME
##' @return A SurvregME object.
##' @importFrom stats na.omit
##' @section Warning:
##'   Fixing the scale parameter is currently not available.
##' @export
SurvregME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      dist = c("weibull", "logistic", "gaussian", "exponential",
                               "rayleigh", "loggaussian", "lognormal", "loglogistic"),
                      scale = 0,
                      silent = TRUE, nofit = FALSE,
                      optim_control = list(outer = list(), optim = list()),
                      ...) {
  dist <- match.arg(dist)
  if (scale > 0 | dist %in% c("exponential", "rayleigh"))
    stop("Fixed scale is currently not available")
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) {
    class(out) <- c("SurvregME", class(out))
    out$model$name <- c(out$model$name, dist)
  }
  return(out)
}

##' ME version of tram::Lehmann
##' @inheritParams tram::Lehmann
##' @inheritParams .tramME
##' @return A LehmannME object.
##' @importFrom stats na.omit
##' @export
LehmannME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      silent = TRUE, nofit = FALSE,
                      optim_control = list(outer = list(), optim = list()),
                      ...) {
  out <- .tramME(formula, data, na.action = na.action,
                 silent = silent, nofit = nofit, optim_control = optim_control)
  if (inherits(out, "tramME")) class(out) <- c("LehmannME", class(out))
  return(out)
}
