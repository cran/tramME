##' Mixed-effects version of \code{\link[tram]{Survreg}}
##' @inheritParams tram::Survreg
##' @inheritParams LmME
##' @param silent logical, make TMB functionality silent
##' @param resid logical, Should the score residuals also be calculated?
##' @param estinit logical, estimate a vector of initial values for the fixed effects parameters
##'   from a (fixed effects only) mlt model
##' @param initpar named list of initial parameter values, if \code{NULL}, it is ignored
##' @inheritParams mlt::mlt
##' @param nofit logical, if TRUE, creates the model object, but does not run the optimization
##' @param control list with controls for optimization
##' @return A \code{SurvregME} object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Survreg
##' @export
SurvregME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      dist = c("weibull", "logistic", "gaussian", "exponential",
                               "rayleigh", "loggaussian", "lognormal", "loglogistic"),
                      scale = 0,
                      silent = TRUE, resid = FALSE, do_update = FALSE,
                      estinit = TRUE, initpar = NULL,
                      fixed = NULL, nofit = FALSE,
                      control = optim_control(),
                      ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Survreg"
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

  ## -- fixing parameter values (from tram::Survreg)
  dist <- match.arg(dist)
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

  mp <- list()
  if (!is.null(fixed)) {
    idx <- which(names(cf) %in% names(fixed))
    bb <- rep(NA, length(cf))
    bb[-idx] <- seq_along(bb[-idx])
    mp <- list(beta = as.factor(bb))
  }

  ## -- create terms required by tramTMB
  mmlt <- mlt::mlt(mod$ctm, data = dat, offset = model.offset(dat),
              weights = model.weights(dat), ## TODO: offset and weights might not be needed
              fixed = fixed, dofit = estinit)
  fe <- fe_terms(mmlt)
  re <- re_terms(mod$ranef, dat, mod$negative)
  inp <- tramTMB_inputs(mod, fe, re, dat, param = initpar)

  cf[names(fixed)] <- fixed

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
            class = c("SurvregME", "tramME"))
}


##' Extract the coefficients of the fixed effects terms of an SurvregME model.
##' @param object An \code{SurvregME} object.
##' @param as.survreg If \code{TRUE}, return the transformed coefficients as in a
##'   \code{survival::survreg} object.
##' @inheritParams coef.LmME
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' library("survival")
##' fit <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats)
##' coef(fit, as.survreg = TRUE)
##' @importFrom stats coef
##' @export
coef.SurvregME <- function(object, as.survreg = FALSE, ...) {
  coef.LmME(object, as.lm = as.survreg, ...)
}
