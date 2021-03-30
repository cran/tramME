##' Mixed-effects version of \code{\link[tram]{Aareg}}
##' @inheritParams tram::Aareg
##' @inheritParams LmME
##' @inheritParams mlt::mlt
##' @return A \code{AaregME} object.
##' @importFrom stats na.omit model.offset model.weights
##' @importFrom tram Aareg
##' @export
AaregME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                 silent = TRUE, resid = FALSE, do_update = FALSE,
                 estinit = TRUE, initpar = NULL,
                 fixed = NULL, nofit = FALSE,
                 control = optim_control(),
                 ...) {
  cl <- match.call()

  ## -- create intial model structure
  fc <- cl
  fc[[1L]] <- quote(tramME_model)
  fc$tram  <-  "Aareg"
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

  ## -- fixing parameter values (from tram::Aareg)
  nm <- names(cf)
  nm <- nm[grep("Bs1", nm)]
  fix <- numeric(length(nm))
  names(fix) <- nm
  fixed <- c(fixed, fix)

  mp <- list()
  idx <- which(names(cf) %in% names(fixed))
  bb <- rep(NA, length(cf))
  bb[-idx] <- seq_along(bb[-idx])
  mp <- list(beta = as.factor(bb))

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
            class = c("AaregME", "tramME"))
}


