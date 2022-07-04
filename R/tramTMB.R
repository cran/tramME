## Helper function to recursively drop terms from a fomula
## @param f Formula or call.
## @param what Things to remove from \code{f}.
## @note Adapted from \code{\link[lme4]{nobars_}}.
drop_terms <- function(f, what = c("|", "||", "s", "te", "ti", "t2")) {
  chk <- function(ff) any(sapply(what, function(x) ff == as.name(x)))

  if (!any(what %in% all.names(f))) return(f)
  if (is.call(f) && chk(f[[1]])) return(NULL)
  if (length(f) == 2) {
    if (is.null(f[[2]] <- drop_terms(f[[2]], what = what))) return(NULL)
    return(f)
  }
  f2 <- drop_terms(f[[2]], what = what)
  f3 <- drop_terms(f[[3]], what = what)
  if (is.null(f2)) return(f3)
  if (is.null(f3)) return(f2)
  f[[2]] <- f2
  f[[3]] <- f3
  return(f)
}


## Remove terms from a formula
## @param f Formula or call.
## @param what Things to remove from \code{f}.
## @param right_only Remove only from the right side of \code{f}.
## @note Adapted from \code{\link[lme4]{nobars}}.
##' @importFrom stats reformulate
remove_from_formula <- function(f, what = c("|", "||", "s", "te", "ti", "t2"),
                                right_only = TRUE) {
  env <- environment(f)
  if (is(f, "formula") && length(f) == 3) {
    nfr <- drop_terms(f[[3]], what)
    nfl <- if (!right_only) drop_terms(f[[2]], what) else f[[2]]
    if (is.null(nfl) && is.null(nfr))
      nf <- ~1
    else if (is.null(nfl))
      nf <- reformulate(deparse1(nfr))
    else if (is.null(nfr))
      nf <- reformulate("1", response = deparse1(nfl))
    else
      nf <- reformulate(deparse1(nfr), response = deparse1(nfl))
  } else {
    nf <- drop_terms(f, what)
  }
  if (is.null(nf))
    nf <- if (is(f, "formula")) ~1 else 1
  if (is(nf, "formula"))
    environment(nf) <- env
  return(nf)
}


## Create a corresponding ctm model for a tramME model
##
## Takes a tramME formula and generates the FE ctm model (model_only = TRUE)
## @param formula Model formula.
## @param mname tram(ME) model name.
## @param ... Optional arguments passed to \code{\link[tram]{tram}}
##' @export
.tramME2ctm <- function(formula, mname, ...) {
  .nosmooth <- function(trm) mgcv::interpret.gam(trm)$pf
  .nobars   <- function(trm) remove_from_formula(trm, c("|", "||"))
  fc <- match.call()
  mname <- sub("ME", "", mname) ## NOTE: remove suffix if present
  fc[[1L]] <- str2lang(paste0("tram::", mname))
  fc$formula <- .nosmooth(.nobars(formula)) ## NOTE: we have to remove the bars first
  fc$model_only <- TRUE
  fc$mname <- NULL
  ctmod <- eval(fc, parent.frame()) ## ctm with the fixed effects
  return(ctmod)
}


##' Create an object that defines a tramME_model
##'
##' There are two ways of defining tramME models:
##' \enumerate{
##'   \item A ctm model and a formula defining the random effects and smooth terms.
##'   \item A formula combining the notation of \pkg{tram}, \pkg{lme4} and \pkg{mgcv},
##'     a tram function name, and a dataset to set up the bases.
##' }
##' @param formula formula that either describes the whole model or
##'   the random effects specification. If the model contains random effects or
##'   smooth terms \code{formula} has to contain their definition in
##'   \pkg{lme4}-style and \pkg{mgcv}-style notation, respectively.
##' @inheritParams tram::tram
##' @param tram tram model name: Lm, BoxCox, Colr, Polr, Coxph, Survreg, Lehmann,
##'   Aareg, or the suffixed versions of these (e.g. ColrME). Ignored when a \code{ctm} model
##'   is also supplied.
##' @param ctm A \code{ctm} model
##' @param smooth Optional pre-defined smooth specification of the class \code{tramME_smooth}.
##'   If present, the smooth terms in the formula are ignored.
##' @param negative an optional parameter that defines whether the random effects have
##'   a positive or a negative sign in the model when the fixed effecst part is defined
##'   through a ctm
##' @param ... optional arguments passed to \pkg{tram} when the model is defined by the formula
##' @return A tramME_model object that defines the mixed effects transfromation model.
##' @note Similarly to \pkg{mlt}, the offsets and the weights are not part of the model,
##'   but they are data and they are not saved in the returned object.
##' @export
tramME_model <- function(formula = NULL, data = NULL, tram = NULL, ctm = NULL,
                         smooth = NULL, negative = NULL, ...) {
  ## -- TODO: add tensor products (t2) and remove this
  if (!is.null(formula) &&
      !identical(formula, remove_from_formula(formula, c("te", "ti", "t2"))))
    stop("Tensor product splines are not implemented yet.")
  ## --
  re <- lme4::findbars(formula)
  if (!is.null(smooth) && inherits(smooth, "tramME_smooth")) {
    sm <- smooth
  } else {
    sm <- if (is.null(formula)) NULL else mgcv::interpret.gam(formula)$smooth.spec
    sm <- if (length(sm)) sm else NULL
  }
  out <- structure(
    list(ctm = ctm, formula = formula, negative = negative, ranef = re,
         smooth = sm),
    class = "tramME_model")
  if (is.null(ctm)) {
    ## no ctm model, the mixed-effects model is defined by the formula
    stopifnot(!is.null(tram), !is.null(formula))
    fc <- match.call()
    fc[[1L]] <- quote(.tramME2ctm)
    fc$mname <- tram
    fc[c("tram", "ctm", "smooth", "negative")] <- NULL
    fc$data <- data ## NOTE: needed for passing the inputs correctly
    out$ctm <- eval(fc, parent.frame())
  } else {
    ## There is a ctm, that describes the FE part,
    ## formula could describe the RE par.
    stopifnot(inherits(ctm, "ctm"))
    out$formula <- fake_formula(structure(list(model = out), class = "tramME"))
  }
  if (is.null(negative)) {
    out$negative <- .mod_negative(ctm, tram)
  }
  return(out)
}


## Helper function to figure out if negative = TRUE in a given model
## @inheritParams tramME_model
.mod_negative <- function(ctm, tram = NULL) {
  if (is.null(tram)) {
    if (!is.null(ctm$bases$shifting)) {
      neg <- get("negative", envir = environment(ctm$bases$shifting),
                 inherits = FALSE)
    } else {
      neg <- mget("negative", envir = environment(ctm$bases$interacting),
                  ifnotfound = FALSE)$negative
    }
    stopifnot(!is.null(neg))
    return(neg)
  } else {
    tram <- sub("ME$", "", tram)
    return(switch(tram, Coxph = , Aareg = , Colr = FALSE,
                  Survreg = , Polr = , Lm = , BoxCox = ,
                  Lehmann = TRUE,
                  stop("Unknown model type")))
  }
}

## Create fixed effects data and initial parameters
## @param mod a mlt model
fe_terms <- function(mod) {
  par <- coef(mod)
  ## -- data
  dat <- mget(c("iY", "eY", "offset"),
    envir = environment(mod$logliki), ifnotfound = list(NULL), inherits = TRUE)
  out <- list()
  ## === Constraints & parameter types
  if (!is.null(dat$eY)) {
    out$constr <- attr(dat$eY$Y, "constraint")
    assign <- attr(dat$eY$Y, "Assign")
  } else {
    out$constr <- attr(dat$iY$Yleft, "constraint")
    assign <- attr(dat$iY$Yleft, "Assign")
  }
  out$pargroup <- apply(assign, 2, function(x) {
    if (any(grepl("shifting", x))) {
      out <- "shift"
    } else {
      out <- "baseline"
    }
    out
  })
  ## === Setting up blank values
  if (is.null(dat$iY)) {
    nm <- colnames(dat$eY$Y)
    dat$iY$Yleft <- dat$iY$Yright <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$iY$Yleft) <- colnames(dat$iY$Yright) <- nm
    dat$iY$which <- integer(0)
    dat$iY$trunc$left <- dat$iY$trunc$right <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$iY$trunc$left) <- colnames(dat$iY$trunc$right) <- nm
  }
  if (is.null(dat$eY)) {
    nm <- colnames(dat$iY$Yleft)
    dat$eY$Y <- dat$eY$Yprime <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$eY$Y) <- colnames(dat$eY$Yprime) <- nm
    dat$eY$which <- integer(0)
    dat$eY$trunc$left <- dat$eY$trunc$right <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$eY$trunc$left) <- colnames(dat$eY$trunc$right) <- nm
  }
  ## === Censoring
  idxr <- which(is.finite(dat$iY$Yleft[, 1]) & !is.finite(dat$iY$Yright[, 1]))
  idxl <- which(!is.finite(dat$iY$Yleft[, 1]) & is.finite(dat$iY$Yright[, 1]))
  idxi <- which(is.finite(dat$iY$Yleft[, 1]) & is.finite(dat$iY$Yright[, 1]))
  out$censl <- list(ay = dat$iY$Yright[idxl, , drop = FALSE],
                    which = dat$iY$which[idxl])
  out$censr <- list(ay = dat$iY$Yleft[idxr, , drop = FALSE],
                    which = dat$iY$which[idxr])
  out$censi <- list(ayl = dat$iY$Yleft[idxi, , drop = FALSE],
                    ayr = dat$iY$Yright[idxi, , drop = FALSE],
                    which = dat$iY$which[idxi])
  ## === Exact observations
  out$exact <- list(ay = dat$eY$Y, aypr = dat$eY$Yprime, which = dat$eY$which)
  ## === Offsets, weights, error distribution, etc
  out$offset <- dat$offset
  ## out$weights <- mod$weights
  ## out$negative <- mod$negative
  out$errdistr <- mod$todistr$name
  ## === Truncation
  nm <- colnames(dat$iY$Yleft)
  if (is.null(dat$eY$trunc$left)) {
    dat$eY$trunc$left <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$eY$trunc$left) <- nm
  }
  if (is.null(dat$eY$trunc$right)) {
    dat$eY$trunc$right <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$eY$trunc$right) <- nm
  }
  if (is.null(dat$iY$trunc$left)) {
    dat$iY$trunc$left <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$iY$trunc$left) <- nm
  }
  if (is.null(dat$iY$trunc$right)) {
    dat$iY$trunc$right <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$iY$trunc$right) <- nm
  }
  ## Left
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) &
                !is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) &
                !is.finite(dat$eY$trunc$right[, 1]))
  out$truncl <- list(ay = rbind(dat$iY$trunc$left[idxi, , drop = FALSE],
                                dat$eY$trunc$left[idxe, , drop = FALSE]),
                     which = c(dat$iY$which[idxi], dat$eY$which[idxe]))
  ## Right
  idxi <- which(!is.finite(dat$iY$trunc$left[, 1]) &
                is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(!is.finite(dat$eY$trunc$left[, 1]) &
                is.finite(dat$eY$trunc$right[, 1]))
  out$truncr <- list(ay = rbind(dat$iY$trunc$right[idxi, , drop = FALSE],
                                dat$eY$trunc$right[idxe, , drop = FALSE]),
                     which = c(dat$iY$which[idxi], dat$eY$which[idxe]))
  ## Interval
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) &
                is.finite(dat$iY$trunc$right[, 1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) &
                is.finite(dat$eY$trunc$right[, 1]))
  out$trunci <- list(ayl = rbind(dat$iY$trunc$left[idxi, , drop = FALSE],
                                 dat$eY$trunc$left[idxe, , drop = FALSE]),
                     ayr = rbind(dat$iY$trunc$right[idxi, , drop = FALSE],
                                 dat$eY$trunc$right[idxe, , drop = FALSE]),
                     which = c(dat$iY$which[idxi], dat$eY$which[idxe]))
  ## -- parameters
  out$offset <- NULL ## FIXME: check this
  out$beta <- par
  return(out)
}

## Create random effects data and initial paramaters
## @param ranef a list of random effects formulas from \code{\link[lme4]{findbars}}
## @param data data.frame containing the variables of the model
## @param negative logical value that indicates whether the random effects have
##   a negative sign
## @param drop.unused.levels
## @return A list containing data and parameter values to be used in the TMB model.
re_terms <- function(ranef, data, negative, drop.unused.levels = TRUE) {
  if (is.null(ranef)) {
    out <- list(Zt = Matrix::Matrix(0, nrow = 0, ncol = nrow(data), doDiag = FALSE),
                termsize = integer(0), blocksize = integer(0),
                ui = Matrix::Matrix(0, nrow = 0, ncol = 0),
                ci = numeric(0),
                gamma = numeric(0), theta = numeric(0))
  } else {
    rt <- lme4::mkReTrms(ranef, data, drop.unused.levels = drop.unused.levels)
    out <- list()
    out$Zt <- rt$Zt
    if (negative) out$Zt <- -out$Zt
    ## --- Structure of the RE covariance matrix
    out$termsize <- sapply(rt$Ztlist, NROW)
    out$blocksize <- sapply(rt$cnms, length)
    out$npar <- sum(out$blocksize * (out$blocksize + 1) / 2)
    out$names <- rt$cnms
    out$levels <- lapply(rt$flist, levels)
    ## --- Constraints
    out$ui <- Matrix::Diagonal(out$npar)
    out$ci <- rep(-Inf, out$npar)
    out$gamma <- rep(0, nrow(out$Zt))
    out$theta <- rep(0, out$npar)
  }
  return(out)
}

## Create data and initial parameter values required for the smooth terms
##
## If \code{smooth} is of \code{tramME_smooth} class, it has an additional data attribute
## that contains the observations to be used to set up the smooth terms.
## @param smooth a list of smooth terms from \code{\link[mgcv]{interpret.gam}} or
##   an object of \code{tramME_smooth} class.
## @param data data.frame containing the variables of the model
## @param negative logical value that indicates whether the random effects have
##   a negative sign
## @return A list containing data and parameter values to be used in the TMB model.
sm_terms <- function(smooth, data, negative) {
  out <- list()
  if (is.null(smooth)) {
    out$X <- matrix(0, nrow = nrow(data), ncol = 0)
    out$Z <- Matrix::Matrix(0, nrow = nrow(data), ncol = 0, doDiag = FALSE)
    out$re_dims <- numeric(0)
  } else {
    if (inherits(smooth, "tramME_smooth")) {
      sm <- .sm_data(smooth, data = attr(smooth, "data"), newdata = data)
    } else {
      sm <- .sm_data(smooth, data)
    }
    out$X <- do.call("cbind", sm$X)
    if (length(sm$Z)) {
      out$Z <- as(do.call("cbind", sm$Z), "dgTMatrix")
      out$re_dims <- sapply(sm$Z, ncol)
    } else {
      out$Z <- Matrix::Matrix(0, nrow = nrow(out$X), ncol = 0,
                              doDiag = FALSE)
      out$re_dims <- numeric(0)
    }
  }
  if (negative) {
    out$X <- -out$X
    out$Z <- -out$Z
  }
  out$uiX <- Matrix::Diagonal(ncol(out$X))
  out$ciX <- rep(-Inf, ncol(out$X))
  out$uiZ <- Matrix::Diagonal(length(out$re_dims))
  out$ciZ <- rep(-Inf, length(out$re_dims))
  out$beta <- rep(0, ncol(out$X))
  out$theta <- rep(0, length(out$re_dims))
  out$gamma <- rep(0, ncol(out$Z))
  return(out)
}


## Create inputs of the tramTMB model
##
## The function generates random values for \code{NA} values in the list of initial parameter
## values.
## @param model a list describing the structure of the model as returned by \code{tramME_model}
## @param ft fixed effects terms as returned by the function \code{fe_terms}
## @param rt random effects terms as returned by the function \code{re_terms}
## @param st smooth terms as returned by the function \code{sm_terms}
## @param data model frame containing offsets and weights
## @param fixed optional named vector with \code{c("name" = value)} pairs of
##   parameters to be fixed
## @param param optional named list of initial parameter values
## @return A list with data matrices and initial parameter values
##' @importFrom stats model.offset model.weights
## FIXME: update man
tramTMB_inputs <- function(model, ft, rt, st, data, param, initpar= NULL) {
  os <- model.offset(data) ## NOTE: get offset and weights from the model frame
  if (is.null(os)) os <- rep(0, nrow(data))
  we <- model.weights(data)
  if (is.null(we)) we <- rep(1, nrow(data))
  errdist <-  which(model$ctm$todistr$name ==
    c("normal", "logistic", "minimum extreme value", "maximum extreme value", "exponential")) - 1L
  bl <- ft$pargroup == "baseline"
  X <- rbind(ft$exact$ay[, !bl, drop = FALSE], ft$censl$ay[, !bl, drop = FALSE],
             ft$censr$ay[, !bl, drop = FALSE], ft$censi$ayl[, !bl, drop = FALSE])
  idx <- c(ft$exact$which, ft$censl$which, ft$censr$which, ft$censi$which)
  X <- X[order(idx), ]
  out <- list(
    data = list(
      errdist = rep(errdist,
        length(ft$censl$which) + length(ft$censr$which) +
        length(ft$censi$which) + length(ft$exact$which)),
      Yl = ft$censl$ay[, bl, drop = FALSE], Yr = ft$censr$ay[, bl, drop = FALSE],
      Yil = ft$censi$ayl[, bl, drop = FALSE], Yir = ft$censi$ayr[, bl, drop = FALSE],
      Ye = ft$exact$ay[, bl, drop = FALSE], Yeprime = ft$exact$aypr[, bl, drop = FALSE],
      whichl = ft$censl$which - 1L, whichr = ft$censr$which - 1L,
      whichi = ft$censi$which - 1L, whiche = ft$exact$which - 1L,
      Ytl = ft$truncl$ay[, bl, drop = FALSE], Ytr = ft$truncr$ay[, bl, drop = FALSE],
      Ytil = ft$trunci$ayl[, bl, drop = FALSE], Ytir = ft$trunci$ayr[, bl, drop = FALSE],
      whichtl = ft$truncl$which - 1L, whichtr = ft$truncr$which - 1L,
      whichti = ft$trunci$which - 1L,
      X = cbind(X, st$X),
      Z = cbind(Matrix::t(rt$Zt), st$Z),
      re_termsize = c(rt$termsize, st$re_dims),
      re_blocksize = c(rt$blocksize, rep(1L, length(st$re_dims))),
      offset = os, weights = we
    ),
    parameters = list(beta0 = ft$beta[bl], beta = c(ft$beta[!bl], st$beta),
      gamma = c(rt$gamma, st$gamma), theta = c(rt$theta, st$theta)),
    constraint = list(ui = as.matrix(Matrix::bdiag(ft$constr$ui, st$uiX, rt$ui, st$uiZ)),
                      ci = c(ft$constr$ci, st$ciX, rt$ci, st$ciZ)),
    negative = model$negative
  )
  if (is.list(initpar)) {
    if (!is.null(initpar$beta)) {
      ## NOTE: initpar is expected to have the same parameter structure as a tramME
      ## (list(beta, theta)), while in tramTMB, we split up beta to beta0 (baseline & interacting)
      ## beta (shifting FE) vectors. Convert from the fist to the second
      initpar$beta0 <- initpar$beta[bl]
      initpar$beta <- initpar$beta[!bl]
    }
    for (n in names(initpar)) {
      stopifnot(length(out$parameters[[n]]) == length(initpar[[n]]))
      out$parameters[[n]][] <- initpar[[n]]
    }
  }
  ## -- missing values are substituted with random initial values
  out$parameters$beta0[is.na(out$parameters$beta0)] <-
    sort(runif(sum(is.na(out$parameters$beta0))))
  out$parameters$beta[is.na(out$parameters$beta)] <-
    runif(sum(is.na(out$parameters$beta)))
  out$parameters$gamma[is.na(out$parameters$gamma)] <-
    runif(sum(is.na(out$parameters$gamma)))
  out$parameters$theta[is.na(out$parameters$theta)] <-
    runif(sum(is.na(out$parameters$theta)))
  ## -- convert fixed parameters as input to TMB (map)
  mp <- list()
  ## beta0 & beta
  bl <- attr(param$beta, "type") == "bl"
  fx <- attr(param$beta, "fixed")
  fx0 <- fx[bl]
  fx1 <- fx[!bl]
  if (any(fx0)) {
    mp$beta0 <- rep(NA, length(fx0))
    mp$beta0[!fx0] <- seq(sum(!fx0))
    mp$beta0 <- as.factor(mp$beta0)
  }
  if (any(fx1)) {
    mp$beta <- rep(NA, length(fx1))
    mp$beta[!fx1] <- seq(sum(!fx1))
    mp$beta <- as.factor(mp$beta)
  }
  for (n in c("theta", "gamma")) {
    fx <- attr(param[[n]], "fixed")
    if (any(fx)) {
      out$parameters[[n]][fx] <- param[[n]][fx]
      mp[[n]] <- rep(NA, length(fx))
      mp[[n]][!fx] <- seq(sum(!fx))
      mp[[n]] <- as.factor(mp[[n]])
    }
  }
  out$map <- mp
  return(out)
}


##' Create a tramTMB object
##'
##' @useDynLib tramME, .registration = TRUE
##' @param constraint list describing the constarints on the parameters
##' @param negative logical, whether the model is parameterized with negative values
##' @param map same as map argument of \code{TMB::MakeADFun}
##' @inheritParams TMB::MakeADFun
##' @param resid logical, indicating whether the score residuals are calculated
##'   from the resulting object
##' @param do_update logical, indicating whether the model should be set up with
##'   updateable offsets and weights
##' @param check_const Logical; if \code{TRUE} check the parameter constraints before
##'   evaluating the returned functions.
##' @param no_int Logical; if \code{FALSE} skip the numerical integration step.
##' @param ... optional parameters passed to \code{TMB::MakeADFun}
##' @return A tramTMB object.
##' @importFrom utils tail
##' @note The post-estimation parameters are supplied as a part of \code{data}
##' @export
## FIXME: check_const to manual
## FIXME: no_int to manual
## FIXME: consolidate no_int & part
tramTMB <- function(data, parameters, constraint, negative, map = list(),
                    resid = FALSE, do_update = FALSE, check_const = TRUE,
                    no_int = FALSE, ...) {
  ## --- ME or FE
  random <- NULL
  if (length(parameters$gamma) > 0 && !no_int)
    random <- "gamma"
  if (no_int) check_const <- FALSE
  ## TODO: REML option by adding beta to random
  ## --- add auxiliary parameters for residuals
  if (resid) {
    nn <- length(data$offset)
    parameters$alpha0 <- rep(0, nn)
  } else {
    parameters$alpha0 <- numeric(0)
  }
  if (do_update) {
    data$do_update <- 1
  } else {
    data$do_update <- 0
  }
  ## --- Adding missing post-estimation data
  if (is.null(data$postest_scale) || data$postest_scale == 0L) {
    data$Ype <- matrix(0, nrow = 0, ncol = length(parameters$beta0))
    data$Xpe <- matrix(0, nrow = 0, ncol = length(parameters$beta))
    ## data$Zpe <- Matrix::Matrix(0, nrow = 0, ncol = length(parameters$gamma),
    ##                            doDiag = FALSE)
    data$Zpe <- as(matrix(0, nrow = 0, ncol = length(parameters$gamma)), "dgTMatrix")
    data$postest_scale <- 0
  }
  if (is.null(data$as_lm)) data$as_lm <- 0
  ## --- Evaluate separate parts of the nll
  if (is.null(data$part))  data$part <- 0
  ## NOTE: if we want to evaluate the unpenalized part or the penalty
  ## of the nll, we usually want to treat all parameters as random
  ## We are not doing an optimization in these cases, just evaluation.
  if (data$part > 0) {
    random <- NULL
    check_const <- FALSE
  }
  ## --- create the TMB object
  obj <- TMB::MakeADFun(data = data, parameters = parameters, random = random,
                        DLL = "tramME", map = map, ...)
  fn <- obj$fn
  gr <- obj$gr
  he <- obj$he
  if (resid) {
    resid_idx <- tail(seq(length(obj$par)), nn)
    out <- list(
      fn = function(par, ...) {
        ## check_par(par)
        fn(c(par, rep(0, length(resid_idx))), ...)
      },
      gr = function(par, ...) {
        gr(c(par, rep(0, length(resid_idx))), ...)[-resid_idx]
      },
      he = function(par, ...) {
        he(c(par, rep(0, length(resid_idx))), ...)[-resid_idx, -resid_idx]
      },
      resid = function(par, ...) {
        res <- gr(c(par, rep(0, length(resid_idx))), ...)[resid_idx]
        if (!negative)
          res <- -res
        return(res)
      }
    )
  } else {
    resid_idx <- NULL
    out <- list(
      fn = function(par, ...) {
        ## check_par(par)
        fn(par, ...)
      },
      gr = function(par, ...) gr(par, ...),
      he = function(par, ...) he(par, ...),
      resid = function(par, ...) {
        stop("Residuals are not calculated in this model.")
      }
    )
  }

  out <- c(out, obj[!(names(obj) %in% c("fn", "gr", "he"))])

  if (resid)
    out$par <- obj$par[-resid_idx]

  ## --- add some info to out$env
  out$env$negative <- negative
  out$env$do_update <- do_update
  out$env$resid <- resid
  out$env$resid_idx <- resid_idx
  out$env$check_const <- check_const

  ## --- adjust constraints to map
  out$env$constraint <- .constr_adj(par = parameters, constr = constraint, map = map)
  ## --- Parameter checking: constraints
  out$env$par_checked <- NULL
  stopifnot(.check_par(out, out$par)) ## check initial parameters
  class(out) <- c("tramTMB", class(out))
  rm(list = setdiff(ls(), c("fn", "gr", "he", "resid_idx", "negative", "out")))
  return(out)
}

## Helper function to check parameter constarints
## @param obj A \code{tramTMB} object
## @param par A parameter vector
## @param eps Tolearnce level
## @param ... optional arguments
.check_par <- function(obj, par, eps = 1e-7, ...) {
  ## NOTE: tolerance (eps) is implied by the default value in auglag which is also
  ## used by tram
  ## FIXME: clean this up
  if (!obj$env$check_const) return(invisible(TRUE))
  if (nrow(obj$env$constraint$ui) > 0)
    out <- all(obj$env$constraint$ui %*% par - obj$env$constraint$ci > -eps)
  else out <- TRUE
  if (isTRUE(out)) {
    obj$env$par_checked <- par ## FIXME: par_checked[] to keep names?
  }
  invisible(out)
}

## Helper function to extract formatted parameters
## @param obj A \code{tramTMB} object
## @param par A parameter vector to be formatted
## @param fixed Logical; get fixed parameters, too
## @param full Should the unused parameters also be returned?
.get_par <- function(obj, par = obj$env$par_checked, fixed = TRUE, full = FALSE) {
  res <- obj$gr(par) ## NOTE: to update last.par
  if (any(obj$env$resid)) {
    par <- c(par, rep(0, length(obj$env$resid_idx)))
    out <- obj$env$parList(x = par)
    out$alpha0 <- NULL ## remove residuals
  } else {
    out <- obj$env$parList(x = par)
  }
  if (!full) {
    nz <- sapply(out, function(x) length(x) > 0)
    out <- out[nz]
  }
  map <- obj$env$map
  if (!fixed && length(map) > 0) {
    for (n in names(map)) {
      out[[n]] <- out[[n]][!is.na(map[[n]])]
    }
  }
  out
}

## Returns an extended map list with the same structure as par
## @param map Named and possibly incomplete list of parameters that shows
##   the fixed parameters. (See \code{TMB})
## @param par A complete named list of parameters
.expand_map <- function(map, par) {
  out <- lapply(par, function(x) as.factor(seq_along(x)))
  for (n in names(map)) {
    out[[n]] <- map[[n]]
  }
  return(out)
}

## Returns the number of parameters in an object
## @param obj A \code{tramTMB} object
## @param names A character vector of names to restrict
## @param fixed Logical, count the fixed parameters, too.
.npars <- function(obj, names = NULL, fixed = TRUE) {
  pr <- .get_par(obj, fixed = fixed)
  if (!is.null(names))
    pr <- pr[names]
  length(unlist(pr))
}

## Helper function to adjust the constraints consistently with the parameter
## restrictions
## @param par list containing the vector of parameters
## @param constr list containing the the constraints
## @inheritParams TMB::MakeADFun
## @return A list with adjusted constraints.
.constr_adj <- function(par, constr, map) {
  uin <- constr$ui
  cin <- constr$ci
  stopifnot(ncol(uin) == length(unlist(par[c("beta0", "beta", "theta")])))
  if (length(map) > 0) {
    map <- .expand_map(map, par)
    ## Fixing parameter values
    par <- c(par$beta0, par$beta, par$theta)
    mp <- which(c(is.na(map$beta0), is.na(map$beta), is.na(map$theta)))
    if (length(mp) > 0) {
      cin <- cin - c(uin[, mp, drop = FALSE] %*% par[mp])
      uin <- uin[, -mp, drop = FALSE]
    }
    ## Equality of parameter values (only within the same parameter vector)
    map <- lapply(map, function(x) x[!is.na(x)])
    if (length(map$beta0) > 0)
      map$beta0 <- paste0("b0", map$beta0)
    if (length(map$beta) > 0)
      map$beta <- paste0("b", map$beta)
    if (length(map$theta) > 0)
      map$theta <- paste0("th", map$theta)
    mp <- as.factor(c(map$beta0, map$beta, map$theta))
    uin <- do.call("cbind", lapply(unique(mp), function(x) {
      rowSums(uin[, mp == x, drop = FALSE])
    }))
  }
  ## Remove all zero rows
  re <- apply(uin, 1, function(x) all(x == 0))
  uin <- uin[!re, , drop = FALSE]
  cin <- cin[!re]
  ## Eliminate -Inf constraints
  re <- is.infinite(cin)
  uin <- uin[!re, , drop = FALSE]
  cin <- cin[!re]
  ## remove duplicate rows
  mm <- unique(cbind(uin, cin))
  uin <- mm[, -ncol(mm), drop = FALSE]
  cin <- mm[, ncol(mm)]
  return(list(ui = unname(uin), ci = unname(cin)))
}


## Optimize the tramTMB object
##
## Currently only with \code{alabama::auglag} with either \code{nlminb} or \code{optim}
## in the case of constrained optimization and \code{nlminb} if there are no constraints.
## @param obj a tramTMB object
## @param par optional vector of initial parameter values
## @param method the method used by \code{alabama::auglag}
## @param control a list of control parameters
## @param trace logical, whether the trace should be printed during the optimization
## @param ntry number of restarts with perturbed initial values when not converged
## @param scale Logical, if \code{TRUE}, the fixed effects design matrices are scaled
##   to improve convergence
## @param ... optional arguments, currently not in use
##' @importFrom stats nlminb
##' @importFrom stats optim
## FIXME: optional final check, w/ pdHess, mgc
optim_tramTMB <- function(obj, par = NULL, method = "nlminb", control = list(),
                          trace = FALSE, ntry = 5, scale = TRUE, ...) {
  if (!is.null(par)) {
    if (!.check_par(obj, par))
      warning(paste("The supplied initial values do not satisfy the constraints.",
                    "Falling back to the value par_checked."))
  }
  par <- obj$env$par_checked
  stopifnot(!is.null(par))
  ## --- scale
  if (scale) {
    mp <- .expand_map(obj$env$map, .get_par(obj, fixed = TRUE))
    fix <- is.na(mp$beta0)
    X <- do.call("rbind",
      obj$env$data[c("Yr", "Yl", "Yil", "Yir", "Ye",
                     "Ytl", "Ytr", "Ytil", "Ytir")])
    sc0 <- apply(abs(X[, !fix, drop = FALSE]), 2, max)
    if (.npars(obj, "beta") > 0) {
      fix <- is.na(mp$beta)
      sc1 <- apply(abs(obj$env$data$X[, !fix, drop = FALSE]), 2, max)
      sc <- c(sc0, sc1)
    } else {
      sc <- sc0
    }
    lt1 <- sc < 1.1
    gt1 <- sc >= 1.1
    sc[gt1] <- 1 / sc[gt1]
    sc[lt1] <- 1
    sc <- c(sc, rep(1, .npars(obj, "theta", fixed = FALSE)))
    par <- par / sc
    fn <- function(par) obj$fn(sc * par)
    gr <- function(par) obj$gr(sc * par) * sc
    ui <- obj$env$constraint$ui
    ci <- obj$env$constraint$ci
    if (!is.null(ui)) {
      ui <- t(t(ui) * sc)
    }
  } else {
    fn <- obj$fn
    gr <- obj$gr
    ui <- obj$env$constraint$ui
    ci <- obj$env$constraint$ci
  }
  ## ---
  warn <- NULL
  opt_time <- system.time( ## FIXME: decide if the timing is needed
    for (i in 1:ntry) {
      opt <- withCallingHandlers(
        if (!is.null(ui) && nrow(ui) > 0) {
          try(alabama::auglag(par = par, fn = fn, gr = gr,
                hin = function(par) ui %*% par - ci, hin.jac = function(par) ui,
                control.outer = list(method = method, kkt2.check = FALSE, trace = trace),
                control.optim = control)[c("par", "convergence", "value")],
              silent = !trace)
        } else {
          switch(method,
            nlminb = {
              try({
                control$trace <- trace
                op <- nlminb(par, objective = fn, gradient = gr,
                  control = control)[c("par", "convergence", "objective")]
                op$value <- op$objective
                op$objective <- NULL
                op}, silent = !trace)
            },
            try({
                control$trace <- trace
                op <- optim(par, fn = fn, gr = gr, method = method,
                  control = control)[c("par", "convergence", "value")]
                op}, silent = !trace))
        },
        warning = function(w) {
          warn <<- append(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      if (!inherits(opt, "try-error") && opt$convergence == 0) break
      par <- .optim_start(obj, par = par)
      warn <- NULL
      warn <- paste0("Number of optimization restarts: ", i)
      if (inherits(opt, "try-error"))
        opt <- list(par = par, convergence = 1)
    }
  , gcFirst = FALSE)
  opt$time <- opt_time
  opt$warnings <- warn
  if (scale) {
    opt$par <- opt$par * sc
  }
  ## NOTE: final sanity check may be redundant
  if (opt$convergence == 0) {
    if (!.check_par(obj, opt$par))
      warning("The optimum does not satisfy the parameter constraints!")
  }
  return(opt)
}


##' Set up and control optimization parameters
##' @param method Optimization procedure.
##' @param scale Logical; if \code{TRUE} rescale the fixed effects design matrix to improve
##'   convergence.
##' @param trace Logical; print trace of the optimization.
##' @param ntry Number of restarts with new random initialization if optimization
##'   fails to converge.
##' @param ... Optional arguments passed to \code{\link[alabama]{auglag}},
##'   \code{\link[stats]{nlminb}} or \code{\link[stats]{optim}} as a list of control
##'   parameters.
##' @export
optim_control <- function(method = c("nlminb", "BFGS", "CG", "L-BFGS-B"),
                          scale = TRUE, trace = FALSE, ntry = 5, ...) {
  method <- match.arg(method)
  list(method = method, scale = scale, trace = trace,
       ntry = ntry, control = list(...))
}


## Get starting values for the fixed effects parameter vector
## NOTE: adapted from .mlt_start
##' @importFrom stats qlogis qnorm qexp lm.fit rnorm runif
.optim_start <- function(obj, par = NULL, resp = NULL) {
  stopifnot(!(is.null(par) && is.null(resp)))
  dat <- obj$env$data
  ## 1) First try: use the strategy similar to mlt
  if (is.null(par)) {
    if (inherits(resp, "response"))
      resp <- resp$approxy
    ## -- NOTE: crude weighted ECDF
    we <- dat$weights
    rwe <- round(we)
    rresp <- rep(resp, rwe)
    if (inherits(resp, "factor")) {
      ra <- rank(rresp)
    } else {
      ra <- xtfrm(rresp)
    }
    pstart <- ra / max(ra)
    pstart <- pstart[cumsum(rwe)]
    pstart <- pmax(.01, pmin(pstart, .99))
    ## -- NOTE: alternative
    ## we <- dat$weights
    ## y <- mlt::R(object = resp)
    ## pstart <- attr(y, "prob")(we)(y$approxy)
    ## pstart <- pmax(.01, pmin(pstart, .99))
    ## --
    ## -- constraints
    nb <- .npars(obj, names = c("beta0", "beta"), fixed = FALSE)
    ui <- obj$env$constraint$ui[, 1:nb, drop = FALSE]
    ci <- obj$env$constraint$ci + sqrt(.Machine$double.eps)
    ## --
    X <- matrix(0, nrow = length(resp),
                ncol = .npars(obj, names = "beta0", fixed = TRUE))
    X[dat$whiche+1, ] <- dat$Ye
    X[dat$whichl+1, ] <- dat$Yl
    X[dat$whichi+1, ] <- dat$Yil
    X[dat$whichr+1, ] <- 0
    if (.npars(obj, "beta", fixed = TRUE) > 0) {
      X <- cbind(X, dat$X[c(dat$whiche+1, dat$whichl+1, dat$whichi+1, dat$whichr+1), ])
    }
    ## -- fixed
    mp <- unlist(.expand_map(obj$env$map, .get_par(obj, fixed = TRUE))[c("beta0", "beta")])
    fix <- is.na(mp)
    os <- dat$offset
    os <- os + X[, fix, drop = FALSE] %*% unlist(.get_par(obj)[c("beta0", "beta")])[fix]
    X <- X[, !fix, drop = FALSE]
    ## --
    ed <- dat$errdist
    z <- numeric(length(pstart))
    z[ed == 0] <- qnorm(pstart[ed == 0])
    z[ed == 1] <- qlogis(pstart[ed == 1])
    z[ed == 2] <- log(-log1p(-pstart[ed == 2]))
    z[ed == 3] <- -log(-log(pstart[ed == 3]))
    z[ed == 4] <- qexp(pstart[ed == 4])
    z <- z - os

    X <- X * sqrt(we)
    z <- z * sqrt(we)
    dvec <- crossprod(X, z)
    Dmat <- crossprod(X)
    diag(Dmat) <- diag(Dmat) + 1e-08

    if (!is.null(ui) && nrow(ui) > 0) {
      bb <- try(c(coneproj::qprog(Dmat, dvec, ui, ci, msg = FALSE)$thetahat),
                 silent = TRUE)
      if (inherits(bb, "try-error")) {
        diag(Dmat) <- diag(Dmat) + 1e-3
        bb <- c(coneproj::qprog(Dmat, dvec, ui, ci, msg = FALSE)$thetahat)
      }
    } else {
      bb <- lm.fit(x = X, y = z)$coef
    }
    out <- rep(0, length(obj$par))
    out[1:nb] <- bb
  }

  ## 2) Draw a random vector of initial parameter values
  if (!is.null(par) || any(is.na(out))) {
    ui <- obj$env$constraint$ui
    ci <- obj$env$constraint$ci
    ## NOTE: use the generalized inverse to invert the constraint matrix
    uii <- tcrossprod(MASS::ginv(crossprod(ui)), ui)
    bb <- NULL
    for (i in 1:100) {
      bb <- c(uii %*% exp(rnorm(ncol(uii), mean = 0, sd = 0.5)))
      if (all(ui %*% bb > ci)) {
        break
      }
    }
    if (is.null(bb)) {
      out <- sort(runif(length(par))) ## Last resort
    } else {
      out <- bb
    }
  }

  return(out)
}


##' Variance-covariance matrix of the parameters
##' @param object A \code{tramTMB} object.
##' @param par An optional vector of parameter values.
##' @param method Method for calculating the covariance matrix.
##' @param control Optional named list of controls to be passed to the specific methods.
##' @param ... Optional arguments (ignored)
##' @importFrom stats vcov optimHess
##' @export
## FIXME: might not be needed when .Hessian is also available
vcov.tramTMB <- function(object, par = object$env$par_checked,
                         method = c("optimHess", "numDeriv", "analytical"),
                         control = list(), ...) {
  method <- match.arg(method)
  if (!.check_par(object, par))
    stop("The supplied parameter vector does not satisfy the constraints.")
  he <- switch(method,
    optimHess = optimHess(par, object$fn, object$gr, control = control),
    numDeriv = {
      if (!is.null(control$method)) {
        meth <- control$method
        control$method <- NULL
      } else {
        meth <- "Richardson"
      }
      numDeriv::jacobian(func = object$gr, x = par,
                         method = meth, method.args = control)
    },
    analytical = {
      stopifnot(is.null(object$env$random))
      object$he(par)
    })
  vc <- .robustInv(he)
  rownames(vc) <- colnames(vc) <- names(par)
  return(vc)
}

## Hessian of the negative log-likelihood function
## @param object A \code{tramTMB} object.
## @param par An optional vector of parameter values.
## @param method Method for calculating the covariance matrix.
## @param control Optional named list of controls to be passed to the specific methods.
## @param joint If \code{TRUE}, calculate joint precision.
## @param ... Optional arguments (ignored)
##' @importFrom stats optimHess
## FIXME: should I make it a proper method of tramTMB?
.Hessian <- function(object, par = object$env$par_checked,
                     method = c("optimHess", "numDeriv", "analytical"),
                     control = list(), joint = FALSE, ...) {
  method <- match.arg(method)
  if (!.check_par(object, par))
    stop("The supplied parameter vector does not satisfy the constraints.")
  he <- switch(method,
    optimHess = optimHess(par, object$fn, object$gr, control = control),
    numDeriv = {
      if (!is.null(control$method)) {
        meth <- control$method
        control$method <- NULL
      } else {
        meth <- "Richardson"
      }
      numDeriv::jacobian(func = object$gr, x = par,
                         method = meth, method.args = control)
    },
    analytical = {
      stopifnot(is.null(object$env$random))
      object$he(par)
    })
  rownames(he) <- colnames(he) <- names(object$par)
  if (joint) {
    sdr <- TMB::sdreport(object, par.fixed = par, hessian.fixed = he,
                         getJointPrecision = TRUE)
    he <- sdr$jointPrecision
  }
  return(he)
}

## Invert (a block) of the Hessian using Schur complements
## @param he The Hessian.
## @param block Index vector of the block we want to invert.
## @param ... Optional arguments passed to \code{.robustInv}
.invHess <- function(he, block = NULL, ...) {
  if (!is.null(block)) {
    h1 <- he[block, block]
    h2 <- he[-block, -block]
    h3 <- he[-block, block]
    he2 <- try(h1 - crossprod(h3, solve(h2, h3)), silent = TRUE) ## w/ crossprod?
    if (inherits(he2, "try-error")) {
      return(.robustInv(he, ...)[block, block])
    } else {
      he <- he2
    }
  }
  .robustInv(he, ...)
}

## Trying harder to invert the Hessian (same as in \code{vcov.mlt})
## @param he The Hessian matrix
## @param lam Adjustmet factor. \code{lam = 0} switches off the robust option.
## @param ret_Hess Return the adjusted Hessian
## @return The variance-covariance matrix
.robustInv <- function(he, lam = 1e-6, ret_Hess = FALSE, ...) {
  step <- 0
  while((step <- step + 1) <= 3) {
    he2 <- he + (step - 1) * lam * diag(nrow(he))
    vc <- try(solve(he2), silent = TRUE)
    if (!inherits(vc, "try-error")) {
      if (step > 1) warning("Hessian could not be inverted, an approximation is used.")
      break
    }
    if (lam == 0) break
  }
  if (inherits(vc, "try-error")) vc <- he * NaN
  if (ret_Hess) return(he2)
  return(vc)
}

##' @importFrom stats update
##' @export
## FIXME: updating map this way will very likely cause errors because constraints are adjusted
update.tramTMB <- function(object, ...) {
  argn <- intersect(union(names(formals(tramTMB)), names(formals(TMB::MakeADFun))),
                    ls(object$env))
  argn <- setdiff(argn, c("random", "DLL"))
  args <- as.list(object$env)[argn]
  newargs <- list(...)
  args[names(newargs)] <- newargs
  do.call("tramTMB", args)
}

##' Generic for copying objects that are (partly) modified in place
##' @param object An object.
##' @param ... Optional parameters.
##' @export
duplicate <- function(object, ...) {
  UseMethod("duplicate")
}

##' Create a duplicate of the tramTMB object
##' @param object A \code{tramTMB} object.
##' @param ... Optional parameters (not used).
##' @importFrom utils lsf.str
##' @export
duplicate.tramTMB <- function(object, ...) {
    unserialize(serialize(object,NULL))
}

##' Post-estimation calculations in a tramTMB model
##' @param object A \code{tramTMB} object
##' @param parameters A named list of parameter values
##' @param scale The scale on which the post-estimation calculations are done
##' @param newdata A named list with elements Y, X and Z (not all necessary)
##' @param cov Logical; If \code{TRUE}, calculate the full covariance matrix
##'   of the calculated values
##' @param as.lm Logical; reparameterize as a LMM
##' @param ... Optional arguments (ignored).
##' @importFrom stats predict
##' @export
predict.tramTMB <- function(object, newdata,
                            parameters = .get_par(object, full = TRUE),
                            scale = c("lp", "trafo"), cov = FALSE,
                            as.lm = FALSE, ...) {
  scale <- as.numeric(match(match.arg(scale), eval(formals(predict.tramTMB)$scale)))
  if (scale > 1 && as.lm)
    warning("Reparametrization is only avalable for 'lp'. as.lm is set to FALSE.")
  data <- object$env$data
  names(newdata) <- paste0(names(newdata), "pe")
  data[names(newdata)] <- newdata
  data$postest_scale <- scale
  data$as_lm <- as.numeric(as.lm)
  newobj <- update(object, data = data, parameters = parameters)
  ## -- FIXME: this part to a separate function to standardize sdr calls later
  sdr <- TMB::sdreport(newobj, getReportCovariance = cov) ## try w/ default setup
  pr <- names(sdr$value) == "pred"
  if (any(is.nan(sdr$sd[pr]))) { ## NOTE: problems with inverting the Hessian. Try a bit harder
    prf <- newobj$env$last.par
    if (!is.null(newobj$env$random)) prf <- prf[-newobj$env$random]
    he <- numDeriv::jacobian(func = newobj$gr, x = prf)
    he <- .robustInv(he, ret_Hess = TRUE)
    sdr <- TMB::sdreport(newobj, par.fixed = prf, hessian.fixed = he)
  }
  ## --
  out <- list(pred = sdr$value[pr], se = sdr$sd[pr])
  if (cov)
    out$cov <- sdr$cov[pr, pr]
  return(out)
}
