##' Remove random effects terms from full formulas or calls
##' @param term Call or formula
##' @importFrom stats as.formula
## FIXME: also allow half formulas?
.nobars <- function(term) {
  stopifnot(length(term) >= 3)
  dp <- paste(deparse(term, width.cutoff = 500L), collapse = " ")
  rh <- strsplit(dp, "~")[[1]][2]
  lh <- strsplit(dp, "~")[[1]][1]
  rh <- gsub(" [\\+]*[ ]*\\([^\\)]+[\\|]+[^\\)]+\\)", "", rh)
  rh <- strsplit(rh, "\\+")[[1]]
  if (length(rh) > 0) {
    rh <- paste(rh[!(rh %in% c(" ", ""))], collapse = "+")
    fc <- paste(lh, "~", rh)
  } else {
    fc <- paste(lh, "~ 1")
  }
  env <- environment(term)
  if (is.null(env))
    return(str2lang(fc))
  return(as.formula(fc, env = env))
}

##' Check whether formula contains bars (RE parts)
##' @param f formula
.isbars <- function(f) {
  !identical(.nobars(f), f)
}

##' Create a corresponding ctm model for a tramME model
##'
##' Takes a tramME formula and generates the FE ctm model (model_only = TRUE)
##' @param formula Model formula.
##' @param mname tram(ME) model name.
##' @param ... Optional arguments passed to \code{\link[tram]{tram}}
.tramME2ctm <- function(formula, mname, ...) {
  fc <- match.call()
  mname <- sub("ME", "", mname) ## NOTE: remove suffix if present
  fc[[1L]] <- str2lang(paste0("tram::", mname))
  fc$formula <- .nobars(formula)
  fc$model_only <- TRUE
  fc$mname <- NULL
  ctmod <- eval(fc, parent.frame()) ## ctm with the fixed effects
  return(ctmod)
}

##' Create a dummy fromula from a ctm object
##' @param ctm A ctm model from wich the fromula is created
##' @importFrom stats variable.names
.ctm2formula <- function(ctm) {
  rv <- variable.names(ctm, "response")
  vn <- setdiff(variable.names(ctm), rv)
  eval(parse(text = paste(rv, "~", paste(vn, collapse = "+"))))
}

##' Combine a set of formulas into one (similar to nlme::asOneFormula)
##' @param formula the first formula, if it contains a response that will be the
##'   response of the resulting formula
##' @param ...  objects from which a formula can be extracted
##' @param omit parameter vector with variable names to be omitted
##' @importFrom stats as.formula
.combine_formulas <- function(formula, ..., omit = ".") {
  if (length(formula) >= 3) {
    rv <- paste(deparse(formula[[2]], width.cutoff = 500L), collapse = " ")
    formula <- formula[[3]]
    omit <- c(omit, rv) ## NOTE: remove response from the RHS
  } else {
    rv <- NULL
  }
  names <- unique(unlist(lapply(c(formula, list(...)), all.vars)))
  names <- names[is.na(match(names, omit))]
  if (length(names))
    eval(parse(text = paste(rv, "~", paste(names, collapse = "+"))))
  else
    eval(parse(text = paste(rv, "~ 1")))
}

##' Create an object that defines a tramME_model
##'
##' There are two ways of defining tramME models:
##' \enumerate{
##'   \item A ctm model and a formula defining the random effects.
##'   \item A formula combining the notation of \pkg{tram} and \pkg{lme4},
##'     a tram function name, and a dataset to set up the bases.
##' }
##' @param formula formula that either describes the whole model or
##'   the random effects specification. If the model contains random effects,
##'   \code{formula} has to contain their definition in \code{lme4}-style notation.
##' @inheritParams tram::tram
##' @param tram tram model name: Lm, BoxCox, Colr, Polr, Coxph, Survreg, Lehmann,
##'   Aareg, or the suffixed versions of these (e.g. ColrME). Ignored when a \code{ctm} model
##'   is also supplied.
##' @param ctm A \code{ctm} model
##' @param negative an optional parameter that defines whether the random effects have
##'   a positive or a negative sign in the model when the fixed effecst part is defined
##'   through a ctm
##' @param ... optional arguments passed to \pkg{tram} when the model is defined by the formula
##' @return A tramME_model object that defines the mixed effects transfromation model.
##' @note Similarly to \pkg{mlt}, the offsets and the weights are not part of the model,
##'   but they are data and they are not saved in the returned object.
##' @export
tramME_model <- function(formula = NULL, data = NULL, tram = NULL, ctm = NULL,
                         negative = NULL, ...) {
  re <- lme4::findbars(formula)
  if (is.null(ctm)) {
    ## no ctm model, the mixed-effects model is defined by the formula
    stopifnot(!is.null(tram), !is.null(formula))
    ## --- ctm
    fc <- match.call()
    fc[[1L]] <- quote(.tramME2ctm)
    fc$mname <- tram
    fc$tram <- NULL
    fc$data <- quote(data) ## NOTE: needed for passing the inputs correctly
    ctm <- eval(fc)
  } else {
    ## there is a ctm, that describes the FE part, and a formula and negative the RE part
    stopifnot(inherits(ctm, "ctm"))
    env <- environment(formula)
    formula <- .combine_formulas(.ctm2formula(ctm), formula)
    environment(formula) <- env
    class(formula) <- c(class(formula), "dummy_formula")
  }
  if (is.null(negative)) {
    negative <- .mod_negative(ctm, tram)
  }
  ## Return model information and data
  structure(list(ctm = ctm, formula = formula, negative = negative, ranef = re),
            class = "tramME_model")
}

##' Helper function to figure out if negative = TRUE in a given model
##' @inheritParams tramME_model
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

##' Create fixed effects data and initial parameters
##' @param mod a mlt model
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

##' Create random effects data and initial paramaters
##' @param ranef a list of random effects formulas from \code{\link[lme4]{findbars}}
##' @param data data.frame containing the variables of the model
##' @param negative logical value that indicates whether the random effects have
##'   a negative sign
##' @return A list containing data and parameter values to be used in the TMB model.
re_terms <- function(ranef, data, negative) {
  if (is.null(ranef)) {
    out <- list(Zt = Matrix::Matrix(0, nrow = 0, ncol = 0, doDiag = FALSE),
                termsize = integer(0), blocksize = integer(0),
                ui = Matrix::Matrix(0, nrow = 0, ncol = 0),
                ci = numeric(0),
                gamma = numeric(0), theta = numeric(0))
  } else {
    rt <- lme4::mkReTrms(ranef, data)
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


##' Create inputs of the tramTMB model
##'
##' The function generates random values for \code{NA} values in the list of initial parameter
##' values.
##' @param model a list describing the structure of the model as returned by \code{tramME_model}
##' @param ft fixed effects terms as returned by the function \code{fe_terms}
##' @param rt random effects terms as returned by the function \code{re_terms}
##' @param data model frame containing offsets and weights
##' @param param optional named list of initial parameter values
##' @importFrom stats model.offset model.weights
##' @return A list with data matrices and initial parameter values
tramTMB_inputs <- function(model, ft, rt, data, param = NULL) {
  os <- model.offset(data) ## NOTE: get offset and weights from the model frame
  if (is.null(os)) os <- rep(0, nrow(data))
  we <- model.weights(data)
  if (is.null(we)) we <- rep(1, nrow(data))
  errdist <-  which(model$ctm$todistr$name ==
    c("normal", "logistic", "minimum extreme value", "maximum extreme value", "exponential")) - 1L
  out <- list(
    data = list(
      errdist = rep(errdist,
        length(ft$censl$which) + length(ft$censr$which) +
        length(ft$censi$which) + length(ft$exact$which)),
      MMl = ft$censl$ay, MMr = ft$censr$ay,
      MMil = ft$censi$ayl, MMir = ft$censi$ayr,
      MMe = ft$exact$ay, MMeprime = ft$exact$aypr,
      whichl = ft$censl$which - 1L, whichr = ft$censr$which - 1L,
      whichi = ft$censi$which - 1L, whiche = ft$exact$which - 1L,
      MTl = ft$truncl$ay, MTr = ft$truncr$ay,
      MTil = ft$trunci$ayl, MTir = ft$trunci$ayr,
      whichtl = ft$truncl$which - 1L, whichtr = ft$truncr$which - 1L,
      whichti = ft$trunci$which - 1L,
      Z = Matrix::t(rt$Zt), re_termsize = rt$termsize, re_blocksize = rt$blocksize,
      offset = os, weights = we
    ),
    parameters = list(beta = ft$beta, gamma = rt$gamma, theta = rt$theta),
    constraint = list(ui = as.matrix(Matrix::bdiag(ft$constr$ui, rt$ui)),
                      ci = c(ft$constr$ci, rt$ci)),
    negative = model$negative
  )
  if (is.list(param)) {
    for (n in names(param)) {
      stopifnot(length(out$parameters[[n]]) == length(param[[n]]))
      out$parameters[[n]][] <- param[[n]]
    }
  }
  ## -- missing values are substituted with random initial values
  out$parameters$beta[is.na(out$parameters$beta)] <-
    sort(runif(sum(is.na(out$parameters$beta))))
  out$parameters$gamma[is.na(out$parameters$gamma)] <-
    runif(sum(is.na(out$parameters$gamma)))
  out$parameters$theta[is.na(out$parameters$theta)] <-
    runif(sum(is.na(out$parameters$theta)))
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
##' @param ... optional parameters passed to \code{TMB::MakeADFun}
##' @return A tramTMB object.
##' @importFrom utils tail
##' @export
tramTMB <- function(data, parameters, constraint, negative, map = list(),
                    resid = FALSE, do_update = FALSE, ...) {
  ## --- ME or FE
  random <- NULL
  if (length(parameters$gamma) > 0)
    random <- "gamma"
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
  ## --- create the TMB object
  obj <- TMB::MakeADFun(data = data, parameters = parameters, random = random,
                        DLL = "tramME", map = map, ...)
  ## ## --- Dummy definitions to remove notes in R CMD check
  ## f <- ADreport <- tracepar <- resid_idx <- value.best <- last.par.best <- NULL
  ## check_par <- tracemgc <- usingAtomics <- ADGrad <- reportenv <- DLL <- NULL
  ## registerFinalizer <- ff <- silent <- par_checked <- gr <- parList <- NULL
  ## par_checked <- parList <- NULL
  ## ---
  fn <- obj$fn
  gr <- obj$gr
  he <- obj$he
  if (resid) {
    resid_idx <- tail(seq(length(obj$par)), nn)
  } else {
    resid_idx <- NULL
  }
  ## ---
  if (resid) {
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

  ## --- adjust constraints to map
  out$env$constraint <- .constr_adj(par = parameters, constr = constraint, map = map)
  ## --- Parameter checking: constraints
  out$env$par_checked <- NULL
  stopifnot(.check_par(out, out$par)) ## check initial parameters
  class(out) <- c("tramTMB", class(out))
  rm(list = setdiff(ls(), c("fn", "gr", "he", "resid_idx", "negative", "out")))
  return(out)
}

##' Helper function to check parameter constarints
##' @param obj A \code{tramTMB} object
##' @param par A parameter vector
##' @param eps Tolearnce level
##' @param ... optional arguments
.check_par <- function(obj, par, eps = 1e-7, ...) {
  ## NOTE: tolerance (eps) is implied by the default value in auglag which is also
  ## used by tram
  if (nrow(obj$env$constraint$ui) > 0)
    out <- all(obj$env$constraint$ui %*% par - obj$env$constraint$ci > -eps)
  else out <- TRUE
  if (isTRUE(out)) {
    obj$env$par_checked <- par ## FIXME: par_checked[] to keep names?
  }
  invisible(out)
}

##' Helper function to extract formatted parameters
##' @param obj A \code{tramTMB} object
##' @param par A parameter vector to be formatted
##' @param fixed Logica; print fixed parameters, too
.get_par <- function(obj, par = obj$env$par_checked, fixed = TRUE) {
  res <- obj$gr(par) ## NOTE: to update last.par
  if (any(obj$env$resid)) {
    par <- c(par, rep(0, length(obj$env$resid_idx)))
    out <- obj$env$parList(x = par)
    out$alpha0 <- NULL ## remove residuals
  } else {
    out <- obj$env$parList(x = par)
  }
  nz <- sapply(out, function(x) length(x) > 0)
  out <- out[nz]
  map <- obj$env$map
  if (!fixed && length(map) > 0) {
    for (n in names(map)) {
      out[[n]] <- out[[n]][!is.na(map[[n]])]
    }
  }
  out
}

##' Helper function to adjust the constraints consistently with the parameter
##' restrictions
##' @param par list containing the vector of parameters
##' @param constr list containing the the constraints
##' @inheritParams TMB::MakeADFun
##' @return A list with adjusted constraints.
.constr_adj <- function(par, constr, map) {
  uin <- constr$ui
  cin <- constr$ci
  if (length(map) > 0) {
    if (is.null(map$beta)) map$beta <- as.factor(seq_along(par$beta))
    if (is.null(map$theta)) map$theta <- as.factor(seq_along(par$theta))
    ## Fixing parameter values
    mp <- which(c(is.na(map$beta), is.na(map$theta)))
    par <- c(par$beta, par$theta)
    if (length(mp) > 0) {
      cin <- cin - c(uin[, mp, drop = FALSE] %*% par[mp])
      uin <- uin[, -mp, drop = FALSE]
    }
    ## Equality of parameter values (only within the same parameter vector)
    map$beta <- map$beta[!is.na(map$beta)]
    map$theta <- map$theta[!is.na(map$theta)]
    if (length(map$beta) > 0)
      map$beta <- paste0("b", map$beta)
    if (length(map$theta) > 0)
      map$theta <- paste0("th", map$theta)
    mp <- as.factor(c(map$beta, map$theta))
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


##' Optimize the tramTMB object
##'
##' Currently only with \code{alabama::auglag} with either \code{nlminb} or \code{optim}
##' in the case of constrained optimization and \code{nlminb} if there are no constraints.
##' @param obj a tramTMB object
##' @param par optional vector of initial parameter values
##' @param method the method used by \code{alabama::auglag}
##' @param control a list of control parameters
##' @param trace logical, whether the trace should be printed during the optimization
##' @param ntry number of restarts with perturbed initial values when not converged
##' @param scale Logical, if \code{TRUE}, the fixed effects design matrices are scaled
##'   to improve convergence
##' @param ... optional arguments, currently not in use
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
    if (!is.null(obj$env$map$beta)) {
      fix <- is.na(obj$env$map$beta)
    } else {
      fix <- rep(FALSE, length(.get_par(obj)$beta))
    }
    X <- do.call(rbind,
                 obj$env$data[c("MMr", "MMl", "MMil", "MMir", "MMe",
                                "MTl", "MTr", "MTil", "MTir")])
    sc <- apply(abs(X[, !fix, drop = FALSE]), 2, max)
    lt1 <- sc < 1.1
    gt1 <- sc >= 1.1
    sc[gt1] <- 1 / sc[gt1]
    sc[lt1] <- 1
    sc <- c(sc, rep(1, length(.get_par(obj, fixed = FALSE)$theta)))
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
##' @importFrom stats qlogis qnorm qexp lm.fit rnorm
.optim_start <- function(obj, par = NULL, resp = NULL) {
  stopifnot(!(is.null(par) && is.null(resp)))

  ## 1) First try: use the strategy similar to mlt
  if (is.null(par)) {
    if (inherits(resp, "response"))
      resp <- resp$approxy
    ## -- NOTE: crude weighted ECDF
    we <- obj$env$data$weights
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
    ## -- FIXME: alternative
    ## we <- obj$env$data$weights
    ## y <- mlt::R(object = resp)
    ## pstart <- attr(y, "prob")(we)(y$approxy)
    ## pstart <- pmax(.01, pmin(pstart, .99))
    ## --
    ## -- constraints
    nb <- length(.get_par(obj, fixed = FALSE)$beta)
    ui <- obj$env$constraint$ui[, 1:nb, drop = FALSE]
    ci <- obj$env$constraint$ci + sqrt(.Machine$double.eps)
    ## --
    nbf <- length(.get_par(obj, fixed = TRUE)$beta) ## with fixed
    X <- matrix(0, nrow = length(resp), ncol = nbf)
    X[obj$env$data$whiche+1, ] <- obj$env$data$MMe
    X[obj$env$data$whichl+1, ] <- obj$env$data$MMl
    X[obj$env$data$whichi+1, ] <- obj$env$data$MMil
    X[obj$env$data$whichr+1, ] <- 0
    ## -- fixed
    if (!is.null(obj$env$map$beta)) {
      fix <- is.na(obj$env$map$beta)
    } else {
      fix <- rep(FALSE, length(.get_par(obj)$beta))
    }
    os <- obj$env$data$offset
    os <- os + X[, fix, drop = FALSE] %*% .get_par(obj)$beta[fix]
    X <- X[, !fix, drop = FALSE]
    ## --
    ed <- obj$env$data$errdist
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
  ## NOTE: try harder to invert the Hessian (same as in vcov.mlt)
  step <- 0
  lam <- 1e-6
  while((step <- step + 1) <= 3) {
        vc <- try(solve(he + (step - 1) * lam * diag(nrow(he))), silent = TRUE)
        if (!inherits(vc, "try-error")) {
          if (step > 1) warning("Hessian could not be inverted, an approximation is used.")
          break
        }
  }
  if (inherits(vc, "try-error")) vc <- he * NaN
  rownames(vc) <- colnames(vc) <- names(par)
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
  newobj <- object
  newobj$env <- as.environment(as.list(object$env))
  parent.env(newobj$env) <- parent.env(object$env)
  env <- as.environment(as.list(environment(newobj$fn)))
  parent.env(env) <- parent.env(environment(newobj$fn))
  nm <- c("fn", "gr", "he")
  for (n in nm) {
    environment(env[[n]]) <- newobj$env
  }
  nm <- c("fn", "gr", "he", "resid")
  for (n in nm) {
    environment(newobj[[n]]) <- env
  }
  nm <- c("report", "retape", "simulate")
  for (n in nm) {
    environment(newobj[[n]]) <- newobj$env
  }
  nm <- as.character(lsf.str(envir = newobj$env))
  for (n in nm) {
    environment(newobj$env[[n]]) <- newobj$env
  }
  ## NOTE: spHess is special (see its def in TMB)
  env2 <- environment(object$env$spHess)
  environment(newobj$env$spHess) <- env2
  return(newobj)
}
