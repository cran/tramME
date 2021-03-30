##' @import methods
NULL


##' Generator for param closure
##' @param par Named list of initial parameters (beta and theta).
##' @param fe Necessary information about fixed effects.
##' @param re Necessary infromation about random effects.
##' @param varnames Names of the variables in the model.
.gen_param <- function(par, fe, re, varnames) {

  beta <- par$beta
  names(beta) <- fe$names

  theta <- .nm2vec(par$theta, re$names)

  if (length(theta) > 0) {
    varcov <- .th2vc(theta, re$blocksize)
    varcov <- .nm2mat(varcov, re$names)
  } else {
    varcov <- matrix(0, nrow = 0, ncol = 0)
  }

  out <- list(beta = beta, theta  = theta, varcov = varcov)
  mostattributes(out) <- list(names = names(out), fe = fe, re = re,
                              varnames = varnames)
  return(out)
}

##' Convert from theta vector to vc matrix
##' @param th Vector of theta parameters (reparametrization of the covariance matrices)
##' @param rbs Vector of random effects matrix dimensions
.th2vc <- function(th, rbs) {
  if (any(is.na(th)))
    return(lapply(rbs, function(x) matrix(NA, nrow = x, ncol = x)))
  ths <- rep(seq(length(rbs)), rbs * (rbs + 1) / 2)
  tapply(th, ths, function(x) {
    n <- (sqrt(1 + 8 * length(x)) - 1) / 2
    v <- exp(x[1:n])^2
    cr <- diag(n)
    cr[upper.tri(cr)] <- x[(n+1):length(x)]
    cr <- t(cr) %*% cr
    s <- sqrt(diag(cr))
    ss <- diag(1/s, nrow = length(s), ncol = length(s))
    cr <- ss %*% cr %*% ss
    sd <- diag(sqrt(v), nrow = length(v), ncol = length(v))
    sd %*% cr %*% sd
  }, simplify = FALSE)
}

##' Convert from vc matrix to theta vector
##' @param vc Covariance matrix of random effects
##' @param rbs Vector of random effects matrix dimensions
##' @importFrom stats cov2cor
.vc2th <- function(vc, rbs) {
  if (any(is.na(vc)))
    return(rep(NA, sum(rbs * (rbs + 1) / 2)))
  out <- lapply(vc, function(x) {
    stopifnot(is.matrix(x))
    n <- nrow(x)
    stopifnot(n == ncol(x))
    sd <- sqrt(diag(x))
    cr <- cov2cor(x)
    cc <- chol(cr)
    cc <- cc %*% diag(1 / diag(cc))
    th <- c(log(sd), cc[upper.tri(cc)])
    th
  })
  unlist(out)
}

##' Add names to the vc matrix
##' @param m List of covariance matrices without names
##' @param rnms Named list of random effects names
##' @return The same list of covariance matrices with proper names
.nm2mat <- function(m, rnms) {
  m <- mapply(function(x, y) { rownames(x) <- colnames(x) <- y; x },
              x = m, y = rnms, SIMPLIFY = FALSE)
  names(m) <- names(rnms)
  m
}

##' Add names to the vector of theta parameters
##' @param v theta vector without names
##' @param rnms Named list of random effects names
##' @return The same vector with proper names
.nm2vec <- function(v, rnms) {
  rn <- mapply(FUN = function(g, v) {
      crn <- outer(v, v, paste, sep = ".")
      paste0(g, "|", c(v, crn[upper.tri(crn)]))
  }, g = names(rnms), v = rnms, SIMPLIFY = FALSE)
  rn <- unlist(rn, use.names = FALSE)
  names(v) <- rn
  v
}

##' Format random effects
##' @param x Vector of random effects
##' @param rts Vector of length of random effects terms for each grouping factor
##' @param rnms Named list of random effects names
##' @param rbs Vector of random effects matrix dimensions
##' @param rlev List of the levels of each grouping factor
.re_format <- function(x, rts, rnms, rbs, rlev) {
  stopifnot(length(x) == sum(rts))
  gr <- rep(seq_along(rts), rts)
  spp <- split(x, gr)
  out <- mapply(FUN = function(sp, bs, rnm, gnm) {
    g <- matrix(sp, ncol = bs, byrow = TRUE)
    g <- as.data.frame(g)
    colnames(g) <- rnm
    rownames(g) <- rlev[[gnm]]
    g
  },
  sp = spp, bs = rbs, rnm = rnms, gnm = names(rnms),
  SIMPLIFY = FALSE)
  names(out) <- names(rnms)
  return(out)
}

##' Get the coefficient vector
##' @param obj The tramME object
.get_cf <- function(obj) {
  obj$param$beta
}

##' Set the coefficient vector
##' @param obj The tramME object
##' @param val The vector of new values
##' @return A new list of parameters with updated beta part
.set_cf <- function(obj, val) {
  stopifnot(length(val) == length(.get_par(obj$tmb_obj, fixed = FALSE)$beta))
  if (all(!is.na(val)) && all(!is.na(obj$param$theta))) { ## NOTE: check constraints and update tmb
    if (!.check_par(obj$tmb_obj, c(val, obj$param$theta))) {
      stop(paste("The assigned parameter values do not satisfy the constraints",
                 "implied by the model.\n\t",
                 "Please check BOTH coef and varcov."))
    }
  }
  obj$param$beta <- .get_par(obj$tmb_obj, c(val, obj$param$theta))$beta ## NOTE: to handle fixed values
  return(obj$param)
}


.get_vc <- function(obj, as.theta = FALSE) {
  if (as.theta)
    obj$param$theta
  else obj$param$varcov
}

##' Set the parameters of the random effect covariance matrix
##' @param obj The tramME object
##' @param val The vector of new values
##' @param as.theta The input is given according to the reparameterization used by tramTMB
##' @return A new list of parameters with updated beta part
## FIXME: It does not support fixed parameter values atm
.set_vc <- function(obj, val, as.theta = FALSE) {
  att <- attributes(obj$param)
  if (!as.theta) {
    stopifnot(identical(lapply(val, dim), lapply(obj$param$varcov, dim)))
    stopifnot(all(sapply(val, is.pd)))
    th_ <- .vc2th(val, att$re$blocksize)
  } else {
    stopifnot(length(val) == length(.get_par(obj$tmb_obj, fixed = FALSE)$theta))
    th_ <- val
  }
  if (all(!is.na(obj$param$beta)) && all(!is.null(th_))) { ## NOTE: check constraints and update tmb
    if (!is.null(obj$tmb_obj$env$map$beta))
      b_ <- obj$param$beta[!is.na(obj$tmb_obj$env$map$beta)]
    else b_ <- obj$param$beta
    if (!.check_par(obj$tmb_obj, c(b_, th_))) {
      stop(paste("The assigned parameter values do not satisfy the constraints",
                 "implied by the model.\n\t",
                 "Please check BOTH coef and varcov."))
    }
  }
  vc_ <- obj$param$varcov
  if (as.theta) {
    obj$param$theta[] <- th_
    obj$param$varcov <- mapply(function(old, new) {old[] <- new[]; old}, ## NOTE: to keep the names
                               vc_, .th2vc(val, att$re$blocksize), SIMPLIFY = FALSE)
  } else {
    obj$param$theta[] <- th_
    obj$param$varcov <- mapply(function(old, new) {old[] <- new[]; old}, ## NOTE: to keep the names
                                vc_, val, SIMPLIFY = FALSE)
  }
  return(obj$param)
}

##' Get parameter indices of various structures
##'
##' If fixed is logical, it inidcates that the indices refer to the extended
##' parameter vector
##' pargroup = c("fixef", "ranef", "shift", "all")
##' @param obj A tramME object
##' @param fixed Logical; should the indices of fixed parameters also be returned?
##' @param pargroup Parameter group
##' @param which Parameter names or indices within groups
##' @param pmatch Is partial matching allowed for \code{which}
##' @param altpar Alternative parameterizations (currently only "lm" possible)
.idx <- function(obj, fixed = NULL, pargroup = "all", which = NULL, pmatch = FALSE,
                 altpar = NULL) {
  att <- attributes(obj$param)
  if (is.null(obj$tmb_obj$env$map$beta) || is.logical(fixed))
    fnm <- names(obj$param$beta)
  else fnm <- names(obj$param$beta)[!is.na(obj$tmb_obj$env$map$beta)]

  if (is.null(obj$tmb_obj$env$map$theta) || is.logical(fixed))
    rnm <- names(obj$param$theta)
  else rnm <- names(obj$param$theta)[!is.na(obj$tmb_obj$env$map$theta)]

  ## -- identify shift parameter names
  snm <- fnm[!grepl(att$varnames[1], fnm, fixed = TRUE)]
  snm <- snm[!grepl("(Intercept)", snm, fixed = TRUE)]

  ## -- grouped names
  if (isTRUE(altpar == "lm")) {
    fnm <- c("(Intercept)", fnm[fnm %in% snm], "(Sigma)")
    gnm <- switch(pargroup, fixef = fnm[-length(fnm)], ranef = rnm,
                  all = c(fnm, rnm),
                  shift = snm, baseline = setdiff(fnm, snm))
  } else {
    gnm <- switch(pargroup, fixef = fnm, ranef = rnm, all = c(fnm, rnm),
                  shift = snm, baseline = setdiff(fnm, snm))
  }

  if (is.numeric(which)) {
    nm <- gnm[which]
    nm <- nm[!is.na(nm)]
  } else if (is.character(which)) {
    if (pmatch) {
      i <- unlist(lapply(which, grep, gnm))
      nm <- gnm[i]
    } else {
      i <- match(which, gnm)
      nm <- gnm[i[!is.na(i)]]
    }
  } else if (is.null(which)) {
    nm <- gnm
  } else {
    stop("Parameters are either identified by their indices or by their names.")
  }
  out <- match(nm, c(fnm, rnm))
  names(out) <- nm

  if (is.logical(fixed) && !fixed) {
    if (length(bfi <- which(is.na(obj$tmb_obj$env$map$beta))))
      out <- out[match(out, bfi, nomatch = 0) == 0]
    ## if (!is.null(map$theta)) { ## FIXME: does not work with fixed theta values atm
    ##   bti <- which(is.na(map$theta))
    ##   out <- out[match(out, bti, nomatch = 0) == 0]
    ## }
  }
  out
}


##' Check positive definiteness
##' @param m a matrix
##' @return logical
##' @export
is.pd <- function(m) {
  if (!isSymmetric(m))
    return(FALSE)
  all(eigen(m, symmetric = TRUE, only.values = TRUE)$values > 0)
}


##' Boilerplate parallel-handling function, modified from \code{glmmTMB}
##' @param parallel Parallel backend
##' @param ncpus Number of cores/cpus
.parallel_default <- function(parallel = c("no","multicore","snow"), ncpus = 1L) {
  if (missing(parallel)) parallel <- getOption("profile.parallel", "no")
  parallel <- match.arg(parallel)
  do_parallel <- (parallel != "no" && ncpus > 1L)
  if (do_parallel && parallel == "multicore" && .Platform$OS.type == "windows") {
    warning("No multicore on Windows, falling back to non-parallel.")
    parallel <- "no"
    do_parallel <- FALSE
  }
  if (do_parallel && parallel == "snow") {
    warning("Snow support is not implemented yet, falling back to non-parallel.")
    parallel <- "no"
    do_parallel <- FALSE
  }
  return(list(parallel=parallel,do_parallel=do_parallel))
}


##' Generates proper model name for the tramME model
##' @param obj A \code{tramME} object.
.model_name <- function(obj) {
  nm <- sub("ME$", "", obj$call[[1L]])
  str <- if (is.null(obj$model$ctm$bases$interacting)) "" else "Stratified "
  me <- if (is.null(obj$model$ranef)) "" else "Mixed-effects "
  mnm <- switch(nm,
    Lm = paste0(me, "Normal Linear Regression Model"),
    BoxCox = paste0("Non-normal (Box-Cox-Type) Linear ", me, "Regression Model"),
    Colr = paste0(me, "Continuous Outcome Logistic Regression Model"),
    Coxph = paste0(me, "Parametric Cox Regression Model"),
    Polr = {
      distr <- obj$model$ctm$todistr$name
      if (distr == "normal") {
        paste0(me, "Ordered Probit Regression Model")
      } else {
        di <- c("logistic" = "Odds ",
                "minimum extreme value" = "Hazards ",
                "maximum extreme value" = "Reverse-time Hazards ")
        paste0("Proportional ", di[distr], me, "Regression Model")
      }
    },
    Survreg = {
      dist <- obj$call[["dist"]]
      if (is.null(dist)) dist <- "weibull"
      dist <- sub("log", "Log-", dist)
      dist <- paste0(toupper(substring(dist, 1, 1)), substring(dist, 2), " ")
      paste0(me, dist, "Linear Regression Model")
    },
    Lehmann = "Mixed-effects Lehmann-alternative Linear Regression Model",
    Aareg = "Mixed-effects Parametric Linear Aalen Regression Model")
  return(paste0(str, mnm))
}


##' Create a 'conditional' ctm model with random effects as offsets
##' @param mod A \code{ctm} model.
##' @param coef Coefficient vector for the dummy ctm model.
##' @param negative The sign of the random effect term in the corresponding
##'   \code{tramME} model.
##' @importFrom variables numeric_var
##' @importFrom basefun as.basis
##' @importFrom mlt ctm "coef<-"
.cctm <- function(mod, coef, negative = FALSE) {
  distr <- switch(mod$todistr$name, "normal" = "Normal",
                  "logistic" = "Logistic",
                  "minimum extreme value" = "MinExtrVal",
                  "maximum extreme value" = "MaxExtrVal")
  re <- numeric_var("re_")
  reb <- as.basis(~ re_ , data = re, remove_intercept = TRUE,
                  negative = negative)
  if (!is.null(mod$bases$shifting)) {
    newshift <- c(fe = mod$bases$shifting, re = reb)
  } else {
    newshift <- reb
  }
  mod <- ctm(response = mod$bases$response, interacting = mod$bases$interacting,
             shifting = newshift, todistr = distr)
  coef(mod) <- c(coef, re_ = 1)
  return(mod)
}


##' Calculates the size of the random effect vector implied by the model and the data
##' @param bs Blocksize vector as returned by the function \code{re_terms}.
##' @param data Dataset containing the required grouping factors.
.re_size <- function(bs, data) {
  uv <- sapply(strsplit(names(bs), ":"), function(vn) {
    nlevels(interaction(data[vn], drop = TRUE))
  })
  list(bsize = bs, nlev = uv)
}
