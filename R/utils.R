##' @import methods
NULL


##' Create RE data
##' @param rhsterm Right-hand-side term of the formula
##' @param data data
##' @param negative should the random effects terms be multiplied with -1?
.re_data <- function(rhsterm, data, negative = FALSE) {
  bars <- lme4::findbars(rhsterm)
  rt <- lme4::mkReTrms(bars, data)
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
  return(out)
}


##' Extract information from an mlt model
##' @param mod mlt model
.mlt_data <- function(mod) {
  stopifnot(inherits(mod, "mlt"))
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
  }
  if (is.null(dat$eY)) {
    nm <- colnames(dat$iY$Yleft)
    dat$eY$Y <- dat$eY$Yprime <- matrix(0, nrow = 0, ncol = length(nm))
    colnames(dat$eY$Y) <- colnames(dat$eY$Yprime) <- nm
    dat$eY$which <- integer(0)
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
  return(out)
}


##' Match parameters with parameter groups and indeces
##' @param mst model structure as in an obj$model
##' @param which parameter name or index within group
##' @param pargroup parameter group
##' @param pmatch partrial matching allowed
##' @param altpar optional alternative parametrization such as LMM (\code{altpar = "lm"})
.paridx <- function(mst, which = NULL, pargroup = "all", pmatch = FALSE, altpar = NULL) {
  fnm <- mst$fixef$parnames
  rnm <- .renames(mst$ranef)
  nfp <- mst$fixef$npar
  nrp <- mst$ranef$npar
  if (isTRUE(altpar == "lm")) { ## linear mixed model parametrization
    fnm <- c("(Intercept)", fnm[3:nfp["all"]], "(Sigma)")
    pgr <- list(fixef = 1:(nfp["all"]-1), ranef = nfp["all"] + 1:nrp,
                shift = 2:(nfp["all"]-1), all = 1:(nfp["all"]+nrp))
  } else { ## general tramME object
    pgr <- list(fixef = 1:nfp["all"], ranef = nfp["all"] + 1:nrp,
                baseline = 1:nfp["baseline"],
                all = 1:(nfp["all"]+nrp))
    if (nfp["shift"] == 0) { ## NOTE: models w/o shift parameters
      pgr$shift <- integer(0)
    } else {
      pgr$shift <- nfp["baseline"] + 1:nfp["shift"]
    }
  }
  stopifnot(pargroup %in% names(pgr))
  anm <- c(fnm, rnm)
  gnm <- anm[pgr[[pargroup]]]
  if (is.numeric(which)) {
    nm <- gnm[which]
    nm <- nm[!is.na(nm)]
  } else if (is.character(which)) {
    if (pmatch) {
      idx <- unlist(lapply(which, grep, gnm))
      nm <- gnm[idx]
    } else {
      idx <- match(which, gnm)
      nm <- gnm[idx[!is.na(idx)]]
    }
  } else if (is.null(which)) {
    nm <- gnm
  } else {
    stop("Parameters are either identified by their indices or by their names.")
  }
  out <- match(nm, anm)
  names(out) <- nm
  return(out)
}


##' Return names for random effects parameters
##' @param rst the random effect structure as saved in obj$model
.renames <- function(rst) {
  nm <- rst$names
  out <- mapply(FUN = function(g, v) {
    crn <- outer(v, v, paste, sep = ".")
    paste0(g, "|", c(v, crn[upper.tri(crn)]))
  }, g = names(nm), v = nm, SIMPLIFY = FALSE)
  unlist(out, use.names = FALSE)
}


##' My own version of lme4::subbar
##' @param term Call or formula
##' @importFrom stats as.formula
.subbars <- function(term) {
  stopifnot(length(term) >= 3)
  fc <- paste(deparse(term), collapse = "")
  fc <- gsub("\\|\\|", "+", fc)
  fc <- gsub("\\|", "+", fc)
  fc <- sub("\\~", "+", fc)
  fc <- sub("\\+", "\\~", fc)
  env <- environment(term)
  if (is.null(env))
    return(str2lang(fc))
  return(as.formula(fc, env = env))
}


##' Remove random effects terms
##' @param term Call or formula
##' @importFrom stats as.formula
.nobars <- function(term) {
  stopifnot(length(term) >= 3)
  fc <- gsub(" \\+[ ]+\\([^\\)]+[\\|]+[^\\)]+\\)", "",
             paste(deparse(term), collapse = ""))
  env <- environment(term)
  if (is.null(env))
    return(str2lang(fc))
  return(as.formula(fc, env = env))
}

##' Dummy ctm model with random effects as offsets
##' @param mst List describing the model structure
##' @param coef Coefficient vector for the dummy ctm model
##' @importFrom variables numeric_var
##' @importFrom basefun as.basis
##' @importFrom mlt ctm "coef<-"
.dummy_ctm <- function(mst, coef) {
  distr <- switch(mst$distr, "normal" = "Normal", "logistic" = "Logistic",
                  "minimum extreme value" = "MinExtrVal",
                  "maximum extreme value" = "MaxExtrVal")
  re <- numeric_var("re_")
  reb <- as.basis(~ re_ , data = re, remove_intercept = TRUE,
                  negative = mst$negative)
  if (!is.null(mst$fixef$bases$shifting)) {
    newshift <- c(fe = mst$fixef$bases$shifting, re = reb)
  } else {
    newshift <- reb
  }
  mod <- ctm(response = mst$response$basis, interacting = mst$fixef$bases$interacting,
             shifting = newshift, todistr = distr)
  coef(mod) <- c(coef, re_ = 1)
  return(mod)
}


##' Reshape unformatted random effects vectors
##' @param rst the random effect structure as saved in obj$model
##' @param x Unformatted random effects vector
.re_format <- function(rst, x) {
  stopifnot(length(x) == sum(rst$termsize))
  gr <- rep(seq_along(rst$termsize), rst$termsize)
  spp <- split(x, gr)
  out <- mapply(FUN = function(sp, bs, rnm, gnm) {
    re <- matrix(sp, ncol = bs, byrow = TRUE)
    re <- as.data.frame(re)
    colnames(re) <- rnm
    rownames(re) <- rst$levels[[gnm]]
    re
  },
  sp = spp, bs = rst$blocksize, rnm = rst$names, gnm = names(rst$names),
  SIMPLIFY = FALSE)
  names(out) <- names(rst$names)
  return(out)
}


##' Calculates the size of the random effect vector implied by the model and the data
##' @param mst list describing the model structure
##' @param data Dataset containing the required grouping factors
.re_size <- function(mst, data) {
  rb <- mst$ranef$blocksize
  uv <- sapply(strsplit(names(rb), ":"), function(vn) {
    nlevels(interaction(data[vn], drop = TRUE))
  })
  list(bsize = rb, nlev = uv)
}


##' Boilerplate parallel-handling function, copied from \code{lme4}
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


##' Check the validity of a vector of coefficients
##' @param x Vector of coefficients
##' @param mst List describing the model structure
.check_coef <- function(x, mst) {
  nfe <- mst$fixef$npar[1]
  nre <- mst$ranef$npar
  if (length(x) != nfe)
    return(FALSE)
  ci <- mst$constraint$ci
  ci <- ci[1:(length(ci)-nre)]
  ui <- mst$constraint$ui
  ui <- ui[1:(nrow(ui)-nre), 1:nfe]
  if (all(ui %*% x - ci > -sqrt(.Machine$double.eps))) return(TRUE) else return(FALSE)
}


##' Check the validity of a vector of coefficients
##' @param x List of covariance matrices
##' @param vc varcov from the tramME object
.check_varcov <- function(x, vc) {
  if (length(vc) != length(x))
    return(FALSE)
  if (any(sapply(vc, dim) != sapply(x, dim)))
    return(FALSE)
  if (!all(sapply(x, isSymmetric)))
    return(FALSE)
  if (all(sapply(x, function(x) all(eigen(x, TRUE)$values > 0)))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


##' Create name for a specific tramME model
##' @param mst model structure list
.model_name <- function(mst) {
  nm <- mst$name[1L]
  str <- if (is.null(mst$fixef$bases$interacting)) "" else "Stratified "
  mnm <- switch(nm,
    Lm = "Mixed-effects Normal Linear Regression Model",
    BoxCox = "Non-normal (Box-Cox-Type) Linear Mixed-effects Regression Model",
    Colr = "Mixed-effects Continuous Outcome Logistic Regression Model",
    Coxph = "Mixed-effects Parametric Cox Regression Model",
    Polr = {
      distr <- mst$distr
      if (distr == "normal") {
        "Mixed-effects Ordered Probit Regression Model"
      } else {
        di <- c("logistic" = "Odds",
                "minimum extreme value" = "Hazards",
                "maximum extreme value" = "Lehmann-alternative")
        paste("Proportional", di[distr], "Mixed-effects Regression Model")
      }
    },
    Survreg = {
      dist <- sub("log", "Log-", mst$name[2L])
      dist <- paste0(toupper(substring(dist, 1, 1)), substring(dist, 2))
      paste("Mixed-effects", dist, "Linear Regression Model")
    },
    Lehmann = "Mixed-effects Lehmann-alternative Linear Regression Model")
  return(paste0(str, mnm))
}
