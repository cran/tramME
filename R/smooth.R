## Smooth shift terms

## Calculate design matrix elements for a specific smooth
## Adapted from the original code by Simon Wood.
## @param sm a smooth object as returned from \code{\link[mgcv]smoothCon}
## @param data a data.frame
## @return A list with the same structure as \code{\link[mgcv]smoothCon}
.s2rPred <- function(sm, data) {
  re <- mgcv::smooth2random(sm, names(data), type = 2)
  ## prediction matrix for new data
  X <- mgcv::PredictMat(sm, data)
  if (sm$fixed) { ## If not penalized, we only need the FE matrix
    return(list(Xf = X))
  }
  ## transform to RE parameterization
  if (!is.null(re$trans.U))
    X <- X %*% re$trans.U
  X <- t(t(X) * re$trans.D)
  ## re-order columns according to random effect re-ordering
  X[, re$rind] <- X[, re$pen.ind != 0]
  ## re-order penalization index in same way
  pen.ind <- re$pen.ind
  pen.ind[re$rind] <- pen.ind[pen.ind > 0]
  ## start returning the object
  Xf <- X[, which(re$pen.ind == 0), drop = FALSE]
  out <- list(rand = list(), Xf = Xf)
  for (i in seq_along(re$rand)) {
    ## loop over random effect matrices
    out$rand[[i]] <- X[, which(pen.ind == i), drop = FALSE]
    attr(out$rand[[i]], "s.label") <- attr(re$rand[[i]], "s.label")
  }
  names(out$rand) <- names(re$rand)
  out
}

## Generate design matrices for smooth terms
## @param smooth list of smooth specification terms
## @param data \code{data.frame} to be used to set up smooth terms
## @param newdata optional \code{data.frame} to set up design matrices
##   for prediction
## @return A list of named lists of FE and RE design matrices for each smooth term
.sm_data <- function(smooth, data, newdata = NULL) {
    bs <- lapply(smooth, mgcv::smoothCon, data = data, absorb.cons = TRUE) ## TODO: check other options
    Xs <- list()
    Zs <- list()
    re_nms <- character(0)
    bs_nms <- character(0)
    for (i in 1:length(bs)) {
      for (j in 1:length(bs[[i]])) {
        sm <- bs[[i]][[j]]
        if (is.null(newdata)) { ## same data for X,Y as setting up the smooth
          rp <- mgcv::smooth2random(sm, names(data), 2)
        } else { ## for new data
          rp <- .s2rPred(sm, newdata)
        }
        Xf <- rp$Xf
        colnames(Xf) <- rep(sm$label, ncol(Xf))
        Xs <- c(Xs, list(Xf))
        if (!is.null(Xr <- rp$rand$Xr)) {
          colnames(Xr) <- rep(sm$label, ncol(Xr))
          Zs <- c(Zs, list(Xr))
          re_nms <- c(re_nms, sm$label)
        }
        bs_nms <- c(bs_nms, sm$label)
      }
    }
    names(Zs) <- re_nms
    names(Xs) <- bs_nms
    list(X = Xs, Z = Zs)
}

## Evaluate the smooth terms of a model.
## @param object An object.
## @param ... Optional parameters.
##' @export
smooth_terms <- function(object, ...) {
  UseMethod("smooth_terms")
}

## Create grid for plotting smooth terms
##
## @param smooth list of smooth terms
## @param data data.frame with the variables for the smooths. Used to define
##   the limits of the grid.
## @param newdata optional data.frame with the values of the necessary variables
##   at which we want to evaluate the smooths. Smooths for which the supplied information
##   in this input is incomplete will be ignored.
## @param k size of the grid (only used when newdata = NULL)
## @return List of data.frames with a grid using the variable and the smooth
##   evaluated on it.
## FIXME: are both data and newdata needed?
.make_grid <- function(smooth, data, newdata = NULL, k = 100) {
  ## NOTE: sorted unique elements of factor as factor with original levels
  fau <- function(f) sort(factor(as.character(unique(f)), levels(f)))

  out <- list()
  for (i in 1:length(smooth)) {
    sm <- smooth[[i]]
    vn <- sapply(sm$term, function(n) all.vars(str2lang(n)))
    if (is.null(newdata)) {
      kk <- floor(k ^ (1/length(vn)))
      gr <- lapply(vn, function(n) {
        vv <- data[[n]]
        if (is.factor(vv)) {
          out <- fau(vv)
        } else {
          rn <- range(vv, na.rm = TRUE)
          out <- seq(rn[1], rn[2], length.out = kk)
        }
        return(out)
      })
      gr <- expand.grid(gr, KEEP.OUT.ATTRS = FALSE)
      colnames(gr) <- vn
    } else {
      if (all(vn %in% names(newdata))) {
        gr <- newdata[vn]
      } else {
        next ## NOTE: if a term is missing, the smooth won't be evaluated
      }
    }
    attr(gr, "label") <- sm$label
    attr(gr, "term") <- sm$term
    if (sm$by != "NA") { ## by
      fn <- sm$by
      if (is.factor(data[[fn]])) { ## factor
        ## FIXME: restrict the grid to levels that are present in newdata
        for (j in 1:nlevels(data[[fn]])) {
          fl <- sort(unique(data[[fn]]))[j]
          gr[[fn]] <- fl
          out <- c(out, list(gr))
          names(out)[length(out)] <- paste0(sm$label, ":", fn, fl)
        }
      } else { ## continuous
        gr[[fn]] <- 1
        out <- c(out, list(gr))
        names(out)[length(out)] <- paste0(sm$label, ":", fn)
      }
    } else {
      out <- c(out, list(gr))
      names(out)[length(out)] <- sm$label
    }
  }
  return(out)
}

##' Extract and evaluate the smooth terms of a tramME model
##' @param object A \code{tramME} object.
##' @param k Integer, the number of points to be used to evaluate the smooth terms.
##'   Ignored when \code{newdata} is supplied.
##' @param newdata A \code{data.frame} with new values for the smooth terms.
##'   If \code{NULL}, the new data is set up based on the \code{model.frame} and
##'   \code{k}. Smooths for which the supplied information in this input is incomplete
##'   will be ignored.
##' @param ... Optional arguments. \code{as.lm} is passed through this when it is necessary.
##' @return A list of results from evaluating the smooth terms of the model.
##' @examples
##' data("mcycle", package = "MASS")
##' fit <- LmME(accel ~ s(times), data = mcycle)
##' plot(smooth_terms(fit))
##' @export
##' @aliases smooth_terms
## TODO: options for extracting SEs & covariance matrices or bands
## TODO: add a fun_x option to change the covariate scale
## right now, when s(fun(x)) is the term, the function return x vs s(fun(x))
## with fun_x, one could supply monotonic (?) functions, to extract e.g.
## fun(x) vs s(fun(x)). Make sure that plot method knows what's happening.
smooth_terms.tramME <- function(object, k = 100, newdata = NULL, ...) {
  if (is.null(as.lm <- list(...)$as.lm))
    as.lm <- FALSE
  ## if no smooth term or not fitted -> empty object
  if (is.null(object$model$smooth))
    return(structure(list(), class = c("smooth.tramME", "list")))
  b <- object$param$beta
  th <- object$param$theta
  if (any(is.na(b[attr(b, "type") == "sm"])) || any(is.na(th[attr(th, "type") == "sm"])))
    return(structure(list(), class = c("smooth.tramME", "list")))
  ## -- NOTE: Make it robuts for different RE structures and scenarios
  ## 1) To avoid generating model matrices w/ columns for unused factors.
  ## 2) Use all levels that are present in the dataset, even when grouping variable is
  ## numeric.
  ## 3) Not all combinations of levels are used in nested structures:
  ## make sure we have the same length of random effects vector as in the original dataset.
  mf <- model.frame(object, drop.unused.levels = TRUE)
  gl <- length(object$param$gamma)
  gn <- names(object$param$gamma)
  ## --
  grs <- .make_grid(object$model$smooth, mf, newdata = newdata, k = k)
  lns <- sapply(grs, nrow)
  idxs <- split(1:sum(lns), rep(1:length(lns), lns))
  Xs <- list()
  Zs <- list()
  out <- list()
  ## NOTE: to speed up calculations, generate model matrices, concatenate, and do the
  ## predict.tramTMB call in one go
  for (i in 1:length(grs)) {
    gr <- grs[[i]]
    lab <- attr(gr, "label")
    dat <- mf[rep(1, nrow(gr)), ]
    dat[names(gr)] <- gr
    mm <- model.matrix(object, data = dat, type = c("X", "Zt"), keep_sign = FALSE)
    X <- mm$X
    X[, !grepl(lab, colnames(X), fixed = TRUE)] <- 0
    ## -- XXX: see explanation in 3)
    Z1 <- Matrix::t(mm$Zt)
    idx1 <- which(grepl(lab, colnames(Z1), fixed = TRUE))
    Z <- nullTMatrix(nrow = nrow(Z1), ncol = gl)
    idx <- which(grepl(lab, gn, fixed = TRUE))
    Z[, idx] <- Z1[, idx1]
    ## --
    Xs[[i]] <- X
    Zs[[i]] <- Z
  }
  X <- do.call("rbind", Xs)
  Z <- as(do.call("rbind", Zs), "TsparseMatrix")
  pr <- predict(object$tmb_obj, newdata = list(X = X, Z = Z), scale = "lp", as.lm = as.lm)
  for (i in 1:length(grs)) {
    gr <- grs[[i]]
    lab <- attr(gr, "label")
    gr[[lab]] <- pr$pred[idxs[[i]]]
    gr$se <- pr$se[idxs[[i]]]
    out[[names(grs)[i]]] <- gr
  }
  class(out) <- c("smooth.tramME", class(out))
  return(out)
}

## Subsetting smooth.tramME objects
##' @export
"[.smooth.tramME" <- function(x, i) {
  out <- unclass(x)[i]
  class(out) <- c("smooth.tramME", class(out))
  out
}

##' Evaluate smooth terms of a \code{LmME} model.
##' @inheritParams smooth_terms.tramME
##' @param as.lm Logical; if \code{TRUE} return the rescaled values according to a LMM
##'   parametrization.
##' @return A list of results from evaluating the smooth terms of the model.
##' @examples
##' data("mcycle", package = "MASS")
##' fit <- LmME(accel ~ s(times), data = mcycle)
##' plot(smooth_terms(fit, as.lm = TRUE))
##' @export
smooth_terms.LmME <- function(object, as.lm = FALSE, k = 100, newdata = NULL, ...) {
  smooth_terms.tramME(object, k = k, newdata = newdata, as.lm = as.lm, ...)
}


##' Plot smooth terms of a tramME model.
##' @param x A \code{smooth.tramME} object.
##' @param which Select terms to be printed by their indices
##' @param col Line color for the point estimates.
##' @param fill Fill color for the confidence intervals.
##' @param trafo Monotonic transformation to be applied on the smooth terms
##' @param add Add the plot to an existing figure.
##' @param ... Optional parameters passed to the plotting functions.
##' @examples
##' data("mcycle", package = "MASS")
##' fit <- LmME(accel ~ s(times), data = mcycle)
##' plot(smooth_terms(fit, as.lm = TRUE))
##' @importFrom graphics par plot.default
##' @importFrom grDevices grey
##' @export
## TODO: same y limits
## TODO: a dedicated function to calculate CIs
plot.smooth.tramME <- function(x, which = seq_along(x), col = 1, fill = grey(0.5, 0.25),
                               trafo = I, add = FALSE, ...) {
  if (length(x) == 0)
    return(invisible(x))
  if ((n <- length(x)) > 1 && !add && length(which) > 1) {
    pp <- par(mfrow = c(nn <- floor(sqrt(n)), ceiling(n / nn)))
    on.exit(par(pp))
  }
  call <- match.call()
  for (i in which) {
    vn <- attr(x[[i]], "term")
    stopifnot(length(vn) == 1)
    vn <- all.vars(str2lang(vn))
    ln <- attr(x[[i]], "label")
    xx <- x[[i]][[vn]]
    yy <- x[[i]][[ln]]
    ci <- trafo(yy + qnorm(0.975) * x[[i]]$se %o% c(-1, 1))
    yy <- trafo(yy)
    if (!add) {
      fc <- call
      fc[c("which", "col", "fill", "trafo", "add")] <- NULL
      fc[[1L]] <- quote(plot)
      fc$x <- 0
      fc$type <- "n"
      fc$ylab <- if (is.null(fc$ylab)) names(x)[i] else fc$ylab
      fc$xlab <- if (is.null(fc$xlab)) names(x[[i]])[1] else fc$xlab
      fc$ylim <- if (is.null(fc$ylim)) range(ci) else fc$ylim
      fc$xlim <- if (is.null(fc$xlim)) range(xx) else fc$xlim
      eval(fc)
    }
    fc <- call
    fc[c(formalArgs(plot.default), "which", "col", "fill", "trafo", "add")] <- NULL
    fc[[1L]] <- quote(lines)
    fc$x <- xx
    fc$y <- yy
    fc$col <- col
    eval(fc)
    fc[[1L]] <- quote(polygon)
    fc$x <- c(xx, rev(xx))
    fc$y <- c(ci[, 1], rev(ci[, 2]))
    fc$border <- NA
    fc$col <- fill
    eval(fc)
  }
  invisible(x)
}


## Calculate the effective degrees of freedom for smooth terms
## @param object An object.
## @param ... Optional parameters.
##' @export
edf_smooth <- function(object, ...) {
  UseMethod("edf_smooth")
}

##' EDFs of smooth shift terms
##'
##' Returns an estimate of effective degrees of freedom associated with each
##' smooth term.
##'
##' @details
##'
##' The EDFs are calculated by summing up the elements of
##'   \deqn{diag(V_{\vartheta}I)}{diag(VI)} term-by-term.
##'   \eqn{V_{\vartheta}}{V} is the joint covariance matrix of fixed and random
##'   parameters (the inverse of the joint precision, i.e., Hessian of the
##'   negative log-likelihood), and \eqn{I} is the joint precision of the
##'   unpenalized negative log-likelihood function. See Wood et al. (2016) or
##'   Wood (2017, Chapter 6) for references.
##'
##' @references
##'
##' Wood, Simon N., Natalya Pya, and Benjamin Saefken (2016).  "Smoothing
##'   Parameter and Model Selection for General Smooth Models."  Journal of the
##'   American Statistical Association 111, <doi:10.1080/01621459.2016.1180986>
##'
##' Wood, Simon N. (2017). Generalized Additive Models: An Introduction with R.
##'   Second edition. Chapman & Hall/CRC Texts in Statistical Science.
##'
##' @param object A \code{tramME} object.
##' @param ... Optional arguments passed to the Hessian calculations.
##' @return A named vector with the edf values.
##' @examples
##' data("mcycle", package = "MASS")
##' fit <- LmME(accel ~ s(times), data = mcycle)
##' edf_smooth(fit)
##' @export
##' @aliases edf_smooth
## TODO: The current approach is very slow with largeish datasets (~10k)
## maybe not integrating out the random effects is not the way to go
## (large, non-sparse Hessians) especially when there are many
## This will force us to give up 'analytical'
## Other idea: fix everything that is not interesting in the auxiliary models
## TODO: extend this to other terms (eg. REs) -- only after sorting out the
## performance issue
##' @importFrom Matrix solve rowSums
edf_smooth.tramME <- function(object, ...) {
  if (is.null(object$model$smooth)) return(NULL)
  ## Get indices
  pg <- c("beta", "gamma")
  pt <- unlist(sapply(pg, function(x) attr(object$param[[x]], "type")))
  pn <- unlist(sapply(pg, function(x) names(object$param[[x]])))
  idx <- which(pt == "sm")
  names(idx) <- sub("\\|.*$", "", pn[idx])
  ## Get penalized (neg log-posterior) Hessian
  newobj <- update(object$tmb_obj, no_int = TRUE)
  args <- list(...)
  args$object <- newobj
  args$par <- c(object$param$beta, object$param$gamma, object$param$theta)
  args$joint <- FALSE
  args$method <- "analytical"
  args$sparse <- TRUE
  Hp <- do.call("Hess", args)
  ## Get unpenalized (neg log-likelihood) Hessian
  data <- object$tmb_obj$env$data
  data$part <- 1
  newobj <- update(object$tmb_obj, data = data)
  args <- list(...)
  args$object <- newobj
  args$par <- c(object$param$beta, object$param$gamma, object$param$theta)
  args$joint <- FALSE
  args$method <- "analytical"
  args$sparse <- TRUE
  Hl <- do.call("Hess", args)
  dMM <- rowSums(solve(Hp) * Hl)[idx]
  ## Cumulate diagonal elements
  nm <- names(idx)
  nm <- factor(nm, levels = unique(nm))
  idx <- split(seq_along(idx), nm)
  sapply(idx, function(i) sum(dMM[i]))
}
