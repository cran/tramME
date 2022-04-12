## Functions for generating the various terms of a tramME model
## FIXME: decide is these functions should create initial parameter vectors


## Create data and other required values for the parametric shift terms
##
## @param \code{ctm} model describing the baseline transformation and parametric shifts
## @param data \code{data.frame} containing the variables of the model
## @param ... Optional arguments, passed to \code{model.matrix}.
## @return A matrix containing the data for the parametric shift terms. Its attributes
## contain the parameter constraints.
## @note The \code{ctm} already contains information on the sign of the shift terms;
## no \code{negative} argument is necessary.
.shift_trms <- function(ctm, data, ...) {
  bs <- ctm$model$bshifting
  if (is.null(bs)) {
    out <- matrix(0, nrow = nrow(data), ncol = 0)
    attr(out, "constraint") <- list(ui = Matrix::Diagonal(0),
                                    ci = rep(-Inf, 0))
  } else {
    out <- model.matrix(bs, data = data, ...)
  }
  attr(out, "parnames") <- if (length(nm <- colnames(out))) nm else character(0)
  return(out)
}


## Create data and other required values for the smooth terms
##
## @details
## If \code{smooth} is of \code{tramME_smooth} class, it has an additional \code{data}
## attribute that contains the observations to be used to set up the smooth terms.
## @param smooth a list of smooth terms from \code{\link[mgcv]{interpret.gam}} or
##   an object of \code{tramME_smooth} class.
## @param data \code{data.frame} containing the variables of the model
## @param negative logical value that indicates whether shift terms of the \code{tramME}
## model are parameterized with negative signs
## @return A list containing data and parameter values to be used in the TMB model.
## @examples
## ## When the model is already set up and we want to use the originally set-up
## ## smooth terms and evaulate on a new dataset
## sm <- .smooth_trms(.tramME_smooth(object), data, negative = FALSE)
.smooth_trms <- function(smooth, data, negative) {
  unms <- function(nms, pf) {
    if (is.null(nms)) return(character(0))
    unlist(lapply(split(nms, factor(nms, levels = unique(nms))),
                  function(n) { paste0(n, pf, seq_along(n)) }))
  }
  out <- list()
  setup_data <- if (inherits(smooth, "tramME_smooth")) attr(smooth, "data")
                else NULL
  if (length(smooth)) {
    if (length(setup_data))
      bs <- lapply(smooth, mgcv::smoothCon, data = setup_data, absorb.cons = TRUE)
    else
      bs <- lapply(smooth, mgcv::smoothCon, data = data, absorb.cons = TRUE)
    Xs <- list()
    Zs <- list()
    re_nms <- character(0)
    bs_nms <- character(0)
    for (i in 1:length(bs)) {
      for (j in 1:length(bs[[i]])) {
        sm <- bs[[i]][[j]]
        if (is.null(setup_data)) { ## same data for X,Y as setting up the smooth
          rp <- mgcv::smooth2random(sm, names(data), 2)
        } else { ## for new data
          rp <- .s2rPred(sm, data)
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
    X <- do.call("cbind", Xs)
    if (length(Zs)) {
      Z <- as(do.call("cbind", Zs), "dgTMatrix")
      attr(Z, "ranef") <- list(termsize = sapply(Zs, ncol))
    } else {
      Z <- as(matrix(0, nrow = nrow(X), ncol = 0), "dgTMatrix")
      attr(Z, "ranef") <- list(termsize = integer(0))
    }
  } else {
    X <- matrix(0, nrow = nrow(data), ncol = 0)
    Z <- as(matrix(0, nrow = nrow(X), ncol = 0), "dgTMatrix")
    attr(Z, "ranef") <- list(termsize = integer(0))
  }
  if (negative) {
    X <- -X
    Z <- -Z
  }
  ts <- attr(Z, "ranef")$termsize
  attr(X, "constraint") <- list(
    ui = Matrix::Diagonal(ncol(X)),
    ci = rep(-Inf, ncol(X))
  )
  attr(X, "parnames") <- as.vector(unms(colnames(X), "|FE"))
  attr(Z, "constraint") <- list(
    ui = Matrix::Diagonal(length(ts)),
    ci = rep(-Inf, length(ts))
  )
  attr(Z, "parnames") <- as.vector(unms(colnames(Z), "|RE"))
  return(list(X = X, Z = Z))
}


## Create random effects data and other required values
##
## @param ranef a list of random effects formulas from \code{\link[lme4]{findbars}}
## @param data data.frame containing the variables of the model
## @param negative logical value that indicates whether the random effects have
##   a negative sign
## @param ... Passed to \code{\link[lme4]{mkReTrms}}.
## @return A list containing data and parameter values to be used in the TMB model.
.ranef_trms <- function(ranef, data, negative, ...) {
  ## XXX: a safer version of lme4::mkReTrms
  mkReTrms <- function(bars, fr, ...) {
    fc <- match.call()
    m <- match(c("bars", "fr", "drop.unused.levels", "reorder.terms",
                 "reorder.vars"), names(fc), 0L)
    fc <- fc[c(1L, m)]
    fc[[1L]] <- quote(lme4::mkReTrms)
    eval(fc, parent.frame())
  }

  if (length(ranef)) {
    rt <- mkReTrms(ranef, data, ...)
    Zt <- rt$Zt
    ri <- list(
      termsize = sapply(rt$Ztlist, NROW),
      blocksize = sapply(rt$cnms, length),
      names = rt$cnms,
      levels = lapply(rt$flist, levels)
    )
    ri$npar <- sum(ri$blocksize * (ri$blocksize + 1) / 2)
    pn <- unlist(mapply(function(n, g) {
      gnl <- expand.grid(g, "|", n, ":", ri$levels[[g]],
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      apply(gnl, 1, paste, collapse = "")
    }, n = ri$names, g = names(ri$names),
    SIMPLIFY = FALSE, USE.NAMES = FALSE))
  } else {
    Zt <- as(matrix(0, nrow = 0, ncol = nrow(data)), "dgTMatrix")
    ri <- list(termsize = integer(0), blocksize = integer(0),
               npar = 0L)
    pn <- character(0)
  }
  if (negative) Zt <- -Zt
  attr(Zt, "ranef") <- ri
  attr(Zt, "constraint") <- list(
    ui = Matrix::Diagonal(ri$npar),
    ci = rep(-Inf, ri$npar)
  )
  attr(Zt, "parnames") <- pn
  return(Zt)
}

## Create data and other required values for the baseline and interacting terms
##
## @param ctm model describing the baseline transformation and parametric shifts
## @param data \code{data.frame} containing the variables of the model
## @param ... Optional arguments, passed to \code{model.matrix}.
## @return A list of matrices containing data for exact, censored (left, right,
## interval) and truncated (left, right, interval) observations calculated by
## evaluating the corresponding basis functions.
## Indices of the corresponding observations in the original data set are given
## in the \code{which} attribute of these matrices.
## The \code{constraint} attribute of the output list contains the necessary
## parameter constraints.
## FIXME: add beta vector?
.base_trms <- function(ctm, data, ...) {
  vn <- variable.names(ctm, "response")
  rv <- as.data.frame(mlt::R(data[[vn]]))
  ## -- XXX: In some cases the truncation bounds are outside the bounds
  ## of the response. Ex: 0 when log_first = TRUE. This is handled in
  ## mlt pragmatically, and the rows with NaNs are removed when the ll
  ## is calculated. Here we run a check on the truncation bounds and
  ## modify them when necessary.
  bs <- attr(ctm$bases$response, "variables")$bounds
  rv$tleft[rv$tleft < bs[1]] <- NA
  rv$tright[rv$tright > bs[2]] <- NA
  ## --
  bs <- if (is.null(ctm$model$bresponse)) ctm$model$binteracting
        else ctm$model$bresponse
  cases <- list(e = quote(is.finite(exact)),
                r = quote(is.finite(cleft) & is.infinite(cright)),
                l = quote(is.finite(cright) & is.infinite(cleft)),
                i = quote(is.finite(cleft) & is.finite(cright)),
                tl = quote(is.finite(tleft) & is.na(tright)),
                tr = quote(is.finite(tright) & is.na(tleft)),
                ti = quote(is.finite(tleft) & is.finite(tright)))
  ## XXX: in tl & tr, is.nas work for Inf & -Inf because of the previous
  ## bounds adjustment
  names <- list(e = "e", r = "r", l = "l", i = c("il", "ir"), tl = "tl",
                tr = "tr", ti = c("til", "tir"))
  ## --
  out <- list()
  cons <- NULL
  nc <- NULL
  for (i in seq_along(cases)) {
    vars <- all.vars(cases[[i]], functions = TRUE, unique = FALSE)
    vars <- vars[which(vars == "is.finite") + 1]
    idx <- eval(as.call(cases[[i]]), rv)
    if (sum(idx)) {
      dat <- data[idx, , drop = FALSE]
      for (j in seq_along(vars)) {
        nm <- paste0("Y", names[[i]][j])
        dat[[vn]] <- rv[idx, vars[j]]
        Y <- model.matrix(bs, dat, ...)
        if (is.null(cons)) cons <- attr(Y, "constraint")
        if (is.null(nc)) nc <- ncol(Y)
        attr(Y, "constraint") <- attr(Y, "Assign") <- NULL
        attr(Y, "which") <- which(idx)
        out[[nm]] <- Y
      }
      if (names(cases[i]) == "e") {
        deriv <- 1L
        names(deriv) <- vn
        Y <- model.matrix(bs, dat, deriv = deriv)
        attr(Y, "constraint") <- attr(Y, "Assign") <- NULL
        attr(Y, "which") <- which(idx)
        out$Yeprime <- Y
      }
    } else {
      for (j in seq_along(vars)) {
        nm <- paste0("Y", names[[i]][j])
        out[[nm]] <- NA
      }
      if (names(cases[i]) == "e") {
        out$Yeprime <- NA
      }
    }
  }
  empty <- !sapply(out, is.matrix)
  Y <- out[[which(!empty)[1L]]][0, ] ## an empty matrix with correct structure
  attr(Y, "which") <- integer(0)
  out[which(empty)] <- rep.int(list(Y), sum(empty))
  attr(out, "constraint") <- cons
  return(out)
}


##' Extract model frame from a tramME model
##'
##' @details
##' In \code{\link[mlt]{mlt}}, the basis functions expect the response variables
##' in the data to be evaluated, i.e. instead of \code{x} and \code{y} columns
##' we should have a \code{`Surv(x, y)`} column when the response is a
##' \code{\link[survival]{Surv}} object. \code{model.frame.tramME} builds the
##' model frame accordingly, assigning to the resulting object the class
##' \code{tramME_data} to indicate this structure to other functions that use
##' its results. If the input \code{data} is a \code{tramME_data} is also expects
##' this structure.
##' @param formula A \code{tramME} object.
##' @param group_as_factor Logical; If \code{TRUE}, automatically convert the
##'   grouping variables of the random effects to factors. (not used, might not be needed) ## FIXME
##' @param ignore_response Logical; If \code{TRUE}, the response is not added to the
##'   result. In this case the function won't look for it in \code{data}.
##' @param ... Optional arguments, passed to \code{\link[stats]{model.frame}}.
##' @return A \code{tramME_data} object, which is also a \code{data.frame}.
##' @importFrom stats model.frame
##' @inheritParams stats::model.frame
##' @examples
##' data("sleepstudy", package = "lme4")
##' mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
##' model.frame(mod)
##' @export
model.frame.tramME <- function(formula, data = NULL, group_as_factor = FALSE,
                               ignore_response = FALSE, ...) {
  dots <- match.call(expand.dots = FALSE)$...
  fc <- formula$call
  if (nodata <- missing(data)) data <- formula$data
  ff <- fake_formula(formula, drop_response = ignore_response)
  ## -- NOTE: In the case of tramME_data, the lhs of the formula is already
  ## evaluated: in the case of Surv(y, c) ~, we have a column `Surv(y, c)`
  ## instead of x and c
  if (inherits(data, "tramME_data")) {
    if (nodata) {
      ## NOTE: using the data saved in the object
      ff <- fc <- NULL
      trms <- attr(data, "terms")
    }
    if (length(ff) == 3L) ff[[2L]] <- as.name(deparse1(ff[[2L]]))
    ## XXX: terms attached to the original data can cause mix-ups
    ## safer to remove them and add at the end
    attr(data, "terms") <- NULL
  }
  if (is.null(fc)) {
    fc <- match.call()
    fc$group_as_factor <- fc$ignore_response <- NULL
  }
  fc$formula <- ff
  fc$data <- data
  m <- match(c("formula", "data", "subset", "weights",
               "na.action", "offset"), names(fc), 0L)
  fc <- fc[c(1L, m)]
  fc[names(dots)] <- dots
  env <- environment(formula$model$formula)
  if (is.null(env)) env <- parent.frame()
  fc[[1L]] <- quote(stats::model.frame)
  mf <- eval(fc, env)
  if (group_as_factor && length(formula$model$ranef)) {
    gv <- variable.names(formula, "grouping")
    mf[gv] <- lapply(mf[gv], as.factor)
  }
  ## NOTE: add back terms when reusing the data saved in the object
  if (nodata && inherits(data, "tramME_data")) attr(mf, "terms") <- trms
  ## -- NOTE: to recognize later the special formatting
  class(mf) <- c("tramME_data", class(mf))
  ## --
  return(mf)
}


##' Model matrices for \code{tramME} models
##'
##' Model matrix for fixed effects, random effects, and baseline transformations
##' (with interacting terms if present).
##'
##' @details
##' Creates model matrices for fixed effects (\code{type = "X"}) and random
##' effects (\code{type = "Zt"}) and baseline transfromation (\code{type = "Y"}),
##' by evaluating the respective basis functions given a new dataset.
##'
##' The response values may be exact, censored (left, right, interval) and
##' truncated (left, right, interval), and the function returns several,
##' potentially empty, model matrices:
##'   \itemize{
##'     \item Ye: Exact observations.
##'     \item Yeprime: The model matrix corresponding to the first derivative
##'           of the baseline transformation, evaluated at exact observations.
##'     \item Yl: Left-censored observations.
##'     \item Yr: Rigt-censored observations.
##'     \item Yil and Yir: Interval-censored observations evaluated at the
##'           left and right bounds of the interval.
##'     \item Ytl: Left-truncated observations.
##'     \item Ytr: Rigt-truncated observations.
##'     \item Ytil and Ytir: Interval-truncated observations evaluated at the
##'           left and right bounds of the interval.
##'   }
##' for the baseline transfromations (unless \code{simplify = TRUE}).
##'
##' @param object A \code{tramME} object.
##' @param data A \code{data.frame} containing the variable values.
##' @param type "X": Fixed effects model matrix. "Zt": Random effects model matrix
##'   (transposed). "Y": Model matrices for the baseline transfromations.
##' @param drop_unused_groups Logical; remove unused levels of the random effects.
##'   (see \code{drop.unused.levels} argument of \code{\link[lme4]{mkReTrms}})
##' @param keep_sign Logical; the terms will have the same sign as in the
##'   \code{tramME} model if \code{TRUE}.
##' @param simplify Logical; Remove empty \code{Y} matrices.
##' @param ... Optional arguments.
##' @return List of requested model matrices.
##' @note The model matrix of the random effects is a sparse matrix and it is transposed
##'   to be directly used with \code{Matrix::crossprod} which is faster than transposing
##'   and multiplying ("Zt" instead of "Z").
##' @importFrom stats model.matrix
##' @export
##' @examples
##' library("survival")
##' rats$litter <- factor(rats$litter)
##' m <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats,
##'              log_first = TRUE, nofit = TRUE)
##' mm <- model.matrix(m)
##' nd <- model.frame(m)[rep(1, 100), ]
##' nd[[1]] <- seq(1, 120, length.out = 100)
##' mm2 <- model.matrix(m, data = nd, simplify = TRUE)
##' mm3 <- model.matrix(m, data = nd, simplify = TRUE, drop_unused_groups = TRUE)
##' ## compare mm2$Zt & mm3$Zt
model.matrix.tramME <- function(object, data = model.frame(object),
                                type = c("Y", "X", "Zt"),
                                drop_unused_groups = FALSE,
                                keep_sign = TRUE, simplify = FALSE,
                                ...) {
  stopifnot(all(type %in% c("Y", "X", "Zt")))
  if (!inherits(data, "tramME_data")) {
    fc <- match.call()
    fc$formula <- object
    fc$object <- fc$type <- fc$drop_unused_groups <- fc$keep_sign <-
      fc$simplify <- NULL
    fc[[1L]] <- quote(model.frame)
    data <- eval(fc, parent.frame())
  }
  neg <- if (keep_sign) object$model$negative else FALSE
  ## -- NOTE: we use different data for setting up the spline bases and
  ## for evaluating them
  sm <- .smooth_trms(.tramME_smooth(object), data = data, negative = neg)
  ## --
  out <- list()
  if ("Zt" %in% type) {
    Zt <- .ranef_trms(object$model$ranef, data = data, negative = neg,
                      drop.unused.levels = drop_unused_groups, ...)
    out$Zt <- rbind(Zt, Matrix::t(sm$Z))
    attr(out$Zt, "parnames") <- c(attr(Zt, "parnames"), attr(sm$Z, "parnames"))
    attr(out$Zt, "type") <- c(rep("re", nrow(Zt)), rep("sm", ncol(sm$Z)))
    ## TODO: select other attributes
  }
  if ("X" %in% type) {
    X <- .shift_trms(object$model$ctm, data = data, ...)
    if (!keep_sign && object$model$negative) X <- -X
    out$X <- cbind(X, sm$X)
    attr(out$X, "parnames") <- c(attr(X, "parnames"), attr(sm$X, "parnames"))
    attr(out$X, "type") <- c(rep("sh", ncol(X)), rep("sm", ncol(sm$X)))
  }
  if ("Y" %in% type) {
    Y <- .base_trms(object$model$ctm, data, ...)
    if (simplify) Y <- Y[sapply(Y, nrow) > 0]
    out[names(Y)] <- Y
  }
  return(out)
}

## A helper function to construct a fake formula from the model elements
## for \code{\link[stats]{model.frame.default}} that contains all
## necessary information.
## @param obj A \code{tramME} object.
## @param drop_response If \code{TRUE} only return the right side of the formula.
## @return A formula, with an additional class \code{fake_formula}.
## TODO: keep offset and weights
fake_formula <- function(obj, drop_response = FALSE) {
  vn <- variable.names(obj, which = "all")
  rv <- if (drop_response) NULL else vn[1L]
  lv <- if (length(vn) == 1L) "1" else vn[-1L]
  ff <- reformulate(lv, response = rv,
                    env = environment(obj$model$formula))
  class(ff) <- c("fake_formula", class(ff))
  return(ff)
}

## FIXME:
## tramTMB_inputs use model.matrix.tramME(list(model = model))
