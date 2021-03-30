##' Generic method for extracting baseline transformations
##' @param object A model object
##' @param ... Optional parameters
##' @return The value of the baseline transfromation function at certain points.
##' @export
trafo <- function(object, ...) {
  UseMethod("trafo", object)
}

##' Get the baseline transformation function and its confidence
##' interval
##'
##' For stratified models, it returns a list of data frames for each
##' stratum.
##' @param object A fitted tramME object.
##' @param newdata Values of the interacting terms to be used.
##' @param type The scale on which the transformation function is evaluated.
##' @param confidence Pointwise confidence interval or confidence band.
##' @param level Confidence level.
##' @param K Integer, number of points in the grid the function is
##'   evaluated on.
##' @param ... Additional parameters (for consistency with generic)
##' @return Matrix or list of matrices containing the point estimates and the
##'   confidence intervals.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' tr <- trafo(fit, type = "distribution", confidence = "interval", K = 100)
##' @importFrom stats qnorm model.matrix
##' @importFrom variables mkgrid
##' @export
trafo.tramME <- function(object, newdata = NULL,
                  type = c("trafo", "distribution", "survivor", "cumhazard"),
                  confidence = c("none", "interval", "band"), level = 0.95,
                  K = 50, ...) {
  type <- match.arg(type)
  confidence <- match.arg(confidence)
  rv <- variable.names(object, "response")
  if (is.null(newdata)) {
    newdata <- model.frame(object)
    newdata <- newdata[, variable.names(object, "interacting"), drop = FALSE]
    newdata <- unique(newdata)
  }
  if (nrow(newdata) > 1) {
    out <- lapply(1:nrow(newdata), function(i)
      trafo(object, newdata = newdata[i, , drop = FALSE], type = type,
            confidence = confidence, level = level, K = K))
    names(out) <- apply(newdata, 1, FUN = function(x)  {
      paste(mapply(FUN = function(n, v) paste(n, "=", v),
                   n = names(newdata), v = x), collapse = ", ")})
    class(out) <- c("trafo.tramME", class(out))
    return(out)
  }
  nd <- newdata[rep(1, K), , drop = FALSE]
  nd[[rv]] <- mkgrid(object$model$ctm$bases$response, n = K)[[rv]]
  dist <- switch(object$model$ctm$todistr$name, "normal" = "Normal", "logistic" = "Logistic",
                 "minimum extreme value" = "MinExtrVal",
                 "maximum extreme value" = "MaxExtrVal")
  mod <- ctm(response = object$model$ctm$bases$response,
             interacting = object$model$ctm$bases$interacting,
             todistr = dist)
  X <- model.matrix(mod, data = nd)
  cci <- complete.cases(X)
  X <- X[cci, , drop = FALSE] ## NOTE: drop rows with missing obs
  cf <- coef(object, with_baseline = TRUE)
  cf <- cf[.idx(object, fixed = TRUE, pargroup = "baseline")]
  vc <- vcov(object, pargroup = "baseline")
  if (length(fixn <- setdiff(names(cf), rownames(vc)))) { ## NOTE: handling fixed parameters
    fixi <- match(fixn, names(cf))
    os <- X[, fixi, drop = FALSE] %*% cf[fixi]
    cf <- cf[-fixi]
    X <- X[, -fixi, drop = FALSE]
  } else {
    os <- 0
  }
  out <- switch(confidence,
    none = {
      ci <- X %*% cf
      colnames(ci) <- "trafo"
      ci + os
    },
    interval = {
      ci <- confint(multcomp::glht(multcomp::parm(cf, vc), linfct = X),
                    calpha = multcomp::univariate_calpha(), level = level)$confint
      colnames(ci)[1] <- "trafo"
      ci + os
    },
    band = {
      ci <- confint(multcomp::glht(multcomp::parm(cf, vc), linfct = X),
                    calpha = multcomp::adjusted_calpha(), level = level)$confint
      colnames(ci)[1] <- "trafo"
      ci + os
    })
  dn <- vector("list", 2)
  names(dn) <- c(rv, "estim")
  dn[[rv]] <- nd[cci, rv]
  dn[[2]] <- colnames(out)
  dimnames(out) <- dn
  out <- switch(type,
                trafo = out,
                distribution = mod$todistr$p(out),
                survivor = 1 - mod$todistr$p(out),
                cumhazard = -log(1 - mod$todistr$p(out)))
  attr(out, "scale") <- type
  class(out) <- c("trafo.tramME", class(out))
  return(out)
}



##' Plotting method for \code{trafo.tramME} objects
##'
##' @param x A \code{trafo.tramME} object.
##' @param col Line colors, recycled if shorter than the size of the
##'   \code{trafo.tramME} object.
##' @param fill Fill color for the confidence intervals.
##' @param lty Line types.
##' @param add If \code{TRUE} add to an existing plot.
##' @param ... Additional arguments, passed to \code{plot} or \code{lines}.
##' @return The original \code{trafo.tramME} object, invisibly.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' tr <- trafo(fit, type = "trafo", confidence = "interval", K = 100)
##' plot(tr, col = 2, main = "Trafo")
##' @importFrom graphics plot lines polygon
##' @export
plot.trafo.tramME <- function(x, col = 1, fill = "lightgrey", lty = 1,
                              add = FALSE, ...) {
  ## --- default values
  fcall <- match.call(expand.dots = TRUE)
  if (inherits(x, "list")) {
    if (is.null(fcall$xlim))
      fcall$xlim <- range(sapply(x, function(x) as.numeric(rownames(x))),
                          na.rm = TRUE)
    if (is.null(fcall$ylim))
      fcall$ylim <- range(unlist(x), na.rm = TRUE)
    if (is.null(fcall$xlab))
      fcall$xlab <- names(dimnames(x[[1]]))[1]
    if (is.null(fcall$ylab))
      fcall$ylab <- paste("baseline", attr(x[[1]], "scale"), "function")
  } else {
    xx <- as.numeric(rownames(x))
    if (is.null(fcall$xlim))
      fcall$xlim <- range(xx, na.rm = TRUE)
    if (is.null(fcall$ylim))
      fcall$ylim <- range(x, na.rm = TRUE)
    if (is.null(fcall$xlab))
      fcall$xlab <- names(dimnames(x))[1]
    if (is.null(fcall$ylab))
      fcall$ylab <- paste("baseline", attr(x, "scale"), "function")
  }
  ## --- New plot if not added
  if (!add) {
    m <- match(c("xlim", "ylim", "xlab", "ylab", "main"), names(fcall), 0L)
    fc <- fcall[c(1L, m)]
    fc$type <- "n"
    fc$x <- 0
    fc[[1L]] <- quote(plot)
    out <- eval(fc)
  }
  ## ---
  if (inherits(x, "list")) {
    out <- mapply(FUN = function(p, c, f, t) {
      plot(p, col = c, fill = f, lty = t, add = TRUE, ...)},
      p = x, c = col, f = fill, t = lty)
  }
  if (inherits(x, "matrix")) {
    if (all(c("lwr", "upr") %in% colnames(x))) {
      lwr <- x[, "lwr"]
      upr <- x[, "upr"]
      polygon(c(xx, rev(xx)), c(lwr, rev(upr)), border = NA, col = fill)
    }
    lines(xx, x[, "trafo"], col = col, lty = lty, ...)
  }
  invisible(x)
}
