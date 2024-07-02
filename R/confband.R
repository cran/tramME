## Confidence bands

##' Confidence intervals and bands from a \code{tramME} model
##'
##' Pointwise confidence intervals or multiplicity-adjusted confidence bands for
##'     transformation, distribution, survivor or cumulative hazard functions.
##'
##' @details
##'
##'   Similarly to \code{\link[mlt]{confband}}, this method evaluates the
##'   conditional distribution of the outcome on a selected scale given a number
##'   of grid-points and claculates the corresponding confidence intervals or
##'   bands (adjusting for multiplicity).
##'
##'   The point estimates retured by this function could also be calculated with
##'   \code{\link{predict.tramME}} (when \code{newdata} does not contain
##'   response values and \code{K} is set to the number of grid points).
##'   While \code{\link{predict.tramME}} is designed to calculate a
##'   potentially large number of point estimates on a wider range of available
##'   scales, \code{confband} calculates the asymptotic intervals from the joint
##'   covariance matrix of the fixed and random effects. For technical reasons,
##'   a smaller set of \code{type} options are available, and the calculations
##'   are slower than with \code{\link{predict.tramME}}. The handling of random
##'   effects is currently stricter than in \code{predict.tramME}: No
##'   \code{ranef} option is available, and grouping factors for random effects
##'   supplied in \code{newdata} must have the same levels as the dataset used
##'   to fit the model.
##'
##'   The multiplicity adjustment is done using
##'   \code{\link[multcomp]{confint.glht}}. The \code{cheat} argument reduces
##'   the dimensionality of the multivariate root-finding problem (see
##'   \code{\link[mvtnorm]{qmvt}}) for speed and (occasionally) numerical
##'   stability. The critical value for the confidence bands are obtained for
##'   \code{cheat < K} grid points, but the confidence bands are calculated for
##'   \code{K} grid points. As a result, the nominal level of the returned
##'   confidence band is not maintained, but the deviation is expected to be
##'   small if \code{cheat} is reasonably large. It is the user's responsibility
##'   to set this value, and by default \code{cheat = K}.
##'
##' @section Warning:
##'
##'   This method implements new functionality. Its user interface may
##'   be subject to change.
##'
##' @param object The \code{tramME} object.
##' @param type The scale for which the confidence bands are calculated.
##' @param newdata A data frame of covariate values.
##' @param level Confidence level.
##' @param type The scale on which the condfidence bands are calculated.
##' @param adjust If \code{TRUE}, multiplicity-adjusted confidence bands are
##'   calculated. (see Details)
##' @param K The number of grid points at which the outcome distribution is
##'   evaluated.
##' @param cheat In the case of multiplicity adjustment (\code{adjust = TRUE}),
##'   an option to decrease the number of grid points (\code{cheat < K}), for
##'   faster calculations and increased numerical stability. (see Details)
##' @param q The quantiles at which the model is evaluated.
##' @param baseline_only If \code{TRUE}, only evaluate the baseline
##'   transformation function and ignore the shift terms.
##' @param ... Optional arguments passed to \code{\link[multcomp]{confint.glht}}.
##' @return A matrix (in the case when \code{newdata} has a single row) or a
##'   list of matrices for each row of \code{newdata}.
## @exportS3Method mlt::confband
##' @importFrom mlt confband
##' @importFrom variables mkgrid
##' @export
## TODO: add q argument to control at which quantiles, figure out how to work with cheat
## TODO: other scales
confband.tramME <- function(object, newdata, level = 0.95,
                            type = c("trafo", "distribution", "survivor", "cumhazard"),
                            adjust = FALSE,
                            K = 40, cheat = K, q = NULL,
                            baseline_only = FALSE, ...) {
  stopifnot(!missing(newdata))
  stopifnot(is.data.frame(newdata) || is.null(newdata))
  type <- match.arg(type)
  if (NROW(newdata) > 1) {
    newdata <- split(newdata, seq(nrow(newdata)))
    out <- lapply(newdata, function(nd) {
      confband(object, nd, level, type, adjust,
               K, cheat, q, baseline_only, ...)
    })
    out <- unname(out)
    class(out) <- c("confband.tramME", class(out))
    return(out)
  }
  ## NOTE: stricter than predict and checks whether levels of random effects in
  ## newdata match the ones in data to rule out confusion.
  if (length(ren <- variable.names(object, "ranef"))) {
    stopifnot(sapply(ren, function(n) {
      identical(levels(model.frame(object)[[n]]),
                levels(newdata[[n]]))
    }))
  }
  if (is.null(newdata)) {
    if (!is.null(object$model$ctm$bases$interacting))
      stop(paste0("newdata must be specified when the ",
                  "model contains interacting terms"))
    baseline_only <- TRUE
    newdata <- model.frame(object)[1L, ]
  }
  vn <- variable.names(object, "response")
  if (is.null(q))
    y <- mkgrid(object$model$ctm, n = K)[[vn]]
  else {
    y <- q
    K <- length(y)
  }
  ## -- NOTE: manually adjust lower bound to avoid log(0) when log_first = TRUE
  if (attr(object$model$ctm$bases$response, "log_first"))
    y[1L] <- max(y[1L], sqrt(.Machine$double.eps))
  ## --
  nd <- newdata[rep(1, K), , drop = FALSE]
  nd[[vn]] <- y
  mm <- model.matrix(object, data = nd,
                     drop_unused_groups = FALSE)
  pr <- coef(object, complete = TRUE)
  if (isTRUE(baseline_only)) X <- mm$Ye
  else X <- cbind(mm$Ye, mm$X, t(as.matrix(mm$Zt)))
  pr <- pr[colnames(X)]
  vc <- vcov(object, parm = colnames(X))
  if (isTRUE(adjust)) calpha <- multcomp::adjusted_calpha()
  else calpha <- multcomp::univariate_calpha()
  if (isTRUE(adjust) && (cheat < (nr <- nrow(X)))) { ## -- adjust with cheat
    X_ <- X[seq(1, nr, length.out = cheat), ]
    ci <- confint(
      multcomp::glht(multcomp::parm(pr, vc), linfct = X_),
      level = level, calpha = calpha, ...)$confint
    calpha <- attr(ci, "calpha")
  }
  ci <- confint(multcomp::glht(multcomp::parm(pr, vc), linfct = X),
                level = level, calpha = calpha, ...)$confint
  Fz <- object$model$ctm$todistr$p
  if (type == "distribution") ci <- Fz(ci)
  if (type == "survivor") ci <- (1 - Fz(ci))[, c(1, 3, 2)]
  if (type == "cumhazard") ci <- -log(1 - Fz(ci))
  out <- cbind(y, ci)
  colnames(out) <- c(vn, "est", "lwr", "upr")
  rownames(out) <- NULL
  attr(out, "cond") <- newdata
  attr(out, "scale") <- type
  class(out) <- c("confband.tramME", class(out))
  return(out)
}

##' @param x The object containing the confidence intervals.
##' @param col Color of the point estimates.
##' @param lty Line type of the point estimates.
##' @param fill Fill color for the intervals.
##' @param add If \code{TRUE}, no new plot is created, the interval is added to
##'   the current plot.
##' @param single_plot If \code{TRUE}, a single new plot is created, and all
##'   intervals are plotted on it.
##' @param trafo_x Transform x-axis before plotting.
##' @param trafo_y Transform y-axis before plotting.
##' @param align_xlim If \code{TRUE}, align the x-axis limits across all
##'   subplots.
##' @param align_ylim If \code{TRUE}, align the y-axis limits across all
##'   subplots.
##' @param ... Optional arguments passed to \code{\link[graphics]{plot.default}}
##'   and \code{\link[graphics]{plot.xy}}.
##' @importFrom graphics plot.default plot.xy par polygon lines
##' @importFrom grDevices grey col2rgb n2mfrow
##' @name plot_ci
NULL

## Helper function to plot confidence intervals
## NOTE TODO A slightly modified version of this function is exported in BTmisc
## package
## TODO: titles
plot_CIs <- function(x,
                     col, lty, fill,
                     add = FALSE, single_plot = FALSE,
                     trafo_x = identity, trafo_y = identity,
                     align_xlim = FALSE, align_ylim = FALSE,
                     ...,
                     xlab_attr = NULL, ylab_attr = NULL) {
  if ((n <- length(x)) == 0) return(invisible(x))
  stopifnot(is.list(x))
  plot_args <- list()
  if (align_xlim || single_plot) {
    plot_args$xlim <- range(sapply(x, `[`, , 1), na.rm = TRUE)
  }
  if (align_ylim || single_plot) {
    plot_args$ylim <- range(sapply(x, `[`, , 2:4), na.rm = TRUE)
  }
  opt_args <- as.list(match.call(expand.dots = FALSE))$`...`
  ## Get optional arguments for plot.default
  anm <- intersect(names(opt_args), formalArgs(plot.default))
  plot_args[anm] <- opt_args[anm]
  anm <- intersect(names(opt_args), formalArgs(plot.xy))
  lines_args <- opt_args[anm]
  if (add || single_plot) defs <- seq_along(x)
  else defs <- 1
  if (missing(col)) col <- defs
  if (missing(lty)) lty <- defs
  if (missing(fill)) {
    fill <- sapply(defs, function(x) {
      do.call("rgb", c(as.list(col2rgb(x)[, 1]),
        max = 255, alpha = 0.2 * 255))
    })
  }
  if (is.null(fill)) fill <- NA
  col <- rep(col, length.out = n)
  lty <- rep(lty, length.out = n)
  fill <- rep(fill, length.out = n)
  if (!single_plot && !add && n > 1 && par("page")) {
    pp <- par(mfrow = rev(n2mfrow(length(x))))
    on.exit(par(pp))
  }
  for (i in seq_along(x)) {
    xx <- trafo_x(x[[i]][, 1L])
    yy <- trafo_y(x[[i]][, 2L])
    yl <- trafo_y(x[[i]][, 3L])
    yu <- trafo_y(x[[i]][, 4L])
    if (!add) {
      if (!single_plot || (single_plot && i == 1L)) {
        xlab <- if (!is.null(xlab_attr)) attr(x[[i]], xlab_attr) else NULL
        ylab <- if (!is.null(ylab_attr)) attr(x[[i]], ylab_attr) else NULL
        if (is.null(xlab)) xlab <- colnames(x[[i]])[1L]
        if (is.null(ylab)) ylab <- names(x)[1L]
        if (is.null(xlab)) xlab <- "x"
        if (is.null(ylab)) ylab <- "y"
        if (!((tx <- deparse1(substitute(trafo_x))) %in% c("identity", "I")))
          xlab <- paste(c(tx, "(", xlab, ")"), collapse = "")
        if (!((ty <- deparse1(substitute(trafo_y))) %in% c("identity", "I")))
          ylab <- paste(c(ty, "(", ylab, ")"), collapse = "")
        arg <- list(xlim = range(xx, na.rm = TRUE),
                    ylim = range(yy, yl, yu, na.rm = TRUE),
                    xlab = xlab, ylab = ylab)
        arg[names(plot_args)] <- plot_args
        arg$type <- "n"
        arg$x <- 0
        do.call("plot.default", arg)
      }
    }
    arg <- list(x = c(xx, rev(xx)),
                y = c(yl, rev(yu)),
                border = NA, col = fill[i])
    do.call("polygon", arg)
    arg <- lines_args
    arg$col <- col[i]
    arg$lty <- lty[i]
    do.call("lines.default", c(list(x = xx, y = yy), arg))
  }
  invisible(x)
}

##' Plot confidence bands from \code{tramME} models
##'
##' Plotting method for \code{confband.tramME} objects.
##' @rdname plot_ci
##' @export
plot.confband.tramME <- function(x, col, lty, fill,
                                 add = FALSE, single_plot = FALSE,
                                 trafo_x = identity, trafo_y = identity,
                                 align_xlim = FALSE, align_ylim = FALSE,
                                 ...) {
  if (!is.list(x)) x <- list(x)
  fc <- match.call(expand.dots = TRUE)
  fc$x <- x
  fc$ylab_attr <- "scale"
  fc[[1L]] <- quote(plot_CIs)
  eval(fc)
}
