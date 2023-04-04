##' Response objects
##'
##' Response objects to represent censored and truncated observations
##'
##' @details
##'
##' \code{Resp} extends the functionality of \code{\link[survival]{Surv}} class
##'   by allowing cases that cannot be defined with it. An example is an
##'   interval-censored outcome with left truncatation (see Examples).
##'
##' Censored and exactly observed data can be defined similarly to \code{type =
##'   "interval2"} objects in \code{\link[survival]{Surv}}. \code{NA} values for
##'   left or right censoring borders mean left- or right-censored observations,
##'   respectively. If both borders are \code{NA}, the observation is considered
##'   \code{NA} by \code{is.na()}.  Truncation times (\code{tleft} and
##'   \code{tright} arguments) can be omitted or take \code{NA} values, which
##'   means no truncation. If only the censoring intervals are provided, i.e.,
##'   no trunction is present, the function returns a \code{Surv} object.
##'
##' \code{Resp} also provides a limited interface between \code{tramME} and the
##'   \code{response} class (technically, inherits from it) of \code{mlt} (see
##'   \code{\link[mlt]{R}}), which uses an internal representation that is not
##'   compatible with \code{tramME}.
##'
##' The optional argument \code{open_lwr_bnd} can be used to enforce lower
##'   boundaries of the outcome. Left boundaries in the \code{Resp} object
##'   (\code{cleft} and \code{tleft}) that are equal to the first element of
##'   \code{bounds} will be increased with one \code{tol} value to avoid
##'   downstream numerical problems in \code{mlt}. This adjustment is recorded
##'   and reversed when we print the object.
##'
##' @section Warning:
##'
##' This function is experimental and currently limited to continuous outcome
##'   types. It may be subject to change.
##'
##' @param cleft A vector of left borders of censoring intervals
##' @param cright A vector of right borders of censoring intervals
##' @param tleft A vector of left truncation values
##' @param tright A vector of right truncation values
##' @param bounds An optional numeric vector of two elements (\code{c(a, b)})
##'   that denotes the lower and upper boundaries of the outcome.
##' @param open_lwr_bnd Logical; if \code{TRUE}, the lower boundary of the
##'   outcome is open, and we want to enforce this.
##' @param tol Tolerance level.
##' @return A \code{Resp} object or a \code{Surv} object
##'
##' @examples
##'
##' dat <- data.frame(x1 = 1:10, x2 = c(2:10, NA), x3 = c(NA, 0:8))
##' dat$r <- with(dat, Resp(x1, x2, x3))
##'
##' dat$r
##' dat[1:3, ]$r
##' dat$r[1:3]
##'
##' is.na(dat$r)
##'
##' model.frame(r ~ 1, data = dat, na.action = na.omit)
##'
##' @export
Resp <- function(cleft, cright, tleft, tright,
                 bounds = c(-Inf, Inf), open_lwr_bnd = TRUE,
                 tol = sqrt(.Machine$double.eps)) {
  stopifnot(is.numeric(cleft), is.numeric(cright))
  .lb <- function(x, v) {
    if (length(i <- which((x == v) %in% TRUE))) i
    else NULL
  }
  if (open_lwr_bnd && is.finite(bounds[1])) {
    clwr_bnd <- .lb(cleft, bounds[1])
    cleft[clwr_bnd] <- cleft[clwr_bnd] + tol
    tlwr_bnd <- .lb(tleft, bounds[1])
    tleft[tlwr_bnd] <- tleft[tlwr_bnd] + tol
  } else {
    clwr_bnd <- NULL
    tlwr_bnd <- NULL
  }
  if (missing(tleft) && missing(tright))
    return(survival::Surv(cleft, cright, type = "interval2"))
  stopifnot(is.numeric(tleft))
  out <- cbind(cleft, cright, tleft)
  nms <- c("cleft", "cright", "tleft")
  if (!missing(tright)) {
    stopifnot(is.numeric(tright))
    out <- cbind(out, tright)
    nms <- c(nms, "tright")
  }
  colnames(out) <- nms
  class(out) <- c("Resp", "response", class(out))
  attr(out, "tol") <- tol
  attr(out, "bounds") <- bounds
  attr(out, "open_lwr_bnd") <- open_lwr_bnd
  ## NOTE: these lower bounds are adjusted
  attr(out, "adj_clwr_bnd") <- clwr_bnd
  attr(out, "adj_tlwr_bnd") <- tlwr_bnd
  return(out)
}


##' @param object A \code{Resp} object
##' @param ... Optional arguments
##' @importFrom mlt R
##' @describeIn Resp Converting \code{Resp} objects to \code{response} (from
##'   \code{mlt}) objects (see \code{\link[mlt]{R}})
##' @export
R.Resp <- function(object, ...) {
  bnd <- attr(object, "bounds")
  tol <- attr(object, "tol")
  adj <- if (attr(object, "open_lwr_bnd") && is.finite(bnd[1])) tol
         else 0
  x <- unclass(object)
  ex <- (x[, 1] == x[, 2]) %in% TRUE
  na <- rowSums(is.na(x[, 1:2, drop = FALSE])) == 2
  nc <- ncol(x)
  out <- mlt::R(object = as.numeric(ifelse(ex, x[, 1], NA)),
    cleft  = as.numeric(ifelse(ex, NA,
                        ifelse(is.na(x[, 1]), bnd[1] + adj, x[, 1]))),
    cright = as.numeric(ifelse(ex, NA, ifelse(is.na(x[, 2]), bnd[2], x[, 2]))),
    tleft  = if (nc >= 3) as.numeric(x[, 3]) else NA,
    tright = if (nc == 4) as.numeric(x[, 4]) else NA)
  out$cleft[na] <- out$cright[na] <- NA
  return(out)
}

## Revert adjustments for lower bounds in a Resp object,
## primarily for printing purposes.
.no_adj_Resp <- function(obj) {
  bnd <- attr(obj, "bounds")
  att <- attributes(obj)
  att$open_lwr_bnd <- FALSE ## XXX: to prevent adjusting NAs
  obj <- unclass(obj)
  if (length(idx <- att$adj_clwr_bnd)) obj[idx, "cleft"] <- bnd[1]
  if (length(idx <- att$adj_tlwr_bnd)) obj[idx, "tleft"] <- bnd[1]
  class(obj) <- c("Resp", "response", class(obj))
  anms <- c("tol", "bounds", "open_lwr_bnd", "adj_clwr_bnd",
            "adj_tlwr_bnd")
  attributes(obj)[anms] <- att[anms]
  obj
}

##' @param x A \code{Resp} object
##' @param ... Optional arguments
##' @describeIn Resp Print method for the \code{Resp} class
##' @export
print.Resp <- function(x, ...) {
  invisible(print(mlt::R(.no_adj_Resp(x)), ...))
}

##' @param x A \code{Resp} object
##' @param i Row index (typically the only index)
##' @param j Column index (typically missing)
##' @param drop If \code{TRUE} the result is coerced to the lowest possible dimension
##' @describeIn Resp Subsetting \code{Resp} objects
##' @export
"[.Resp" <- function(x, i, j, drop = FALSE) {
  .adj_idx <- function(i1, i2) {
    if (length(i <- which(i1 %in% i2))) i
    else NULL
  }
  if (missing(j)) {
    idx <- seq(length(x))
    att <- attributes(x)
    x <- unclass(x)[i, , drop = FALSE]
    class(x) <- c("Resp", "response", class(x))
    anms <- c("tol", "bounds", "open_lwr_bnd")
    attributes(x)[anms] <- att[anms]
    ## NOTE: adjust idices according to subsetting
    idx <- idx[i]
    if (length(bi <- att$adj_clwr_bnd))
      attr(x, "adj_clwr_bnd") <- .adj_idx(idx, bi)
    if (length(bi <- att$adj_tlwr_bnd))
      attr(x, "adj_tlwr_bnd") <- .adj_idx(idx, bi)
    return(x)
  } else {## FIXME
    x <- unclass(x)
    NextMethod("[")
  }
}


##' @param x A \code{Resp} object
##' @describeIn Resp Missing values
##' @export
is.na.Resp <- function(x)
  as.vector(rowSums(is.na(unclass(x)[, 1:2, drop = FALSE])) == 2)


##' @param x A \code{Resp} object
##' @describeIn Resp Length of a \code{Resp} object
##' @export
length.Resp <- function(x) nrow(x)

## XXX: this is a hack, because mlt::R does not offer as.character
##' @importFrom utils capture.output
as.character.Resp <- function(x, ...) {
  capt <- capture.output(out <- print.Resp(x, ...))
  out
}

##' @param x A \code{Resp} object
##' @describeIn Resp \code{format} method for a \code{Resp} object
##' @export
format.Resp <- function(x, ...) format(as.character.Resp(x), ...)

