##' Generic method for extracting terms of the linear predictor
##' @param object A model object
##' @param ... Optional parameters
##' @return The value of the baseline transfromation function at certain points.
##' @export
lpterms <- function(object, ...) {
  UseMethod("lpterms", object)
}


##' Get inidvidual terms of the linear predictor and their confidence
##' intervals
##'
##' \code{term} can be variable groups (baseline/interacting or shift) or names of
##' variables.
##' @note Currently it only takes the fixed effects into account when calculating
##'   intervals (either pointwise confidence intervals or confidence bands).
##' @param object A \code{tramME} object.
##' @param term The names or identifiers of the terms we want to evaluate.
##' @param newdata A \code{data.frame} containing the values at which the
##'   functions are evaluated.
##' @param type The scale on which the functions are evaluated.
##' @param confidence Pointwise confidence interval or confidence band.
##' @param level Confidence level.
##' @param K Integer, number of points of the grid the function is
##'   evaluated on.
##' @param ... Additional parameters (for consistency with generic)
##' @return Matrix or list of matrices containing the point estimates and the
##'   confidence intervals.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' tr <- lpterms(fit, type = "distribution", confidence = "interval", K = 100)
##' @importFrom stats model.matrix complete.cases
##' @importFrom variables mkgrid
##' @export
## FIXME: 'prediction interval' that takes random effects into account.
## FIXME: add 'cheat' argument
## FIXME: implement other term options
lpterms.tramME <- function(object, newdata = model.frame(object)[, -1L],
                           term = c("baseline", "shift"),
                           type = c("trafo", "distribution", "survivor", "cumhazard"),
                           confidence = c("none", "interval", "band"),
                           level = 0.95, K = 50, ...) {
  type <- match.arg(type)
  confidence <- match.arg(confidence)
  term <- term[1]
  rv <- variable.names(object, "response")

  if (is.factor(rvar <- model.frame(object)[[rv]]))
    K <- nlevels(rvar)

  if (!(rv %in% names(newdata))) {

    if (nrow(newdata) > 1) {
      out <- lapply(1:nrow(newdata), function(i)
        lpterms(object, newdata = newdata[i, , drop = FALSE], term = term,
                type = type, confidence = confidence, level = level, K = K))
      names(out) <- apply(newdata, 1, FUN = function(x)  {
        paste(mapply(FUN = function(n, v) paste(n, "=", v),
                     n = names(newdata), v = x), collapse = ", ")})
      class(out) <- c("terms.tramME", class(out))
      return(out)
    }

    nd <- newdata[rep(1, K), , drop = FALSE]
    nd[[rv]] <- mkgrid(object$model$ctm$bases$response, n = K)[[rv]]
  } else {
    nd <- newdata
  }

  ## -- identify the terms
  if (term == "baseline") {
    pns <- .idx(object, fixed = TRUE, pargroup = term)
    if (length(variable.names(object, "interacting"))) {
      tms <- unique(sub("^.*:", "", names(pns)))
      pns <- lapply(tms, function(x) pns[grepl(x, names(pns), fixed = TRUE)])
      names(pns) <- tms
    } else {
      pns <- list(baseline = pns)
    }
  } else if (term == "shift") {
    stop("Shift option is not implemented yet.")
  } else {
    pns <- .idx(object, fixed = TRUE, pargroup = "fixef")
    pns <- lapply(term, function(x) pns[grepl(x, names(pns))])
    names(pns) <- term
  }

  if (length(pns) == 0) {
    return(NULL)
  }

  X <- model.matrix(object$model$ctm, data = nd)
  cf <- coef(object, with_baseline = TRUE)
  vc <- vcov(object)
  out <- lapply(pns, function(pn) {
    X_ <- X[, names(pn), drop = FALSE]
    cci <- complete.cases(X_)
    X_ <- X_[cci, , drop = FALSE] ## NOTE: drop rows with missing obs
    cf_ <- cf[names(pn)]
    vc_ <- vc[names(pn), names(pn), drop = FALSE]
    if (length(fixn <- setdiff(names(cf_), rownames(vc_)))) { ## NOTE: handling fixed parameters
      fixi <- match(fixn, names(cf_))
      os <- X_[, fixi, drop = FALSE] %*% cf_[fixi]
      cf_ <- cf_[-fixi]
      X_ <- X_[, -fixi, drop = FALSE]
    } else {
      os <- 0
    }
    res <- switch(confidence,
      none = {
        ci <- X_ %*% cf_
        colnames(ci) <- "est"
        ci + os
      },
      interval = {
        ci <- confint(multcomp::glht(multcomp::parm(cf_, vc_), linfct = X_),
                      calpha = multcomp::univariate_calpha(), level = level)$confint
        colnames(ci)[1] <- "est"
        ci + os
      },
      band = {
        ci <- confint(multcomp::glht(multcomp::parm(cf_, vc_), linfct = X_),
                      calpha = multcomp::adjusted_calpha(), level = level)$confint
        colnames(ci)[1] <- "est"
        ci + os
      })
    dn <- vector("list", 2)
    names(dn) <- c(rv, "estim")
    dn[[rv]] <- nd[cci, rv]
    dn[[2]] <- colnames(res)
    dimnames(res) <- dn
    res <- switch(type,
                  trafo = res,
                  distribution = object$model$ctm$todistr$p(res),
                  survivor = 1 - object$model$ctm$todistr$p(res),
                  cumhazard = -log(1 - object$model$ctm$todistr$p(res)))
    attr(res, "scale") <- type
    res
  })
  class(out) <- c("terms.tramME", class(out))
  return(out)
}

