##' Simulate outcome variable from an estimated model
##'
##' Utilizes the simulation method of mlt. When the vector of random effects is supplied, the
##' simulation is conditional on it.
##' @param object A fitted tramME object.
##' @param ranef If NULL, random effects are simulated from their estimated
##'   distribution for each draw in nsim, i.e. the simulation is from the
##'   marginal/joint distribution of the response (and random effects).
##'   Otherwise the simulation is conditional on the supplied random effects.
##'   When \code{ranef = "zero"}, a vector of zeros with the right size is
##'   substituted.
##' @param what Defaults to \code{'response'}. \code{'ranef'} returns draws from
##'   the random effects distribution, \code{'joint'} results in simulated data
##'   from the joint distribution of random effects and responses. When it is
##'   set to other than 'response', \code{ranef=NULL} and \code{bysim=TRUE} must
##'   be set.
##' @inheritParams mlt::simulate.mlt
##' @param ... Additional arguments, passed to \code{\link[mlt]{simulate.mlt}}.
##' @return A simulate.tramME object with the structure defined by the inputs.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' sim <- simulate(fit, nsim = 10, seed = 123)
##' @importFrom stats simulate runif
##' @importFrom mlt R
##' @export
## TODO: parallel support
simulate.tramME <- function(object, nsim = 1, seed = NULL,
                           newdata = NULL, ranef = NULL,
                           what = c("response", "ranef", "joint"),
                           bysim = TRUE, ...) {

  ## --- Seed: from stats:::simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ## --- Unfitted: checking inputs
  if (!object$fitted) {
    if (is.null(newdata))
      stop("newdata must be supplied for unfitted models")
    if (any(is.na(object$pars$coef)))
      stop("Coefficient values must be set before calling predict")
    if (any(is.na(unlist(object$pars$varcov))) && is.null(ranef))
      stop(paste(c("Random effects parameters must be set for",
                   "simulations from the joint/marginal distributions")))
  }
  ## --- bysim only when random effects are also simulated
  what <- match.arg(what)
  if (what != "response")
    stopifnot(isTRUE(bysim), is.null(ranef))
  ## --- Data
  if (is.null(newdata))
    newdata <- object$data$mf
  ## --- Dummy ctm model NOTE: REs enter as fixed parameter FEs
  fmlt <- .dummy_ctm(object$model, coef(object, with_baseline = TRUE))
  ## --- REs
  rsiz <- .re_size(object$model, newdata) ## required size of the RE vector
  re_struct <- .re_data(object$call$formula[-2L], data = newdata,
                        negative = FALSE) ## NOTE: set in .dummy_ctm
  if (is.null(ranef)) {
    ## --- Sample from the marginal/joint distribution
    re <- replicate(nsim, .sim_re(object$pars$varcov, n = rsiz$nlev), simplify = FALSE)
    if (what == "ranef")
      return(re)
    if (nsim > 1) {
      if (bysim) {
        out <- vector(mode = "list", length = nsim)
      } else {
        out <- replicate(nrow(newdata), vector(mode = "list", length = nsim),
                         simplify = FALSE)
      }
    }
    for (i in seq(nsim)) {
      newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, re[[i]]))
      smpl <- simulate(fmlt, newdata = newdata, nsim = 1, bysim = TRUE, ...)
      if (nsim > 1) {
        if (bysim) {
          out[[i]] <- smpl
        } else { ## NOTE: there might be a more elegant (and faster) way of doing this
          if (is.data.frame(smpl)) {
            for (j in seq(nrow(smpl))) {
              out[[j]][[i]] <- smpl[j, ]
            }
          } else {
            for (j in seq(length(smpl))) {
              out[[j]][[i]] <- smpl[j]
            }
          }
        }
      } else {
        out <- smpl
      }
    }
    if (!bysim && (nsim > 1)) { ## Join elements of lists when not bysim
      for (i in seq(length(out))) {
        if (any(sapply(out[[i]], is.data.frame))) {
          out[[i]] <- do.call("rbind", lapply(out[[i]], R))
        } else {
          out[[i]] <- unlist(out[[i]], use.names = FALSE)
        }
      }
    }
  } else {
    ## --- Sample from the conditional distribution
    if (identical(ranef, "zero"))
      ranef <- rep(0, sum(rsiz$bsize * rsiz$nlev))
    re <- ranef
    stopifnot(sum(rsiz$bsize * rsiz$nlev) == length(ranef)) ## check if the supplied RE vector has the right size
    newdata$re_ <- as.numeric(Matrix::crossprod(re_struct$Zt, re))
    out <- simulate(fmlt, newdata = newdata, nsim = nsim, bysim = bysim, ...)
  }
  if (what == "joint") {
    if (inherits(out, "response") && (length(re) == 1)) {
      out <- list(responses = out, ranef = unlist(re))
    } else {
      out <- mapply(function(y, g) list(responses = y, ranef = g),
                    y = out, g = re, SIMPLIFY = FALSE)
    }
  }
  attr(out, "seed") <- RNGstate
  class(out) <- c("simulate.tramME", class(out))
  return(out)
}


##' Print method for \code{simulate.tramME} objects
##'
##' Automatically hides the seed attribute of the object.
##' @param x A \code{simulate.tramME} object
##' @param suppress_seed Logical, suppress seed if true.
##' @param ... Additional parameters passed to various print methods.
##' @return The input simulate.tramME object, invisibly.
##' @export
print.simulate.tramME <- function(x, suppress_seed = TRUE, ...) {
  pr <- x
  if (suppress_seed)
    attr(pr, "seed") <- NULL
  class(pr) <- class(pr)[-1L]
  print(pr, ...)
  invisible(x)
}


##' Refit the model with a new response vector
##'
##' Useful for parametric bootstrap.
##' @param object A tramME object.
##' @param newresp A vector of new response values.
##' @param ... optional arguments for compatibility
##' @importFrom lme4 refit
refit.tramME <- function(object, newresp, ...) {
  fc <- object$call
  fc$formula[[2L]] <- quote(response_)
  newdata <- object$data$mf
  newdata$response_ <- newresp
  fc$weights <- quote(weights)
  fc$offset <- quote(offset)
  fc$data <- quote(newdata)
  return(eval(fc, envir = object$data))
}


##' Simulates random effects vector from a tramME object
##' @param vc list of RE variance-covariances
##' @param n list of number of values to be simulated for each grouping factor
##' @importFrom stats rnorm
.sim_re <- function(vc, n) {
  stopifnot(length(vc) == length(n))
  gamma <- mapply(FUN = function(sig, ns) {
    stopifnot(isSymmetric(sig, tol = sqrt(.Machine$double.eps)))
    ev <- eigen(sig, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1L]))) {
      stop("Covariance matrix is not positive semidefinite")
    }
    r <- matrix(rnorm(nrow(sig) * ns), ncol = ns)
    c(ev$vectors %*%
      crossprod(diag(sqrt(pmax(ev$values, 0)), nrow(sig)), r))
  }, sig = vc, ns = n, SIMPLIFY = FALSE)
  unlist(gamma, use.names = FALSE)
}
