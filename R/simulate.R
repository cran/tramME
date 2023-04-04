##' Simulate from a \code{tramME} model
##'
##' @section Warning:
##'
##' This method is under active development and may be subject to change. It is
##'   currently limited to simulating random effects.
##'
##' @param object A \code{tramME} object.
##' @param type Defaults to \code{"ranef"}. Currently the only avalable option.
##' @inheritParams mlt::simulate.mlt
##' @param ... Additional arguments, passed to \code{\link[mlt]{simulate.mlt}}.
##' @return A length \code{nsim} list of draws.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' sim <- simulate(fit, nsim = 10, seed = 123)
##' @importFrom stats simulate runif
##' @export
simulate.tramME <- function(object, nsim = 1, seed = NULL,
                            newdata = model.frame(object),
                            type = c("ranef", "response", "joint"),
                            ...) {
  type <- match.arg(type)
  if (!identical(type, "ranef"))
    stop("Currently only for simulating from the random effects distribution.")

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

  re <- .ranef_trms(object$model$ranef, newdata, TRUE)
  nm <- attr(re, "parnames")
  re <- attr(re, "ranef")
  vc <- varcov(object)
  n <- re$termsize / re$blocksize
  ## -- TODO: this could be made more efficient by drawing all realizations
  ## within a single sim_ranef call
  out <- replicate(nsim, sim_ranef(vc, n, nm), simplify = FALSE)
  ## --

  attr(out, "seed") <- RNGstate
  return(out)
}

## Simulates random effects vector from a tramME object
## @param vc list of RE variance-covariances
## @param n list of number of values to be simulated for each grouping factor
## @param nm ptional name vector for the random effects.
sim_ranef <- function(vc, n, nm = NULL) {
  stopifnot(length(vc) == length(n))
  gamma <- mapply(FUN = function(sig, ns) {
    rn <- mvtnorm::rmvnorm(ns, sigma = sig)
    c(t(rn))
  }, sig = vc, ns = n, SIMPLIFY = FALSE)
  out <- unlist(gamma, use.names = FALSE)
  if (length(nm))
    names(out) <- nm
  out
}
