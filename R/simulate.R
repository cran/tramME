##' Simulate from a tramME model
##'
##' Utilizes the simulation method of mlt. When the vector of random effects is supplied, the
##' simulation is conditional on it.
##'
##' In certain settings, the conditional CDF of the outcome cannot be inverted on some
##' some intervals. In these cases, \code{simulate.mlt} returns censored observations.
##' @param object A fitted tramME object.
##' @param ranef If \code{NULL}, random effects are simulated from their estimated
##'   distribution for each of the \code{nsim} draws, i.e. the simulation is from the
##'   marginal/joint distribution of the response (and random effects).
##'   Otherwise the simulation is conditional on the supplied random effects.
##'   When \code{ranef = "zero"}, a vector of zeros with the right size is
##'   substituted.
##' @param what Defaults to \code{"response"}. \code{what = "ranef"} returns draws from
##'   the random effects distribution, \code{what = "joint"} results in simulated data
##'   from the joint distribution of random effects and responses. When it is
##'   set to other than 'response', \code{ranef=NULL} and \code{bysim=TRUE} must
##'   be set.
##' @inheritParams mlt::simulate.mlt
##' @param ... Additional arguments, passed to \code{\link[mlt]{simulate.mlt}}.
##' @return A \code{simulate.tramME} object with the structure defined by the inputs.
##' @examples
##' data("sleepstudy", package = "lme4")
##' fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
##' sim <- simulate(fit, nsim = 10, seed = 123)
##' @importFrom stats simulate runif
##' @importFrom mlt R
##' @export
## TODO: parallel support
## FIXME: strip unnecessary attributes from results returned by mlt
simulate.tramME <- function(object, nsim = 1, seed = NULL,
                           newdata = model.frame(object), ranef = NULL,
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

  ## NOTE: in case ranef is formatted, transfrom it to a vector
  if (is.list(ranef)) {
    ranef <- unname(unlist(lapply(ranef, function(x) c(t(x)))))
  }

  ## --- bysim only when random effects are also simulated
  what <- match.arg(what)
  if (what != "response")
    stopifnot(isTRUE(bysim), is.null(ranef))
  ## --- Dummy ctm model NOTE: REs enter as fixed parameter FEs
  fmlt <- .cctm(object$model$ctm, coef(object, with_baseline = TRUE, fixed = TRUE),
                negative = object$model$negative)
  ## --- REs
  rsiz <- .re_size(attr(object$param, "re")$blocksize, newdata)
  ## FIXME: use model.matrix(, type = "ranef") instead
  re_struct <- re_terms(object$model$ranef, data = newdata,
                        negative = FALSE) ## NOTE: set in .cctm
  if (is.null(ranef)) {
    ## --- Sample from the marginal/joint distribution
    re <- replicate(nsim, .sim_re(varcov(object), n = rsiz$nlev), simplify = FALSE)
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


##' Simulates random effects vector from a tramME object
##' @param vc list of RE variance-covariances
##' @param n list of number of values to be simulated for each grouping factor
.sim_re <- function(vc, n) {
  stopifnot(length(vc) == length(n))
  gamma <- mapply(FUN = function(sig, ns) {
    rn <- mvtnorm::rmvnorm(ns, sigma = sig)
    c(t(rn))
  }, sig = vc, ns = n, SIMPLIFY = FALSE)
  unlist(gamma, use.names = FALSE)
}


parboot <- function(object, ...) {
  UseMethod("parboot")
}


##' Do parametric bootsrap using a tarmME model
##' @param object A \code{tramME} object.
##' @param statistic A function that calculates the statistic of interest.
##' @param nsim Number of draws.
##' @param conditional Logical, if \code{TRUE}, the resampling is conditional on the
##'   fitted vector of random effects.
##' @param ... Optional arguments passed to \code{statistic}.
##' @inheritParams base::sapply
##' @inheritParams simulate.tramME
##' @inheritParams confint.tramME
##' @return A list/vector/array (whichever is consistent with \code{simplify}) of bootstrapped
##'   values returned by \code{statistic}.
##' @export
parboot.tramME <- function(object, statistic, nsim = 1, conditional = FALSE,
                           seed = NULL, ..., simplify = TRUE,
                           parallel = c("no", "multicore", "snow"),
                           ncpus = getOption("profile.ncpus", 1L)) {

  plist <- .parallel_default(parallel, ncpus)

  ## --- Set up random seeds
  if (plist$do_parallel && plist$parallel == "multicore") {
    RNGkind("L'Ecuyer-CMRG")
  }
  if (length(seed) == 1) set.seed(seed)
  if (length(seed) > 1) {
    stopifnot(length(seed) == nsim)
    seeds <- seed ## NOTE: If a vector, pass it to simulate and save
  } else{
    seeds <- seq(nsim)
  }

  dat <- model.frame(object)
  rv <- variable.names(object, "response")
  ip <- list(beta = coef(object, with_baseline = TRUE),
             theta = varcov(object, as.theta = TRUE)) ## NOTE: initial parameters

  if (conditional) {
    re <- ranef(object, raw = TRUE)
  } else {
    re <- NULL
  }

  fun <- function(s) {
    if (length(seed) <= 1) s <- NULL
    out <- NA
    if (!is.null(s)) attr(out, "seed") <- s
    dat[[rv]] <- simulate(object, ranef = re, nsim = 1, seed = s)
    refit <- try(update(object, data = dat, ctm = object$model$ctm,
                        initpar = ip), silent = TRUE)
    if (inherits(refit, "try-error")) {
      attr(out, "status") <- "estimation error"
      return(out)
    }
    if (refit$opt$convergence != 0) {
      attr(out, "status") <- "non-convergence"
      return(out)
    }
    res <- try(statistic(refit, ...), silent = TRUE)
    if (inherits(res, "try-error")) {
      attr(out, "status") <- "statistic error"
      return(out)
    }
    out <- res
    if (!is.null(s)) attr(out, "seed") <- s
    attr(out, "status") <- "ok"
    return(out)
  }

  if (plist$do_parallel) {
      if (plist$parallel == "multicore") {
        boot <- parallel::mclapply(seeds, fun, mc.cores = ncpus)
      } else if (plist$parallel == "snow") {
        ## FIXME: add snow support
        stop("No snow support yet")
      }
  } else {
    boot <- lapply(seeds, fun)
  }

  if (!isFALSE(simplify) && length(boot))
    boot <- simplify2array(boot, higher = (simplify == "array"))
  return(boot)
}
