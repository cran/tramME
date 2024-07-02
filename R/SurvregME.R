##' Mixed-effects Additive Parametric Survival Models
##'
##' Estimates various mixed-effects additive parametric models (not exclusively)
##' for survival analysis.
##'
##' @inheritParams tramME
##' @inheritParams tram::Survreg
##' @details
##'
##' The parameterization is slightly different from
##'   \code{\link[survival:survreg]{survival::survreg}}, see Hothorn et al.
##'   (2018).  The results can be transformed back to the \code{survreg}
##'   parameterization with specific methods provided by \code{tramME}.
##'
##' The model extends \code{\link[tram:Survreg]{tram::Survreg}} with random
##'   effects and (optionally penalized) additive terms. For details on
##'   mixed-effect transformation models, see Tamasi and Hothorn (2021).
##'
##' The elements of the linear predictor are parameterized with negative
##'   parameters (i.e. \code{negative = TRUE} in \code{\link[tram]{tram}}).
##' @inherit tramME references
##' @return A \code{SurvregME} model object.
##' @examples
##' library("survival")
##' rats$litter <- factor(rats$litter)
##' m <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats,
##'                dist = "weibull")
##' summary(m)
##' coef(m, as.survreg = TRUE)
##' @importFrom tram Survreg
##' @export
SurvregME <- function(formula, data, subset, weights, offset, na.action = na.omit,
                      dist = c("weibull", "logistic", "gaussian", "exponential",
                               "rayleigh", "loggaussian", "lognormal", "loglogistic"),
                      scale = 0,
                      silent = TRUE, resid = FALSE, do_update = FALSE,
                      estinit = TRUE, initpar = NULL,
                      fixed = NULL, nofit = FALSE,
                      control = optim_control(),
                      ...) {
  cl <- match.call(expand.dots = TRUE)
  args <- as.list(cl[-1L])
  args$call <- cl
  args$tram <- "Survreg"
  args$dist <- match.arg(dist)
  args$scale <- scale
  out <- do.call("tramME", args = args, envir = environment(formula))
  class(out) <- c("SurvregME", class(out))
  return(out)
}


##' Extract the coefficients of the fixed effects terms of an SurvregME model.
##' @param object An \code{SurvregME} object.
##' @param as.survreg If \code{TRUE}, return the transformed coefficients as in a
##'   \code{survival::survreg} object.
##' @inheritParams coef.LmME
##' @return A numeric vector of the transformed coefficients.
##' @examples
##' library("survival")
##' fit <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats)
##' coef(fit, as.survreg = TRUE)
##' @importFrom stats coef
##' @export
coef.SurvregME <- function(object, as.survreg = FALSE, ...) {
  coef.LmME(object, as.lm = as.survreg, ...)
}
