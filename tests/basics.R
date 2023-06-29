## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")

## -- Model setup with initpar
ip <- list(beta = c(1.5, 0.08, 0.2))
mod <- LmME(dist ~ speed, data = cars, initpar = ip, nofit = TRUE)
## NOTE: initpars are not set as actual model parameters...
chkeq(ip$beta, coef(mod, with_baseline = TRUE), check.attributes = FALSE,
      tol = 0.1, scale = 1, chkdiff = TRUE)
## ... but the tramTMB object is set up with using them
chkeq(ip$beta, mod$tmb_obj$env$par_checked, check.attributes = FALSE)

## -- Data can come from the gobal environment
data("sleepstudy", package = "lme4")
fit_lm1 <- LmME(Reaction ~ Days + (Days || Subject), data = sleepstudy)
attach(sleepstudy)
fit_lm2 <- LmME(Reaction ~ Days + (Days || Subject))
chkeq(logLik(fit_lm1), logLik(fit_lm2))

## -- Check .th2vc and .vc2th helper functions
library("survival")
mod <- CoxphME(
  Surv(tstart, tstop, status) ~ treat + age + weight + height + (age + weight + height |id),
  data = cgd, log_first = TRUE, order = 5, nofit = TRUE)
pr <- mod$tmb_obj$env$last.par
th <- runif(sum(names(pr) == "theta"))
pr[names(pr) == "theta"] <- th
vc1 <- mod$tmb_obj$report(pr) ## NOTE: using REPORT from TMB
vc1 <- diag(vc1$sd_rep[[1]]) %*% vc1$corr_rep[[1]] %*% diag(vc1$sd_rep[[1]])
rs <- attr(mod$param, "re")
vc2 <- tramME:::.th2vc(th, rs$blocksize)
chkeq(vc1, vc2[[1]])
th2 <- tramME:::.vc2th(vc2, rs$blocksize) ## NOTE: check back-transformation
chkeq(th, th2, check.attributes = FALSE)

## -- Setting up models flexibly: TODO: extend a little bit with simulated data
m1 <- BoxCoxME(dist ~ s(speed), data = cars)
m2a <- BoxCox(dist ~ 1, data = cars, model_only = TRUE)
m2 <- tramME(ctm = m2a, formula = ~ s(speed), data = cars)
chkeq(logLik(m1), logLik(m2))

## -- Call models from within lapply, when some inputs are defined in the function
## both for LmME and tramME
chkid({
  lapply(1, function(i) {
    dat <- cars
    inherits(LmME(dist ~ speed, data = dat), "LmME")
  })[[1]]}, TRUE)

chkid({
  lapply(1, function(i) {
    dat <- cars
    m1 <- BoxCox(dist ~ 1, data = dat, model_only = TRUE)
    m2 <- tramME(ctm = m1, formula = ~ s(speed), data = dat)
    inherits(m2, "tramME")
  })[[1]]}, TRUE)

summarize_tests()

options(oldopt)
