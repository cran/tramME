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
chkerr(chkeq(ip$beta, coef(mod, with_baseline = TRUE), check.attributes = FALSE))
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

options(oldopt)
