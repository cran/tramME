## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("tram")
library("lme4")
gamdat <- mgcv::gamSim(6, n = 500, scale = 2, verbose = FALSE)

## -- Fixed effects only model
m01 <- LmME(Reaction ~ Days, data = sleepstudy, fixed = c("Days" = 0.6))
m02 <- Lm(Reaction ~ Days, data = sleepstudy, fixed = c("Days" = 0.6))
stopifnot(m01$opt$convergence == 0)
stopifnot(m02$convergence == 0)
chkeq(logLik(m01), logLik(m02), check.attributes = FALSE)
chkeq(coef(m01, with_baseline = TRUE), coef(m02, with_baseline = TRUE), tol = 1e-5)
vc1 <- vcov(m01, pargroup = "all")
vc2 <- vcov(m01, pargroup = "all", method = "analytical")
chkeq(vc1, vcov(m02, with_baseline = TRUE), tol = 1e-3)
chkeq(vc2, vcov(m02, with_baseline = TRUE), tol = 1e-5)

## -- Compare with lmer
m1 <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
m2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
stopifnot(m1$opt$convergence == 0)

## logLik, fixed effects, sigma
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
chkeq(fixef(m2), coef(m1, as.lm = TRUE), tol = 1e-6)
chkeq(sigma(m1), sigma(m2), tol = 1e-5)

## --- SEs of fixed effects
vc1 <- vcov(m1, as.lm = TRUE, pargroup = "fixef")
vc2 <- vcov(m2)
chkeq(sqrt(diag(vc1)), sqrt(diag(vc2)), tol = 1e-4, check.attributes = FALSE)

## --- Variances and correlations of the random effects
vc1 <- VarCorr(m1, as.lm = TRUE)
vc2 <- lme4::VarCorr(m2)
chkeq(vc1$Subject$var, diag(vc2[[1]]), tol = 1e-4)
chkeq(vc1$Subject$corr[1, 2], as.data.frame(vc2)[3, "sdcor"], tol = 1e-4)
v <- diag(varcov(m1, as.lm = TRUE)[[1]])
chkeq(v, vc1$Subject$var)
## check original parametrization of the covariance matrix and its as.lm version
th1 <- varcov(m1, as.lm = TRUE, as.theta = TRUE)
th2 <- varcov(m1, as.theta = TRUE)
sig <- sigma(m1)
chkeq(th1, th2 + c(log(sig), log(sig), 0))

## --- Confidence intervals for the fixed effects
ci1 <- confint(m1, as.lm = TRUE, pargroup = "fixef")
ci2 <- confint(m2, 5:6, method = "Wald")
chkeq(ci1, ci2, tol = 1e-5, check.attributes = FALSE)

## --- Random effects
re1 <- ranef(m1, as.lm = TRUE)
re2 <- lme4::ranef(m2)
chkid(colnames(re1$Subject), colnames(re2$Subject))
chkeq(re1$Subject, re2$Subject, tol = 1e-4, check.attributes = FALSE)
pr <- list(beta = coef(m1, with_baseline = TRUE), theta = varcov(m1, as.theta = TRUE))
chkerr(ranef(m1, as.lm = TRUE, param = pr), "not supported")

## --- Residuals
res1 <- resid(m1, as.lm = TRUE)
res2 <- resid(m2)
chkeq(res1, res2, tol = 1e-5)

## --- against gam
mod_gm <- LmME(y ~ s(x0)+ x1 + s(x2) + (1|fac), data = gamdat)
mod_gam <- mgcv::gam(y ~ s(x0)+ x1 + s(x2) + s(fac, bs = "re"), data = gamdat, method = "ML")
vc <- varcov(mod_gm, as.lm = TRUE)
chkid(names(vc), "fac")
vc <- varcov(mod_gm, as.lm = TRUE, full = TRUE)
chkid(names(vc), c("fac", "s(x0)", "s(x2)"))
## NOTE: different parameter structure only check intercept and coef of x1
cf1 <- coef(mod_gm, as.lm = TRUE)[1:2]
cf2 <- coef(mod_gam)[1:2]
chkeq(cf1, cf2, tol = 1e-5)
vc1 <- vcov(mod_gm, as.lm = TRUE, parm = c("(Intercept)", "x1"))
vc2 <- vcov(mod_gam)[1:2, 1:2]
chkeq(vc1, vc2, tol = 1e-2) ## NOTE: would it be more precise with numDeriv?
## logLik
chkeq(as.numeric(logLik(mod_gm)), -mod_gam$gcv.ubre,
      check.attributes = FALSE)

summarize_tests()

options(oldopt)
