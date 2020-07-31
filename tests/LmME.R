library("tramME")
source("test_util.R")
suppressPackageStartupMessages(library("lme4"))

oldopt <- options(digits = 5)
chktol <- function(x, y, tol = sqrt(.Machine$double.eps))
  stopifnot(isTRUE(all.equal(x, y, tol = tol, check.attributes = FALSE)))
## chkwarn <- function(expr, wm) tryCatch(expr, warning = function(w) grepl(wm, w))

data("sleepstudy")

fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
fit2 <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)

## coefficients w/ as.lm
fit3 <- fit
chkwarn(coef(fit3) <- coef(fit, with_baseline = TRUE), "unfitted") ## warning
chktol(coef(fit, as.lm = TRUE), coef(fit3, as.lm  = TRUE))
fit3 <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
chktol(names(coef(fit3, as.lm = TRUE))[1], "(Intercept)")

## compare to lme4
chktol(fixef(fit2), coef(fit, as.lm = TRUE), tol = 1e-6)
chktol(sigma(fit2), sigma(fit), tol = 1e-5)
chktol(logLik(fit), logLik(fit2))
## --- SEs of fixed effects
chktol(sqrt(diag(vcov(fit, as.lm = TRUE, pargroup = "fixef"))),
       sqrt(diag(vcov(fit2))), tol = 1e-4)
## --- Variances and correlations of the random effects
vc1 <- VarCorr(fit, as.lm = TRUE)
vc2 <- lme4::VarCorr(fit2)
chktol(vc1$Subject$var, diag(vc2[[1]]), tol = 1e-4)
chktol(vc1$Subject$corr[1, 2], as.data.frame(vc2)[3, "sdcor"], tol = 1e-4)
## --- Random effects
re1 <- ranef(fit, as.lm = TRUE)
re2 <- lme4::ranef(fit2)
chktol(colnames(re1$Subject), colnames(re2$Subject))
chktol(re1$Subject, re2$Subject, tol = 1e-4)
## --- Confidence intervals for the fixed effects
ci1 <- confint(fit, as.lm = TRUE, pargroup = "fixef")
ci2 <- confint(fit2, 5:6, method = "Wald")
chktol(ci1, ci2, tol = 1e-5)
try(confint(fit, as.lm = TRUE, pargroup = "ranef", type = "profile")) ## error

## internal compatibility of VarCorr.LmME
(vc <- VarCorr(fit, as.lm = TRUE))
th <- fit$tmb_sdr$value[4:6]
var <- exp(th[1:2])^2
cor <- matrix(c(1, th[3], 0, 1), ncol = 2)
cor <- cor %*% t(cor)
s <- sqrt(diag(cor))
cor <- diag(1/s) %*% cor %*% diag(1/s)
chktol(vc$Subject$var, var)
chktol(vc$Subject$corr, cor)

options(oldopt)
