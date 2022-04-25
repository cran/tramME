## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("survival")
data("sleepstudy", package = "lme4")
set.seed(100)
gamdat <- mgcv::gamSim(6, n = 500, scale = 2, verbose = FALSE)
data("mcycle", package = "MASS")

## -- Set and get parameters
mod_lm <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
chkerr(coef(mod_lm) <- c(1, 1))
coef(mod_lm) <- c(1, -1, 0.5) ## doesn't raise an error until varcov is defined
vc <- varcov(mod_lm)
vc[[1]][] <- diag(2)
chkerr(varcov(mod_lm) <- vc, em = "constraints")
coef(mod_lm) <- c(-1, 0.5, 1) ## no error
varcov(mod_lm) <- vc
chkeq(varcov(mod_lm)$Subject, diag(2), check.attributes = FALSE)
vc[[1]][] <- matrix(c(1, 0.2, 0.6, 2), ncol = 2)
chkerr(varcov(mod_lm) <- vc)

mod_gm <- LmME(y ~ s(x0)+ x1 + s(x2) + (1|fac), data = gamdat, nofit = TRUE)
vc <- varcov(mod_gm)
chkid(names(vc), "fac")
vc <- varcov(mod_gm, full = TRUE)
chkid(names(vc), c("fac", "s(x0)", "s(x2)"))
cf <- coef(mod_gm, with_baseline = TRUE)
chkid(length(cf), 5L) ## NOTE: 2 baseline + 1 shift + 2 smooth

## -- Log-likelihood
mod_lm <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
stopifnot(is.na(logLik(mod_lm)))
chkerr(logLik(mod_lm, param = list(beta = c(-5, -1, 2), theta = c(0, 0, 0))),
       em = "constraints")
mod_lm2 <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
ss2 <- sleepstudy[sample(1:nrow(sleepstudy)), ] ## Just reshuffle
chkeq(logLik(mod_lm2), logLik(mod_lm2, newdata = ss2))
## NOTE: sometimes there are very small numerical differences
## (it might come from the numerical integration)

## -- vcov
stopifnot(all(is.na(vcov(mod_lm))))
chkid(dim(vcov(mod_lm, pargroup = "fixef")), c(3L, 3L))
chkid(dim(vcov(mod_lm, pargroup = "ranef")), c(3L, 3L))
chkid(dim(vcov(mod_lm, pargroup = "shift")), c(1L, 1L))
mod_lm <- update(mod_lm, fixed = c("Days" = 0.5))
chkeq(dim(vcov(mod_lm, pargroup = "shift")), c(0L, 0L))
chkid(rownames(vcov(mod_lm, pargroup = "fixef")), c("(Intercept)", "Reaction"))

if (!.run_test) {
  mod_sr <- SurvregME(Surv(time, status) ~ rx, data = rats)
  vc1 <- vcov(mod_sr, method = "numDeriv")
  vc2 <- vcov(mod_sr, method = "analytical") ## NOTE: default in this specific model
  vc3 <- vcov(mod_sr, method = "optimHess")
  chkeq(vc1, vc2)
  ## NOTE: w/ optimHess, it's slightly different
  chkeq(vc1, vc3, tol = sqrt(.Machine$double.eps), chkdiff = TRUE)
}

mod_gm <- LmME(y ~ s(x0)+ x1 + s(x2) + (1|fac), data = gamdat)
chkid(dim(vcov(mod_gm, pargroup = "smooth")), c(4L, 4L))

## -- variable names
chkid(variable.names(mod_sr, "grouping"), NULL)
chkid(variable.names(mod_sr, "interacting"), NULL)
chkid(variable.names(mod_sr, "smooth"), NULL)
chkid(variable.names(mod_sr, "response"), "Surv(time, status)")
mod_sr2 <- SurvregME(Surv(time, status) ~ rx + (1 | litter/rx), data = rats,
                 nofit = TRUE)
chkid(variable.names(mod_sr2, "grouping"), c("rx", "litter"))

chkid(variable.names(mod_gm, "smooth"), c("x0", "x2"))
chkid(variable.names(mod_gm), c("y", "x1", "x0", "x2", "fac"))
## NOTE: linear shift term comes first (x1)

## -- VarCorr
chkid(length(VarCorr(mod_sr)), 0L)
chkid(length(VarCorr(mod_sr2)), 2L)
chkid(length(VarCorr(mod_gm)), 1L)

## -- confint
ci <- confint(mod_sr, pargroup = "ranef", type = "profile", estimate = TRUE)
chkid(dim(ci), c(0L, 3L))
ci <- confint(mod_sr2)
chkid(dim(ci), c(5L, 2L))
stopifnot(all(is.na(ci)))
ci <- confint(mod_lm, "foo")
chkid(dim(ci), c(0L, 2L))
ci <- confint(mod_lm, parm = "Subject", pmatch = TRUE)
chkid(dim(ci), c(3L, 2L))

m03 <- LmME(dist ~ speed, data = cars)
m04 <- Lm(dist ~ speed, data = cars)
chkeq(confint(m03, pargroup = "shift"), confint(m04), tol = 1e-5,
      check.attributes = FALSE)

## -- random effects
mod_lm <- update(mod_lm, fixed = NULL)
stopifnot(all(is.na(ranef(mod_lm)[[1]])))
pr <- list(beta = coef(mod_lm2, fixed = FALSE, with_baseline = TRUE),
           theta = varcov(mod_lm2, as.theta = TRUE))
re1 <- ranef(mod_lm, param = pr, condVar = TRUE)
re2 <- ranef(mod_lm2, condVar = TRUE)
chkeq(re1, re2)
chkid(ranef(mod_sr), NULL)

nd <- sleepstudy[1:20, ]
re1 <- ranef(mod_lm2, newdata = nd)
re2 <- ranef(mod_lm2, condVar = FALSE)
chkeq(re1$Subject, re2$Subject[1:2, ])

re1 <- ranef(mod_gm, raw = TRUE)
re2 <- ranef(mod_gm, condVar = TRUE)
chkid(re2$fac[[1]], re1[1:4])
chkid(dim(attr(re2$fac, "condVar")), c(4L, 1L))

## fixing smooth terms, or parts
mod_gm2 <- LmME(accel ~ s(times), data = mcycle)
re1 <- ranef(mod_gm2, raw = TRUE)
re2 <- ranef(mod_gm2, raw = TRUE, newdata = mcycle[1:2, ]) ## fix_smooth is on by default
chkeq(re1, re2)
pr <- list(beta = coef(mod_gm2, fixed = FALSE, with_baseline = TRUE),
           theta = varcov(mod_gm2, as.theta = TRUE, full = TRUE))
pr$gamma <- mod_gm2$param$gamma[1]
re3 <- ranef(mod_gm2, param = pr, raw = TRUE)
chkeq(re1, re3)

pr <- list(beta = coef(mod_gm, fixed = FALSE, with_baseline = TRUE),
           theta = varcov(mod_gm, as.theta = TRUE, full = TRUE))
pr$gamma <- mod_gm$param$gamma[1]
re <- ranef(mod_gm, param = pr, condVar = TRUE, fix_smooth = TRUE)
chkid(is.na(attr(re$fac, "condVar")[[1]]), c(TRUE, rep(FALSE, 3)))
chkeq(re[[1]], ranef(mod_gm)[[1]], check.attributes = FALSE)

## -- Residuals
library("survival")
mod_sr <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats,
                    support = c(1, 90))
stopifnot(mod_sr$opt$convergence == 0)
res1 <- resid(mod_sr, newdata = rats[1:30, ])
res2 <- resid(mod_sr)[1:30]
## NOTE: if the smaller sample doesn't contain complete groups, the returned
## residuals will be different for that specific litter. This is because tramME
## refits the random effects when it creates a new object (as it happens with
## newdata)
chkeq(res1, res2)

if (!.run_test) {
  res1 <- resid(mod_gm, fix_smooth = TRUE, newdata = subset(gamdat, subset = fac == 1))
  res2 <- resid(mod_gm, fix_smooth = FALSE)[gamdat$fac == 1]
  chkeq(res1, res2)
  res1 <- resid(mod_gm, fix_smooth = FALSE, newdata = subset(gamdat, subset = fac == 1))
  chkeq(res1, res2, tol = sqrt(.Machine$double.eps), chkdiff = TRUE)
}

mod_gm_bc <- BoxCoxME(y ~ s(x0)+ x1 + s(x2) + (1|fac), data = gamdat)
res1 <- resid(mod_gm_bc, fix_smooth = TRUE)
res2 <- resid(mod_gm_bc, fix_smooth = FALSE)
chkeq(res1, res2)

pr <- list(gamma = mod_gm_bc$param$gamma[1:4])
res1 <- resid(mod_gm_bc, param = pr, newdata = gamdat)
res2 <- resid(mod_gm_bc)
chkeq(res1, res2)

## -- FIXME: the tests below fail! Why?
## probably because of the non-linearity of the log-likelihood
## the derivative of the integrated ll != the derivative of the
## penalized ll wrt the constant
## pr <- list(gamma = mod_sr$param$gamma[1:2])
## res1 <- resid(mod_sr, param = pr, newdata = rats[1:4, ])
## res2 <- resid(mod_sr)[1:4]

## pr <- list(gamma = mod_sr$param$gamma)
## res1 <- resid(mod_sr, param = pr)
## res2 <- resid(mod_sr)
## all.equal(res1, res2)

## m <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = eortc, log_first = TRUE)
## pr <- list(gamma = m$param$gamma[1:2])
## res1 <- resid(m, param = pr, newdata = eortc[c(1, 65), ])
## res2 <- resid(m)[c(1, 65)]
## all.equal(res1, res2, tol = 1e-3)
## --

m <- CoxphME(Surv(time, status) | celltype ~ trt + s(age) + s(karno), data = veteran,
             log_first = TRUE)
res1 <- resid(m, fix_smooth = TRUE)
res2 <- resid(m, fix_smooth = FALSE)
chkeq(res1, res2, tol = 1e-5)

## -- print & summary
mod_sr3 <- SurvregME(Surv(time, status) | celltype ~ trt + age + karno, data = veteran,
                     dist = "loglogistic", fixed = c("age" = 0.02))
stopifnot(mod_sr3$opt$convergence == 0)
mod_sr4 <- Survreg(Surv(time, status) | celltype ~ trt + age + karno, data = veteran,
                   dist = "loglogistic", fixed = c("age" = 0.02))
chkeq(logLik(mod_sr3), logLik(mod_sr4), check.attributes = FALSE)
ss <- summary(mod_sr3)
stopifnot(grepl("Stratified", ss$name, fixed = TRUE))
chkid(ss$fixed, c("age" = 0.02))
ss <- summary(mod_lm)
stopifnot(grepl("Mixed-Effects", ss$name, fixed = TRUE))
stopifnot(!ss$fitted)
ss <- summary(mod_lm2)
stopifnot(ss$fitted)
f <- dist ~ speed
mm <- LmME(f, data = cars)
chkid(summary(mm)$formula, f)

## -- subsets and na.actions
data("soup", package = "ordinal")
dat <- soup
dat$RESP[dat$AGEGROUP == "18-30"] <- NA
chkerr(mod_polr1 <- PolrME(SURENESS | SOUPFREQ ~ PROD + (1 | RESP/PROD),
                           data = dat, nofit = TRUE, na.action = na.fail))
mod_polr1 <- PolrME(SURENESS | SOUPFREQ ~ PROD + (1 | RESP/PROD),
                    data = dat, nofit = TRUE, na.action = na.omit)
mod_polr2 <- PolrME(SURENESS | SOUPFREQ ~ PROD + (1 | RESP/PROD),
                    data = soup, nofit = TRUE, na.action = na.fail,
                    subset = AGEGROUP != "18-30")
par <- list(beta = coef(mod_polr1, with_baseline = TRUE),
            theta = varcov(mod_polr1, as.theta = TRUE))
chkeq(logLik(mod_polr1, param = par), logLik(mod_polr2, param = par))

data("eortc", package = "coxme")
dat <- eortc
dat[dat$center <= 10, "center"] <- NA
mm1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat,
                log_first = TRUE, nofit = TRUE) ## na.omit is the default
mm2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = eortc,
                log_first = TRUE, subset = center > 10, nofit = TRUE)
par <- list(beta = coef(mm1, with_baseline = TRUE),
            theta = varcov(mm1, as.theta = TRUE))
chkeq(logLik(mm1, param = par), logLik(mm2, param = par))

## -- weights & offsets
## NOTE: many of this functionality is disabled at the moment
## data("eortc", package = "coxme")
## mod_cox1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = eortc,
##                     log_first = TRUE, order = 10, nofit = TRUE, ## w/o do_update = TRUE!
##                     support = c(1, 2500)) ## NOTE: set support explicitly to define same bases
## stopifnot(is.null(weights(mod_cox1)))
## we <- sample(c(1, 3), nrow(eortc), replace = TRUE)
## chkerr(weights(mod_cox1) <- we, "do_update")
## mod_cox1 <- update(mod_cox1, do_update = TRUE)
## weights(mod_cox1) <- we
## chkid(weights(mod_cox1), we)
## mod_cox2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = eortc,
##                     log_first = TRUE, order = 10, nofit = TRUE, weights = we,
##                     support = c(1, 2500))
## par <- mod_cox2$tmb_obj$par
## chkeq(logLik(mod_cox1, param = par), logLik(mod_cox2, param = par))
## dat <- eortc[rep(1:nrow(eortc), we), ]
## mod_cox3 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = dat,
##                     log_first = TRUE, order = 10, nofit = TRUE,
##                     support = c(1, 2500))
## chkeq(logLik(mod_cox1, param = par), logLik(mod_cox3, param = par))

## subsequently updated weights and offsets are carried forward...
## mod_cox1 <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats, log_first = TRUE,
##                     order = 12, nofit = TRUE, do_update = TRUE)
## offset(mod_cox1) <- rep(0.1, nrow(rats))
## mod_cox2 <- update(mod_cox1, resid = TRUE)
## chkid(offset(mod_cox1), offset(mod_cox2))
## ## ...but will lead to errors when the data changes its size (as expected)
## chkerr(mod_cox3 <- update(mod_cox1, data = rats[1:200, ]), "differing number of rows")

## -- weights
data("eortc", package = "coxme")
we <- sample(c(1, 3), nrow(eortc), replace = TRUE)
mod_cox1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = eortc,
                    log_first = TRUE, order = 10, nofit = TRUE, weights = we,
                    support = c(1, 2500)) ## NOTE: set support explicitly to define same bases
dat <- eortc[rep(1:nrow(eortc), we), ]
mod_cox2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = dat,
                    log_first = TRUE, order = 10, nofit = TRUE,
                    support = c(1, 2500))
par <- list(beta = coef(mod_cox2, with_baseline = TRUE),
            theta = varcov(mod_cox2, as.theta = TRUE))
chkeq(logLik(mod_cox1, param = par), logLik(mod_cox2, param = par))

## -- offsets
os <- runif(nrow(sleepstudy))
mod_lm1 <- Lm(Reaction ~ Days, data = sleepstudy, offset = os)
mod_lm2 <- LmME(Reaction ~ Days, data = sleepstudy)
chkeq(logLik(mod_lm1), logLik(mod_lm2), check.attributes = FALSE,
      tol = 0.1, scale = 1, chkdiff = TRUE)
mod_lm2 <- update(mod_lm2, offset = os)
chkeq(logLik(mod_lm1), logLik(mod_lm2), check.attributes = FALSE)

## -- update
## NOTE: When the updated model must have the exact same bases, pass the ctm into update
## (used by e.g. logLik.tramME)
mod_cox1 <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats, log_first = TRUE,
                    order = 12, nofit = TRUE)
mod_cox2 <- update(mod_cox1, data = rats[1:200, ])
chkid(mod_cox1$model$ctm, mod_cox2$model$ctm, chkdiff = TRUE) ## not the same
mod_cox2 <- update(mod_cox1, data = rats[1:200, ], ctm = mod_cox1$model$ctm)
chkid(mod_cox1$model$ctm, mod_cox2$model$ctm) ## same

## -- fitmod
data("neck_pain", package = "ordinalCont")
mod_colr <- ColrME(vas ~ time * laser + (1 | id), data = neck_pain, bounds = c(0, 1),
                   support = c(0, 1), order = 4, nofit = TRUE)
fit <- fitmod(mod_colr)
## NOTE: they do not share the environment in the tramTMB
stopifnot(!identical(mod_colr$tmb_obj$env, fit$tmb_obj$env))
fit2 <- ColrME(vas ~ time * laser + (1 | id), data = neck_pain, bounds = c(0, 1),
               support = c(0, 1), order = 4)
chkeq(logLik(fit), logLik(fit2))

data("mcycle", package = "MASS")
m <- LmME(accel ~ s(times), data = mcycle, nofit = TRUE)
f1 <- fitmod(m)
f2 <- LmME(accel ~ s(times), data = mcycle)
chkeq(f1$param, f2$param, tol = 1e-4) ## NOTE: not exactly equal bec of different starting values

## -- model.frame
mod_cox3 <- CoxphME(Surv(time, status) | celltype ~ trt + s(age) + karno,
                    data = veteran, log_first = TRUE, nofit = TRUE)
chkeq(model.frame(mod_cox3), mod_cox3$data, check.attributes = FALSE)
chkeq(model.frame(mod_cox3, data = veteran[1:10, ]), mod_cox3$data[1:10, ],
      check.attributes = FALSE)
chkeq(model.frame(mod_cox3, data = veteran[1:10, ], subset = karno < 60),
      subset(mod_cox3$data[1:10, ], subset = karno < 60),
      check.attributes = FALSE)
chkeq(model.frame(mod_cox3, data = veteran, subset = karno < 60),
      model.frame(mod_cox3, data = mod_cox3$data, subset = karno < 60),
      check.attributes = FALSE)
chkeq(model.frame(mod_cox3, subset = karno < 60),
      subset(mod_cox3$data, subset = karno < 60),
      check.attributes = FALSE)
chkeq(model.frame(mod_cox1, subset = litter <= 3),
      subset(mod_cox1$data, subset = litter <= 3), check.attributes = FALSE)
chkeq(model.frame(mod_cox3, data = veteran, subset = time > 60)[[1]][, 1],
      veteran$time[veteran$time > 60])
chkeq(nlevels(model.frame(mod_lm,
                          subset = as.numeric(as.character(Subject)) > 310,
                          drop.unused.levels = TRUE)$Subject),
      nlevels(sleepstudy$Subject) - 3L)
## w/ offset
st <- sleepstudy
st$foo <- runif(nrow(st))
mod_os <- update(mod_lm, offset = -log(foo), data = st)
chkeq(mf <- model.frame(mod_os, data = st[1:10, ]), mod_os$data[1:10, ])
chkid("(offset)" %in% colnames(mf), TRUE)
chkerr(model.frame(mod_os, data = sleepstudy[1:10, ]),
       "'foo' not found")
chkeq(model.frame(mod_os), model.frame(mod_os, data = st))
chkeq(model.offset(model.frame(mod_os, subset = Reaction < 250)),
      -log(st$foo[st$Reaction < 250]))

## -- model.matrix
nd <- model.frame(mod_cox3)[rep(1, 100), ]
nd[[1]] <- seq(1, 120, length.out = 100)
mm1 <- model.matrix(mod_cox3, data = nd, simplify = TRUE)
mm2 <- model.matrix(mod_cox3, data = nd, simplify = TRUE, keep_sign = FALSE)
chkid(mm1, mm2) ## equal in the case of CoxphME
mm1 <- model.matrix(mod_cox3, data = veteran, subset = karno > 40)
mm2 <- model.matrix(mod_cox3, data = model.frame(veteran, subset = karno > 40))
chkid(mm1, mm2)
nd <- model.frame(mod_lm)[rep(1, 100), ]
nd[[1]] <- seq(150, 250, length.out = 100)
mm1 <- model.matrix(mod_lm, data = nd, simplify = TRUE)
mm2 <- model.matrix(mod_lm, data = nd, simplify = TRUE, drop_unused_groups = TRUE)
chkid(dim(mm1$Zt), c(36L, 100L))
chkid(dim(mm2$Zt), c(2L, 100L))

## -- Anova
## NOTE: this should be at the end because ordinal will mask ranef and VarCorr,
## which can cause problems
library("ordinal")
fit1a <- PolrME(rating ~ temp + contact + (1 | judge), data = wine, method = "probit")
fit2a <- clmm2(rating ~ temp + contact, random = judge, data = wine,
               Hess = TRUE, nAGQ = 1, link = "probit")
chkeq(logLik(fit1a), logLik(fit2a), tol = 1e-7, check.attributes = FALSE)
fit1b <- PolrME(rating | contact ~ temp + (1 | judge), data = wine, method = "probit")
fit2b <- clmm2(rating ~ temp, nominal = ~ contact, random = judge, data = wine,
               Hess = TRUE, nAGQ = 1, link = "probit")
lrt1 <- anova(fit1a, fit1b)
lrt2 <- anova(fit2a, fit2b)
chkeq(lrt1$Chisq[2], lrt2$`LR stat.`[2], tol = 1e-5)

summarize_tests()

options(oldopt)
