## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("survival")
data("sleepstudy", package = "lme4")

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

## NOTE: strange behavior -- accidentally create new objects, that share the same
## param:
mod_lm <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
pr1 <- coef(mod_lm, with_baseline = TRUE)
par1 <- tramME:::.get_par(mod_lm$tmb_obj)
## Expected behavior: the function call creates a copy of the object,
## setting the parameters in the function shouldn't affect the origiginal
fun <- function(obj) {coef(obj) <- c(-1, 0.05, 0.5); logLik(obj)}
chkwarn(fun(mod_lm), "Removing")
## However, originally, the values in mod_lm object also changed,
## because environments (param$env) are not copied but changed in-place.
pr2 <- coef(mod_lm, with_baseline = TRUE)
par2 <- tramME:::.get_par(mod_lm$tmb_obj)
chkeq(pr1, pr2)
chkeq(par1, par2) ## Check the tmb_object separately

## -- Log-likelihood
mod_lm <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
stopifnot(is.na(logLik(mod_lm)))
chkerr(logLik(mod_lm, param = c(-5, -1, 2, 0, 0, 0)), em = "constraints")
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

mod_sr <- SurvregME(Surv(time, status) ~ rx, data = rats)
vc1 <- vcov(mod_sr, method = "numDeriv")
vc2 <- vcov(mod_sr, method = "analytical") ## NOTE: default in this specific model
vc3 <- vcov(mod_sr, method = "optimHess")
chkeq(vc1, vc2)
chkerr(chkeq(vc1, vc3)) ## NOTE: w/ optimHess, it's slightly different

## -- variable names
chkid(variable.names(mod_sr, "grouping"), NULL)
chkid(variable.names(mod_sr, "interacting"), NULL)
chkid(variable.names(mod_sr, "response"), "Surv(time, status)")
mod_sr2 <- SurvregME(Surv(time, status) ~ rx + (1 | litter/rx), data = rats,
                 nofit = TRUE)
chkid(variable.names(mod_sr2, "grouping"), c("rx", "litter"))

## -- VarCorr
chkid(length(VarCorr(mod_sr)), 0L)
chkid(length(VarCorr(mod_sr2)), 2L)

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
pr <- c(coef(mod_lm2, fixed = FALSE, with_baseline = TRUE),
        varcov(mod_lm2, as.theta = TRUE))
re1 <- ranef(mod_lm, param = pr, condVar = TRUE)
re2 <- ranef(mod_lm2, condVar = TRUE)
chkeq(re1, re2)
chkid(ranef(mod_sr), NULL)

nd <- sleepstudy[1:20, ]
re1 <- ranef(mod_lm2, newdata = nd)
re2 <- ranef(mod_lm2, condVar = FALSE)
chkeq(re1$Subject, re2$Subject[1:2, ])

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

## -- print & summary
data("veteran", package = "survival")
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
stopifnot(grepl("Mixed-effects", ss$name, fixed = TRUE))
stopifnot(!ss$fitted)
ss <- summary(mod_lm2)
stopifnot(ss$fitted)

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
par <- mod_polr1$tmb_obj$par
chkeq(logLik(mod_polr1, param = par), logLik(mod_polr2, param = par))

data("eortc", package = "coxme")
dat <- eortc
dat[dat$center <= 10, "center"] <- NA
mm1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat,
                log_first = TRUE, nofit = TRUE) ## na.omit is the default
mm2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = eortc,
                log_first = TRUE, subset = center > 10, nofit = TRUE)
par <- mm1$tmb_obj$par
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
par <- mod_cox2$tmb_obj$par
chkeq(logLik(mod_cox1, param = par), logLik(mod_cox2, param = par))

## -- offsets
os <- runif(nrow(sleepstudy))
mod_lm1 <- Lm(Reaction ~ Days, data = sleepstudy, offset = os)
mod_lm2 <- LmME(Reaction ~ Days, data = sleepstudy)
chkerr(chkeq(logLik(mod_lm1), logLik(mod_lm2), check.attributes = FALSE))
mod_lm2 <- update(mod_lm2, offset = os)
chkeq(logLik(mod_lm1), logLik(mod_lm2), check.attributes = FALSE)

## -- update
## NOTE: When the updated model must have the exact same bases, pass the ctm into update
## (used by e.g. logLik.tramME)
mod_cox1 <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats, log_first = TRUE,
                    order = 12, nofit = TRUE)
mod_cox2 <- update(mod_cox1, data = rats[1:200, ])
chkerr(chkid(mod_cox1$model$ctm, mod_cox2$model$ctm)) ## not the same
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


options(oldopt)
