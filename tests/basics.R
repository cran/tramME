library("tramME")

chk <- function(x, y) stopifnot(isTRUE(all.equal(x, y)))
chkerr <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

## nofit
library("survival")
data("eortc", package = "coxme")
fit <- CoxphME(Surv(y, uncens) ~ trt, data = eortc, log_first = TRUE, nofit = TRUE)
inherits(fit, "ctm")
names(fit)

## setting parameters
data("sleepstudy", package = "lme4")
mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
mod$fitted
chkerr(coef(mod) <- c(1, -1, 0.5))
chkerr(coef(mod) <- c(1, 1))
chkerr(coef(mod) <- c(-1, 0.5, 1)) ## no error
vc <- varcov(mod)
vc[[1]] <- matrix(c(1, 0.2, 0.6, 2), ncol = 2)
chkerr(varcov(mod) <- vc)
vc[[1]] <- diag(1, 3)
chkerr(varcov(mod) <- vc)
vc[[1]] <- matrix(c(1, 0, 0, 2), ncol = 2)
chkerr(varcov(mod) <- vc) ## no error

## missing values and subsets
## --- when a covariate is missing
set.seed(100)
dat1 <- eortc
idx <- sample(nrow(dat1), 200)
dat1[idx, "trt"] <- NA
dat2 <- eortc
dat2$id <- 1
dat2[idx, "id"] <- 0
fit1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat1, log_first = TRUE) ## na.omit is the default
fit2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat2, subset = id > 0,
                log_first = TRUE)
chk(coef(fit1, with_baseline = TRUE), coef(fit2, with_baseline = TRUE))
chk(logLik(fit1), logLik(fit2))
## --- when the response is missing
dat1 <- eortc
dat1[idx, "y"] <- NA
fit1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat1, log_first = TRUE) ## na.omit is the default
chk(coef(fit1, with_baseline = TRUE), coef(fit2, with_baseline = TRUE))
chk(logLik(fit1), logLik(fit2))
## --- when the RE grouping variable is missing
dat1 <- eortc
dat1[dat1$center <= 10, "center"] <- NA
fit1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat1, log_first = TRUE) ## na.omit is the default
fit2 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = eortc, log_first = TRUE, subset = center > 10)
chk(coef(fit1, with_baseline = TRUE), coef(fit2, with_baseline = TRUE))
chk(logLik(fit1), logLik(fit2))
## --- other na.actions
try(fit1 <- CoxphME(Surv(y, uncens) ~ trt + (1 | center), data = dat1,
                    log_first = TRUE, na.action = na.fail))

## optional parameters
data("sleepstudy", package = "lme4")
mod <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE,
                order = 10, log_first = TRUE)
chk(length(coef(mod, with_baseline = TRUE)), 11+1)
attr(mod$model$response$basis, "log_first")

## variable names
variable.names(fit1, which = "response")
variable.names(fit1, which = "shifting")
variable.names(fit1, which = "grouping")
variable.names(fit1)

## outputs/print methods
fit1
summary(fit1)
mod
summary(mod)
VarCorr(mod)

## weights
set.seed(100)
dat1 <- sleepstudy
dat1$w <- sample(10, nrow(sleepstudy), replace = TRUE)
fit1 <- LmME(Reaction ~ Days + (Days || Subject), weights = w, data = dat1)
dat2 <- dat1[rep(seq(nrow(dat1)), dat1$w), ]
fit2 <- LmME(Reaction ~ Days + (Days || Subject), data = dat2)
chk(coef(fit1, with_baseline = TRUE), coef(fit2, with_baseline = TRUE))
chk(logLik(fit1), logLik(fit2))
VarCorr(fit1)
summary(fit1)
logLik(fit1)

## RE structures
mod1 <- LmME(Reaction ~ Days + (Days || Subject), data = sleepstudy, nofit = TRUE)
mod2 <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
chk(lapply(varcov(mod1), dim), list(Subject = c(1L, 1L), Subject = c(1L, 1L)))
chk(lapply(varcov(mod2), dim), list(Subject = c(2L, 2L)))
