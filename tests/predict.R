## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
data("sleepstudy", package = "lme4")

# -- newdata w/ response
fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
pr1 <- predict(fit, newdata = sleepstudy[1:20, ], type = "trafo")
pr2 <- predict(fit, type = "trafo")[1:20]
chkid(pr1, pr2)

## -- newdata w/o response
nd <- sleepstudy
nd[["Reaction"]] <- NULL
pr1 <- predict(fit, newdata = nd, ranef = ranef(fit, raw = TRUE), type = "trafo", K = 100)
chkid(unname(dim(pr1)), c(100L, nrow(sleepstudy)))
chkerr(pr2 <- predict(fit, newdata = nd, K = 100))

## -- check "zero" option of ranef
nd <- sleepstudy[c(4, 14, 24), ]
pr1 <- predict(fit, newdata = nd, ranef = rep(0, 6), q = c(300, 320), type = "distribution")
pr2 <- predict(fit, newdata = nd, ranef = "zero", q = c(300, 320), type = "distribution")
chkid(pr1, pr2)

## -- various formatting of ranef
pr1 <- predict(fit, ranef = ranef(fit, raw = TRUE), type = "quantile", prob = 0.5)
## ... also with formatted ranefs
pr2 <- predict(fit, ranef = ranef(fit), type = "quantile", prob = 0.5)
chkid(pr1, pr2)

## -- compare w/ lmer
library("lme4")
lm1 <- LmME(Reaction ~ Days + (Days || Subject), data = sleepstudy)
lm2 <- lmer(Reaction ~ Days + (Days || Subject), data = sleepstudy, REML = FALSE)
chkeq(predict(lm1, type = "lp") * sigma(lm1) + coef(lm1, as.lm = TRUE)[1],
      fitted(lm2), tol = 1e-6)

## -- w/ Surv response & w/o shift terms
library("survival")
fit2 <- SurvregME(Surv(time, status) | rx ~ 1 + (1 | litter), data = rats,
                  dist = "rayleigh")
pdf(file = NULL)
pl <- plot(fit2, newdata = rats[1:2, ], ranef = -1, type = "survivor", K = 100,
           col = 1)
dev.off()
chkeq(dim(pl), c(100L, 2L), check.attributes = FALSE)

## -- compare with tram
m1 <- Lm(dist ~ speed, data = cars, fixed = c(speed = 0.3))
m2 <- LmME(dist ~ speed, data = cars, fixed = c(speed = 0.3))
stopifnot(m2$opt$convergence == 0)
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
nd <- cars
nd[["dist"]] <- NULL
pr1 <- predict(as.mlt(m1), newdata = nd, type = "trafo", K = 200)
pr2 <- predict(m2, newdata = nd, type = "trafo", K = 200)
chkeq(pr1, pr2, tol = 1e-6)

data("neck_pain", package = "ordinalCont")
m1 <- ColrME(vas ~ laser * time, data = neck_pain, bounds = c(0, 1), support = c(0, 1))
m2 <- Colr(vas ~ laser * time, data = neck_pain, bounds = c(0, 1), support = c(0, 1))
stopifnot(m2$opt$convergence == 0)
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
nd <- neck_pain[1, ]
nd$vas <- NULL
pr1 <- predict(m1, newdata = nd, K = 5, type = "distribution")
pr2 <- predict(as.mlt(m2), newdata = nd, K = 5, type = "distribution")
chkeq(pr1, pr2, tol = 1e-5)

## -- edge cases
mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
pr <- predict(mod, ranef = "zero", type = "trafo")
stopifnot(all(is.na(pr)))
coef(mod) <- c(-5, 2, 0.05)
pr <- predict(mod, ranef = "zero", type = "trafo")
stopifnot(all(!is.na(pr)))
chkerr(pr <- predict(mod), "supply random effects")

## -- check random effects
m <- CoxphME(Surv(time, status) ~ rx + (1 | litter), data = rats, log_first = TRUE)
nd <- rats[1, ]
pr <- predict(m, newdata = nd, type = "trafo", ranef = "zero", K = 5)
pr2 <- predict(m, newdata = nd, type = "trafo", ranef = 1, K = 5)
chkeq(pr + 1, pr2)
m <- LmME(Reaction ~ Days + (1 | Subject), data = sleepstudy)
nd <- sleepstudy[10, ]
nd$Reaction <- NULL
pr <- predict(m, newdata = nd, type = "trafo", ranef = "zero", K = 5)
pr2 <- predict(m, newdata = nd, type = "trafo", ranef = 2, K = 5)
chkeq(pr - 2, pr2)


## --- w/ smooth terms
## lp
library("mgcv")
fit_sm1 <- LmME(Reaction ~ s(Days) + (1 | Subject), data = sleepstudy)
fit_sm2 <- gam(Reaction ~ s(Days) + s(Subject, bs = "re"), data = sleepstudy, method = "ML")
## fitted REs
pr1a <- (predict(fit_sm1) - coef(fit_sm1, with_baseline = TRUE)[1]) * sigma(fit_sm1)
pr2 <- predict(fit_sm2)
chkeq(pr1a, c(pr2), tol = 1e-5)
## zero REs
pr1b <- (predict(fit_sm1, ranef = "zero") - coef(fit_sm1, with_baseline = TRUE)[1]) * sigma(fit_sm1)
pr2 <- predict(fit_sm2, exclude = "s(Subject)")
chkeq(pr1b, c(pr2), tol = 1e-5)
## new data
nd <- sleepstudy[1:10,  ]
pr1c <- (predict(fit_sm1, newdata = nd) - coef(fit_sm1, with_baseline = TRUE)[1]) * sigma(fit_sm1)
chkeq(pr1a[1:10], pr1c)
pr2 <- predict(fit_sm2, newdata = nd)
chkeq(pr1c, c(pr2), tol = 1e-5)

## density
data("mcycle", package = "MASS")
m <- LmME(accel ~ s(times, bs = "cr"), data = mcycle)
nd <- data.frame(times = seq(5, 30, by = 5))
qs <- seq(-134, 75, length.out = 50)
lp <- predict(m, nd, type = "lp")
lp2 <- (lp - coef(m, with_baseline = TRUE)[1]) * sigma(m)
pr1 <- predict(m, nd, type = "density", q = qs)
pr2 <- sapply(lp2, function(x) dnorm(qs, mean = x, sd = sigma(m)))
chkeq(pr1, pr2, check.attributes = FALSE)

## a complicated plot
mod_su <- SurvregME(Surv(time, status) | celltype ~ age + s(karno, by = trt),
                    data = veteran, dist = "loglogistic")
pdf(file = NULL)
plot(mod_su, type = "density")
dev.off()

options(oldopt)
