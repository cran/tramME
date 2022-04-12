## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("mgcv")
library("tramME")
data("mcycle", package = "MASS")
gamdat <- mgcv::gamSim(6, n = 500, scale = 2, verbose = FALSE)
gamdat2 <- mgcv::gamSim(6, n = 100, scale = 2, verbose = FALSE)

## --
m <- LmME(accel ~ s(times), data = mcycle, nofit = TRUE)
sm <- smooth_terms(m)
chkid(length(sm), 0L)
chkid(class(sm), c("smooth.tramME", "list"))

## -- Check reparametrization
m <- LmME(accel ~ s(times), data = mcycle)
sm <- smooth_terms(m, as.lm = TRUE)
sm2 <- smooth_terms(m)
chkeq(sm[[1]][, 2], sm2[[1]][, 2] * sigma(m))

## -- basic comparison to mgcv w/ MLE
m2 <- mgcv::gam(accel ~ s(times), data = mcycle, method = "ML")
nd <- with(mcycle, data.frame(times = seq(min(times), max(times), length.out = 100)))
sm2 <- predict(m2, newdata = nd, type = "terms")
chkeq(sm[[1]][, 2], sm2[, 1], check.attributes = FALSE, tol = 1e-6)

## -- set up and eval smooth on different data: out-of-sample logLik
ll1 <- logLik(m)
nd <- mcycle[sample(seq(nrow(mcycle))), ]
ll2 <- logLik(m, newdata = nd, type = "integrated") ## only same when no RE is fixed
chkeq(ll1, ll2)

nd <- mcycle[1:100, ]
m3 <- LmME(accel ~ s(times), data = nd, nofit = TRUE)
sm <- tramME:::.tramME_smooth(m) ## set up the smooth term on the whole dataset
m4 <- LmME(accel ~ s(times), data = nd, smooth = sm, nofit = TRUE)
chkerr(chkeq(m3$tmb_obj$env$data$X, m4$tmb_obj$env$data$X,
             tol = 0.1, scale = 1)) ## not the same

mod_gm <- LmME(y ~ s(x0)+ s(x1) + s(x2) + (1|fac), data = gamdat)
## NOTE: by fixing random effects, we change the log-likelihood
chkerr(chkeq(logLik(mod_gm, newdata = gamdat2),
             logLik(mod_gm, newdata = gamdat2, type = "integrated"),
             tol = 0.1, scale = 1)) ## not the same

## -- restrict smooths to be evaluated through newdata argument
chkid(length(sm1 <- smooth_terms(mod_gm)), 3L)
nd <- data.frame(x0 = sm1[[1]]$x0, x2 = sm1[[3]]$x2)
chkid(length(sm2 <- smooth_terms(mod_gm, newdata = nd)), 2L)
chkeq(sm1[c(1, 3)], sm2, check.attributes = FALSE)
nd <- nd[c(1, 50, 100), 1, drop = FALSE]
chkid(length(sm3 <- smooth_terms(mod_gm, newdata = nd)), 1L)
chkeq(sm1[[1]][c(1, 50, 100), ], sm3[[1]], check.attributes = FALSE)

## -- function of variable in smoother
mod <- LmME(accel ~ s(sqrt(times)), data = mcycle)
mod2 <- gam(accel ~ s(sqrt(times)), data = mcycle, method = "ML")
chkeq(-summary(mod2)$sp.criterion, as.numeric(logLik(mod)),
      check.attributes = FALSE, tol = 1e-5)
chkid(variable.names(mod, which = "smooth"), "times")
pr1 <- (predict(mod, newdata = data.frame(times = 1:10), type = "lp") -
  coef(mod, with_baseline = TRUE)[1]) * sigma(mod)
pr2 <- c(predict(mod2, newdata = data.frame(times = 1:10)))
chkeq(pr1, pr2, tol = 1e-5)

## w/ by factor
gamdat <- mgcv::gamSim(4, n = 200)
m1 <- BoxCoxME(y ~ x0 + x1 + s(x2, by = fac, fx = TRUE, k = 8), data = gamdat)
edf1 <- edf_smooth(m1)
m2 <- gam(y ~ x0 + x1 + s(x2, by = fac, fx = TRUE, k = 8), data = gamdat)
edf2 <- summary(m2)$edf
chkeq(edf1, as.vector(edf2), check.attributes = FALSE)
