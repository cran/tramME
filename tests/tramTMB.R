## Testing the internal machinery of the package

## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("mlt")

## -- Models can be set up based either on ctm objects or on formula notation
data("sleepstudy", package = "lme4")

yv <- numeric_var("Reaction", support = c(195, 450), bounds = c(0, Inf))
yb <- Bernstein_basis(yv, order = 5, ui = "increasing")
xb <- as.basis(~ Days, data = sleepstudy, remove_intercept = TRUE, negative = TRUE)
ctmod <- ctm(yb, shifting = xb , todistr = "Normal", data = sleepstudy)

mod1 <- tramME_model(ctm = ctmod, formula = ~ (Days | Subject))
mod2 <- tramME_model(Reaction ~ Days + (Days | Subject), data = sleepstudy,
  tram = "BoxCox", order = 5, support = c(195, 450), bounds = c(0, Inf))

chkid(deparse(mod1$ctm), deparse(mod2$ctm))
chkid(mod1$ctm, mod2$ctm, ignore.environment = TRUE)
chkid(mod1$negative, mod2$negative)
chkid(mod1$ranef, mod2$ranef)
## NOTE: formulas are not the same
stopifnot(inherits(mod1$formula, "dummy_formula"))

## -- Random effects are not necessary
mod1 <- tramME_model(ctm = ctmod)
chkid(mod1$ranef, NULL)
mod2 <- tramME_model(Reaction ~ Days, data = sleepstudy, tram = "BoxCox",
                     order = 5, support = c(195, 450), bounds = c(0, Inf))
chkid(deparse(mod1$ctm), deparse(mod2$ctm))

## -- Data can come from the global environment
mod1 <- tramME_model(Reaction ~ Days + (Days | Subject), tram = "Lm", data = sleepstudy)
attach(sleepstudy)
mod2 <- tramME_model(Reaction ~ Days + (Days | Subject), tram = "Lm")
chkid(mod1, mod2, ignore.environment = TRUE)

## --- update
library("survival")
mod <- SurvregME(Surv(time, status) ~ rx + (1 | litter), data = rats, nofit = TRUE)
pr <- mod$tmb_obj$par
chkerr(mod$tmb_obj$resid(pr), "Residuals")
mod_tmb <- update(mod$tmb_obj, resid = TRUE)
res <- mod_tmb$resid(pr)
## TODO: create a test with updating map to see if it really causes an error

## --- copy
obj2 <- duplicate(mod$tmb_obj)
stopifnot(!identical(mod$tmb_obj$env, obj2$env))
chkeq(mod$tmb_obj$fn, obj2$fn, check.environment = FALSE)
stopifnot(!identical(environment(mod$tmb_obj$fn), environment(obj2$fn)))
chkid(environment(environment(obj2$fn)$fn), obj2$env)
