## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("tram")
library("survival")

m1 <- SurvregME(Surv(time, status) ~ rx, data = rats, dist = "exponential")
m2 <- survreg(Surv(time, status) ~ rx, data = rats, dist = "exponential")
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
chkeq(coef(m1, as.survreg = TRUE), coef(m2), tol = 1e-6)

m1 <- SurvregME(Surv(time, status) ~ rx, data = rats, dist = "lognormal",
                fixed = c("rx" = -0.4))
m2 <- Survreg(Surv(time, status) ~ rx, data = rats, dist = "lognormal",
                fixed = c("rx" = -0.4))
stopifnot(m1$opt$convergence == 0, m2$convergence == 0)
chkeq(logLik(m1), logLik(m2), check.attributes = FALSE)
chkeq(coef(m1, with_baseline = TRUE), coef(m2, with_baseline = TRUE), tol = 1e-4)

options(oldopt)
