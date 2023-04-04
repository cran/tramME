## ===== Resp class =====
## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("survival")

dat <- data.frame(x1 = c(100, 1:10, NA, NA), x2 = c(100, 2:10, NA, 11, NA),
                  x3 = c(NA, NA, 0:10))
dat$r <- with(dat, Resp(x1, x2, x3))

with(dat, chkid(Resp(x1, x2), Surv(x1, x2, type = "interval2")))

chkid(length(dat$r), 13L)
chkid(dat[1:3, ]$r, dat$r[1:3])
chkid(is.na(dat$r), c(rep(FALSE, 12), TRUE))

mf <- model.frame(r ~ 1, data = dat, na.action = na.omit)
chkid(nrow(mf), 12L)

m1 <- CoxphME(Resp(x1, x2) ~ 1, data = dat, nofit = TRUE)
chkid(inherits(d1 <- m1$data[[1]], "Surv"), TRUE)
m2 <- CoxphME(Resp(x1, x2, x3) ~ 1, data = dat, nofit = TRUE,
              support = c(0, 20))
## NOTE: support is needed otherwise mlt throws an error (Report to Torsten?)
chkid(inherits(d2 <- m2$data[[1]], "Resp"), TRUE)
chkid(nrow(d1), nrow(d2))
chkerr(CoxphME(Resp(x1, x2, x3) ~ 1, data = dat, support = c(0, 20),
               nofit = TRUE, na.action = na.fail))
## FIXME: wrong error message?!

## === setting bounds
dat$r2 <- with(dat, Resp(x1, x2, x3, bounds = c(0, Inf)))
chkid(print(dat$r[-12]), print(dat$r2[-12]))
chkid(identical(print(dat$r[12]), print(dat$r2[12])), FALSE)

## === adjust in priniting
chkid(identical(print(R(dat$r2)), print(dat$r2)), FALSE)
## but the saved values are adjusted
chkid(all(unclass(dat$r2) > 0, na.rm = TRUE), TRUE)

summarize_tests()

options(oldopt)

