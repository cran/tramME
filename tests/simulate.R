## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
data("sleepstudy", package = "lme4")

## -- output structure
fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
sim <- simulate(fit)
chkid(ifelse(is.numeric(sim), length(sim), nrow(sim)), nrow(sleepstudy))
sim <- simulate(fit, newdata = sleepstudy[1:4, ], nsim = 3)
chkid(length(sim), 3L)
chkid(ifelse(is.vector(sim[[1]]), length(sim[[1]]), nrow(sim[[1]])), 4L)
sim <- simulate(fit, ranef = ranef(fit, raw = TRUE), nsim = 2, bysim = FALSE)
chkid(length(sim), nrow(sleepstudy))
chkid(ifelse(is.vector(sim[[1]]), length(sim[[1]]), nrow(sim[[1]])), 2L)
chkerr(sim <- simulate(fit, newdata = sleepstudy[c(1, 11), ], ranef = c(0.1, -0.1)))
sim <- simulate(fit, newdata = sleepstudy[1, ], ranef = c(0, 0),
                nsim = 10, bysim = FALSE)
chkid(length(sim), 1L)
chkerr(sim <- simulate(fit, nsim = 2, what = "joint", bysim = FALSE))
sim <- simulate(fit, nsim = 2, what = "joint")
chkid(names(sim[[1]]), c("responses", "ranef"))
sim <- simulate(fit, nsim = 1, what = "joint")
chkid(names(sim), c("responses", "ranef"))
sim <- simulate(fit, nsim = 4, what = "ranef")
chkid(length(sim[[4]]), 36L)
fit2 <- LmME(Reaction ~ Days + (1 | Subject), data = sleepstudy)
sim <- simulate(fit2, newdata = sleepstudy[c(1, 11, 21), ], what = "ranef", nsim = 10)
chkid(length(sim[[5]]), 3L)

## -- setting the seed
s1 <- simulate(fit, seed = 123, what = "joint", nsim = 4)
s2 <- simulate(fit, seed = 123, what = "joint", nsim = 4)
chkid(s1, s2, ignore.environment = TRUE)
s2 <- simulate(fit, seed = 124, what = "joint", nsim = 4)
chkerr(chkid(s1, s2, ignore.environment = TRUE))
set.seed(123)
s1 <- simulate(fit, bysim = FALSE, nsim = 4)
set.seed(123)
s2 <- simulate(fit, bysim = FALSE, nsim = 4)
chkid(s1, s2, ignore.environment = TRUE)

## -- check "zero" option for ranef
nd <- sleepstudy
s1 <- simulate(fit, newdata = nd, ranef = rep(0, 2 * nlevels(sleepstudy$Subject)), seed = 10)
s2 <- simulate(fit, newdata = nd, ranef = "zero", seed = 10)
chkid(s1, s2, ignore.environment = TRUE)

options(oldopt)
