library("tramME")

chk <- function(x, y, ...) stopifnot(isTRUE(all.equal(x, y, ...)))
chkerr <- function(x) inherits(try(x), "try-error")

## structure of the output
data("sleepstudy", package = "lme4")
fit <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
sim <- simulate(fit)
chk(ifelse(is.numeric(sim), length(sim), nrow(sim)), nrow(sleepstudy))
sim <- simulate(fit, newdata = sleepstudy[1:4, ], nsim = 3)
chk(length(sim), 3)
chk(ifelse(is.vector(sim[[1]]), length(sim[[1]]), nrow(sim[[1]])), 4)
sim <- simulate(fit, ranef = ranef(fit, raw = TRUE), nsim = 2, bysim = FALSE)
chk(length(sim), nrow(sleepstudy))
chk(ifelse(is.vector(sim[[1]]), length(sim[[1]]), nrow(sim[[1]])), 2)
chkerr(sim <- simulate(fit, newdata = sleepstudy[c(1, 11), ], ranef = c(0.1, -0.1)))
sim <- simulate(fit, newdata = sleepstudy[1, ], ranef = c(0, 0),
                nsim = 10, bysim = FALSE)
chk(length(sim), 1)
chkerr(sim <- simulate(fit, nsim = 2, what = "joint", bysim = FALSE))
sim <- simulate(fit, nsim = 2, what = "joint")
chk(names(sim[[1]]), c("responses", "ranef"))
sim <- simulate(fit, nsim = 4, what = "ranef")
chk(length(sim[[4]]), 36)
fit2 <- LmME(Reaction ~ Days + (1 | Subject), data = sleepstudy)
sim <- simulate(fit2, newdata = sleepstudy[c(1, 11, 21), ], what = "ranef", nsim = 10)
chk(length(sim[[5]]), 3)

## setting the seed
s1 <- simulate(fit, seed = 123, what = "joint", nsim = 4)
s2 <- simulate(fit, seed = 123, what = "joint", nsim = 4)
chk(s1, s2)
set.seed(123)
s1 <- simulate(fit, bysim = FALSE, nsim = 4)
set.seed(123)
s2 <- simulate(fit, bysim = FALSE, nsim = 4)
chk(s1, s2)

## check "zero" option for ranef
nd <- sleepstudy
s1 <- simulate(fit, newdata = nd, ranef = rep(0, 2 * nlevels(sleepstudy$Subject)), seed = 10)
s2 <- simulate(fit, newdata = nd, ranef = "zero", seed = 10)
chk(s1, s2)

## simulating from an unfitted model
mod <- LmME(Reaction ~ Days + (Days | Subject), data = sleepstudy, nofit = TRUE)
coef(mod) <- c(-10, 0.05, 0.6)
nd <- sleepstudy
nd$Reaction <- NULL
s1 <- simulate(mod, newdata = nd, ranef = rep(0, 2 * nlevels(sleepstudy$Subject)))
chkerr(s1 <- simulate(mod, newdata = sleepstudy))
vc <- varcov(mod)
vc[[1]] <- diag(c(1, 0.05))
varcov(mod) <- vc
s1 <- simulate(mod, newdata = sleepstudy)

## basic properties of simulated data
data("wine", package = "ordinal")
fit <- PolrME(rating | contact ~ temp + (1 | judge), data = wine)
nr1 <- simulate(fit)
nr2 <- simulate(fit)
isTRUE(all.equal(nr1, nr2, check.attributes = FALSE)) ## FALSE
chk(levels(nr1), levels(wine$rating))
is.ordered(nr1)
