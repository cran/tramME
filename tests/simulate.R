## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")

##
dat <- data.frame(g1 = rep(1:40, each = 50),
                  g2 = rep(1:400, each = 5),
                  x  = runif(2000))
v1 <- matrix(c(1.2, -0.6, -0.6, 0.8), nrow = 2)
u1 <- mvtnorm::rmvnorm(40, sigma = v1)
v2 <- 0.4
u2 <- rnorm(400, sd = sqrt(v2))
dat$y <- with(dat, sin(pi * x)
  + rowSums(u1[rep(1:40, each = 50), ] * cbind(1, x))
  + rep(u2, each = 5)
  + rnorm(2000, sd = 1)
  )

fit <- BoxCoxME(y ~ s(x) + (x | g1) + (1 | g2), data = dat)

## -- setting the seed
s1 <- simulate(fit, seed = 123, nsim = 4)
s2 <- simulate(fit, seed = 123, nsim = 4)
chkid(s1, s2, ignore.environment = TRUE)
s2 <- simulate(fit, seed = 124, nsim = 4)
chkid(identical(s1, s2, ignore.environment = TRUE), FALSE)

set.seed(123)
s1 <- simulate(fit, newdata = model.frame(fit)[1:6, ])
set.seed(123)
s2 <- simulate(fit, newdata = model.frame(fit)[1:6, ])
chkid(s1, s2, ignore.environment = TRUE)

## -- random effects
chkid(names(s2[[1]]),
      c("g1|(Intercept):1", "g1|x:1", "g2|(Intercept):1", "g2|(Intercept):2"))
## TODO check if it simulates from the correct RE distributions

summarize_tests()

options(oldopt)
