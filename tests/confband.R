## ===== confband =====
## -- Test utils & settings
source("test_util.R")
.run_test <- identical(Sys.getenv("NOT_CRAN"), "true")
oldopt <- options(digits = 4)
set.seed(100)

library("tramME")
library("survival")

m_cars <- BoxCoxME(dist ~ speed, data = cars)

## -- Consistency under various argument values
pr1 <- confband(m_cars, newdata = data.frame(speed = c(0, 1)))
pr2 <- confband(m_cars, newdata = NULL)
pr3 <- confband(m_cars, newdata = data.frame(speed = c(0, 1)),
                baseline_only = TRUE)

chkeq(pr1[[1]], pr2, check.attributes = FALSE)
chkeq(pr1[[1]], pr3[[2]], check.attributes = FALSE)

## set q
pr1 <- confband(m_cars, newdata = data.frame(speed = c(0, 1)),
                       q = c(20, 40))

## -- Consistency with predict
data("sleepstudy", package = "lme4")
m_sleep <- ColrME(Reaction ~ s(Days) + (Days | Subject), data = sleepstudy)

nd <- model.frame(m_sleep)[c(1, 171), ]
cb <- confband(m_sleep, newdata = nd, K = 100, type = "distribution")
re <- ranef(m_sleep)[[1]][c(1, 18), ]
pr <- predict(m_sleep, newdata = nd[, -1L], ranef = c(t(re)),
              K = 100, type = "distribution")
chkeq(cb[[1]][, 2], pr[, 1], check.attributes = FALSE)
chkeq(cb[[2]][, 2], pr[, 2], check.attributes = FALSE)

## set q
nd <- data.frame(dist = c(20, 40), speed = c(1, 1))
pr2 <- predict(m_cars, newdata = nd, type = "trafo")
chkeq(pr1[[2]][, 2], pr2, check.attributes = FALSE)

## TODO -- Adjust
## pr1 <- confband.tramME(m_cars, newdata = data.frame(speed = c(0, 1)),
##                        adjust = TRUE, K = 100, cheat = 20)

summarize_tests()

options(oldopt)
