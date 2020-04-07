library("tramME")

oldopt <- options(digits = 4)
chk <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE)))

## structure of the output
data("sleepstudy", package = "lme4")
fit <- LmME(Reaction ~ Days + (1 | Subject), data = sleepstudy)
nd <- sleepstudy
nd[["Reaction"]] <- NULL
pr <- predict(fit, newdata = nd, ranef = unlist(ranef(fit)), K = 100)
chk(dim(pr), c(100, nrow(nd)))
pr <- predict(fit)
chk(length(pr), nrow(sleepstudy))
pdf(file = NULL)
pl <- plot(fit)
dev.off()
chk(dim(pl), c(50, nrow(sleepstudy)))

## equality of predicted values
fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
## -- With predict.tramME method: for the first subject at two reaction values
nd <- sleepstudy[4, ]
nd[["Reaction"]] <- NULL
pr <- predict(fit, newdata = nd, ranef = c(0.2, -0.1), q = c(300, 320), type = "distribution")
## -- Manually (NOTE: negative=TRUE)
nd2 <- merge(data.frame(Reaction = c(300, 320)), nd, by = NULL)
cf <- coef(fit, with_baseline = TRUE)
br <- fit$model$response$basis
bs <- fit$model$fixef$bases$shifting ## NOTE: no interacting; already w/ negative = TRUE
re <- tramME:::.re_data(~ (Days | Subject), data = nd2, negative = TRUE)
prtr <- predict(br, newdata = nd2, coef = cf[-length(cf)]) + predict(bs, newdata = nd2, coef = cf[length(cf)]) +
  + Matrix::crossprod(re$Zt, c(0.2, -0.1))
pr2 <- pnorm(as.numeric(prtr))
chk(c(pr), pr2)

## check "zero" option of ranef
nd <- sleepstudy[c(4, 14, 24), ]
pr1 <- predict(fit, newdata = nd, ranef = rep(0, 6), q = c(300, 320), type = "distribution")
pr2 <- predict(fit, newdata = nd, ranef = "zero", q = c(300, 320), type = "distribution")
chk(pr1, pr2)

## Prediction with Surv objects as response
data("eortc", package = "coxme")
suppressPackageStartupMessages(library("survival"))
fit <- CoxphME(Surv(y, uncens) ~ trt + (1 | center/trt), data = eortc, log_first = TRUE)
pdf(file = NULL)
pl <- plot(fit, newdata = eortc[1, ], ranef = c(0, 0), type = "survivor", K = 100)
dev.off()
chk(dim(pl), c(100, 1))

## model w/o shift terms
fit <- BoxCoxME(Reaction ~ 1 + (1 | Subject), data = sleepstudy)
nd <- sleepstudy
nd$Reaction <- NULL
pr <- predict(fit, newdata = nd[1:2, ], ranef = 0, K = 10)
pr

options(oldopt)
