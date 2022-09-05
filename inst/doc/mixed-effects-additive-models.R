## ----setup, echo=FALSE, message=FALSE, results="hide"-------------------------
knitr::opts_chunk$set(size = "small", prompt = TRUE, comment = NA,
                      out.width=".9\\linewidth")
knitr::knit_hooks$set(
  document = function(x) {sub('\\usepackage[]{color}', '\\usepackage{xcolor}',
                              x, fixed = TRUE)}
)
oldpar <- par(no.readonly = TRUE) ## NOTE: for setting back at the end
oldopt <- options()
options(prompt = "R> ", continue = "+  ")
options(width = 80, digits = 3)

## Dependencies
library("tramME")
library("survival")
library("mgcv")
library("glmmTMB")
library("xtable")
library("gamm4")

## ----tramME-exmpl1, eval=FALSE, echo=TRUE-------------------------------------
#  ## tramME is available from CRAN:
#  ## install.packages("tramME")
#  library("tramME")
#  mc1 <- CoxphME(                            ### conditional proportional hazards
#      Time ~                                 ### response time intervals
#             Insects +                       ### fixed effects
#             Habitat +
#             Landscape +
#             s(Temperature, k = 20) +        ### non-linear terms, as in mgcv
#             s(Elevation100, k = 20) +
#             (1 | PlotID),                   ### random intercept, as in lme4
#      data = carrion,                        ### data
#      log_first = TRUE,                      ### log(time) before modeling
#      order = 6                              ### order of Bernstein
#  )

## ----load-carrion, include=FALSE----------------------------------------------
carrion <- read.csv("carrion.csv")
carrion$Time <- with(carrion, Surv(time1, time2, type = "interval2"))
carrion$Time_rc <- with(carrion,
  Surv(ifelse(is.finite(time2), time2, time1), event = is.finite(time2)))
carrion$Insects <- factor(carrion$Insects, labels = c("no", "yes"))
carrion$Habitat <- factor(carrion$Habitat)
carrion$Habitat <- relevel(carrion$Habitat, ref = "forest")
carrion$Landscape <- factor(carrion$Landscape)
carrion$Landscape <- relevel(carrion$Landscape, ref = "seminatural")
carrion$PlotID <- factor(carrion$PlotID)

## ----fig-carrion-data, echo=FALSE, fig.width=4, fig.height=3.5, out.width=".5\\textwidth"----
sv <- survfit(Time ~ Insects, data = carrion)
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
plot(sv, lty = c(1, 2), lwd = 2, xlab = "", ylab = "P(T>t)")
title(xlab = "Follow-up time (days)", line = 2.3)
grid()
legend("topright", c("Insect = no", "Insect = yes"), lty = c(1, 2),
       lwd = 2, bty = "n")

## ----carrion-model------------------------------------------------------------
dcmp <- CoxphME(Time ~ Insects + Habitat + Landscape
                + s(Temperature, k = 20) + s(Elevation100, k = 20)
                + (1 | PlotID), data = carrion,
                log_first = TRUE, order = 6)
summary(dcmp)

## ----eval=FALSE---------------------------------------------------------------
#  plot(smooth_terms(dcmp))

## ----plot-carrion-smooth, echo=FALSE, fig.width=8, fig.height=4.5-------------
plot(smooth_terms(dcmp))

## ----resid, echo=FALSE--------------------------------------------------------
## From Surv object to a matrix with left and right censoring values
get_bounds <- function(x) {
  stopifnot(is.Surv(x))
  lb <- x[, 1]
  ub <- x[, 2]
  ty <- x[, 3]
  ub[ty == 2] <- lb[ty == 2]
  lb[ty == 2] <- -Inf
  ub[ty == 0] <- Inf
  ub[ty == 1] <- lb[ty == 1]
  cbind(lb, ub)
}

## Evaluate the marginal survivor function
marginalize <- function(model, data, type = "survivor",
  lower = -Inf, upper = Inf) {
  ## fun(y|x,g) * density(g)
  joint <- function(re, nd, mod, type) {
    nd <- nd[rep(1, length(re)), ]
    nd[[variable.names(mod, "grouping")]] <- seq(nrow(nd)) ## to take vector-valued REs
    pr <- predict(mod, newdata = nd, ranef = re, type = type) *
      dnorm(re, 0, sd = sqrt(varcov(mod)[[1]][1, 1]))
    c(pr)
  }
  ## integrate out the random effects row by row
  res <- parallel::mclapply(split(data, seq(nrow(data))), function(nd) {
    out <- integrate(joint, lower = lower, upper = upper, nd = nd, mod = model,
                     type = type)
    if (out$message == "OK") return(out$value)
    return(NA)
  }, mc.cores = 8)
  unlist(res)
}

## Cox-Snell residuals
## type:
##   ic: Evaluate the cumulative hazard at the left and right censoring values
##   and return residuals as interval-censored.
##   adj: Adjusted CS resid (Farrington,2000 <doi:10.1111/j.0006-341X.2000.00473.x>),
##   replace the intervals with the expected values under unit exponential on these
##   intervals (under the assumption of correct model fit)
CSresid <- function(model, data = model.frame(model), type = c("both", "adj", "ic")) {
  type <- match.arg(type)
  rv <- variable.names(model, "response")
  if (!is.Surv(data[[rv]])) data[[rv]] <- Surv(data[[rv]])
  bns <- get_bounds(data[[rv]])
  stopifnot(all(is.finite(bns[,1])))
  data[[rv]] <- bns[, 1]
  Sl <- marginalize(model, data, type = "survivor")
  Su <- numeric(length = nrow(bns))
  ie <- bns[,1] == bns[,2]
  Su[ie] <- Sl[ie]
  Su[ir <- is.infinite(bns[,2])] <- 0
  ii <- !(ie | ir)
  data[[rv]] <- bns[, 2]
  Su[ii] <- marginalize(model, data[ii, ], type = "survivor")
  res_ic <- Surv(-log(Sl), -log(Su), type = "interval2")
  if (type == "ic") return(res_ic)
  res <- numeric(length = nrow(bns))
  res[ie] <- -log(Sl[ie])
  res[ir] <- 1 - log(Sl[ir])
  res[ii] <- (Sl[ii] * (1 - log(Sl[ii])) - Su[ii] * (1 - log(Su[ii]))) /
    (Sl[ii] - Su[ii])
  if (type == "adj") return(res)
  data.frame(IC = res_ic, adjusted = res)
}

if (!file.exists("CS-resid.rda")) {
  resids <- CSresid(dcmp, type = "both")
  save(resids, file = "CS-resid.rda")
} else {
  load("CS-resid.rda")
}

## ----resid-plot, echo=FALSE, fig.width=7, fig.height=3.5----------------------
rf_cens <- survfit(IC ~ 1, data = resids) ## Turnbull NPMLE for interval-censored
rf_adj <- survfit(Surv(adjusted) ~ 1, data = resids) ## KM for adjusted

par(mfrow = c(1, 2), cex = 0.9, mar = c(4, 4, 2, 1))
plot(rf_cens$time, -log(rf_cens$surv), type = "s", lwd = 1, las = 1,
     xlab = "Cox-Snell residuals", ylab = "Cumulative hazard")
abline(0, 1, lwd = 2, lty = 2)
grid()
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "A", xpd = TRUE, cex = 1.2)

plot(rf_adj$time, rf_adj$cumhaz, type = "s", lwd = 1, las = 1,
     xlab = "Cox-Snell residuals", ylab = "Cumulative hazard")
abline(0, 1, lwd = 2, lty = 2)
grid()
blx <- grconvertX(0.10, from = "nfc", to = "user")
bly <- grconvertY(0.98, from = "nfc", to = "user")
text(blx, bly, labels = "B", xpd = TRUE, cex = 1.2)

## ----carrion-model2-----------------------------------------------------------
dcmp2 <- CoxphME(Time | Insects ~ Habitat + Landscape
                 + s(Temperature, k = 20) + s(Elevation100, k = 20)
                 + (1 | PlotID), data = carrion,
                 log_first = TRUE, order = 6)

## ----fig-carrion-tve, echo=FALSE, fig.width=4, fig.height=3.5, warning=FALSE, out.width=".5\\textwidth"----
cf <- coef(dcmp2, with_baseline = TRUE)
vc <- vcov(dcmp2, pargroup = "baseline")
cf <- cf[grep("Insectsyes", names(cf), fixed = TRUE)]
idx <- grep("Insectsyes", colnames(vc), fixed = TRUE)
vc <- vc[idx, idx]

ns <- 200
nd <- model.frame(dcmp2)[rep(1, ns), ]
nd[[variable.names(dcmp2, "response")]] <- seq(1, 100, length.out = ns)
X <- model.matrix(dcmp2, data = nd, type = "Y", simplify = TRUE)$Ye
idx <- grep("Insectsyes", colnames(X), fixed = TRUE)
X <- X[, idx]
ci <- confint(multcomp::glht(multcomp::parm(cf, vc), linfct = X),
              calpha = multcomp::univariate_calpha(), level = 0.95)$confint
## use multcomp::adjusted_calpha() for multiplicity adjustments

par(cex = 0.8, mar = c(4, 4, 2, 1), las = 1)
plot(nd$Time, ci[, 1], type = "l", col = 1, lwd = 2, ylim = c(-1, 3),
     panel.first = {grid(); abline(h = 0, lwd = 2, col = "lightgrey")},
     ylab = "Log-hazard ratio", xlab = "")
title(xlab = "Time (days)", line = 2.3)
polygon(c(nd$Time, rev(nd$Time)), c(ci[, 2], rev(ci[, 3])), border = NA,
        col = grey(0.5, 0.2))
ci2 <- confint(dcmp, parm = "Insectsyes", estimate = TRUE)
ci2 <- ci2[rep(1, ns), ]
matlines(nd$Time, ci2, col = 1, lwd = c(1, 1, 2), lty = c(3, 3, 2))
legend("bottomright", c("Time-varying effect", "Proportional effect"),
       lty = 1:2, lwd = 2, bty = "n", cex = 0.9)

## ----ecoli-data, echo=FALSE---------------------------------------------------
ecoli <- read.csv("Hulvey2021.csv")
fs <- c("treatment", "stream", "pasture", "cattle", "rotation")
ecoli[fs] <- lapply(ecoli[fs], factor)

## ----ecoli-est, echo=TRUE-----------------------------------------------------
## specifications w/o random effects
mf <- c(log10(ecoli_MPN) ~ treatment + cattle +
          s(DOY, bs = 'cr', by = treatment),
        log10(ecoli_MPN) ~ treatment + cattle + s(DOY, bs = 'cr'),
        log10(ecoli_MPN) ~ treatment + s(DOY, bs = 'cr', by = treatment),
        log10(ecoli_MPN) ~ cattle + s(DOY, bs = 'cr'),
        log10(ecoli_MPN) ~ treatment + s(DOY, bs = 'cr'),
        log10(ecoli_MPN) ~ s(DOY, bs = 'cr'))
names(mf) <- paste("Model", c(1:5, "Null"))
ecoli_res <- data.frame(matrix(NA, nrow = length(mf), ncol = 3))
colnames(ecoli_res) <- c("gamm", "LmME", "BoxCoxME")
rownames(ecoli_res) <- names(mf)
for (i in seq_along(mf)) {
  m_gamm <- gamm4(mf[[i]], data = ecoli,
                  random = ~ (1 | year:stream:pasture) + (1 | stream),
                  REML = FALSE)
  ecoli_res$gamm[i] <- logLik(m_gamm$mer)
  mf2 <- update(mf[[i]], . ~ . + (1 | year:stream:pasture) + (1 | stream))
  m_LmME <- LmME(mf2, data = ecoli)
  if (m_LmME$opt$convergence == 0) ecoli_res$LmME[i] <- logLik(m_LmME)
  m_BCME <- BoxCoxME(mf2, data = ecoli)
  if (m_BCME$opt$convergence == 0) ecoli_res$BoxCoxME[i] <- logLik(m_BCME)
}

## ----ecoli-tbl, echo=FALSE, results="asis"------------------------------------
names(ecoli_res) <- c("\\multicolumn{1}{m{2cm}}{\\centering GAMM}",
 "\\multicolumn{1}{m{4cm}}{\\centering Additive normal transformation model}",
 "\\multicolumn{1}{m{4cm}}{\\centering Additive non-normal transformation model}")
print(xtable(ecoli_res, align = c("l", "R{2cm}", rep("R{4cm}", 2))),
      floating = FALSE, booktabs = TRUE,
      sanitize.colnames.function = function(x){x})

## ----ecoli-show-m1------------------------------------------------------------
update(mf[[1]], . ~ . + (1 | year:stream:pasture) + (1 | stream))

## ----ecoli-plot-fun, echo=FALSE-----------------------------------------------
plot2cis <- function(x, y, col = c(1, 2),
                     fill = c(grey(0.1, 0.25), rgb(1, 0, 0, 0.25)),
                     xlabs = NULL, ylabs = NULL, mains = NULL,
                     ...) {
  stopifnot(length(x) == length(y))
  for (ii in seq_along(x)) {
    plot(0, type = "n",
         xlab = if (is.null(xlabs)) colnames(x[[ii]])[1] else xlabs[ii],
         ylab = if (is.null(ylabs)) colnames(x[[ii]])[2] else ylabs[ii],
         main = if (is.null(mains)) NULL else mains[ii],
         xlim = range(x[[ii]][, 1], y[[ii]][, 1]),
         ylim = range(x[[ii]][, 2:4], y[[ii]][, 2:4]),
         panel.first = grid(), ...)
    lines(x[[ii]][, 1], x[[ii]][, 2], col = col[1])
    lines(y[[ii]][, 1], y[[ii]][, 2], col = col[2])
    polygon(c(x[[ii]][, 1], rev(x[[ii]][, 1])),
            c(x[[ii]][, 3], rev(x[[ii]][, 4])),
            border = NA, col = fill[1])
    polygon(c(y[[ii]][, 1], rev(y[[ii]][, 1])),
            c(y[[ii]][, 3], rev(y[[ii]][, 4])),
            border = NA, col = fill[2])
  }
}

## NOTE: currently only w/ pnorm inverse link
smooth2ci <- function(sm, PI = FALSE, ilink = "pnorm") {
  PIfun <- switch(ilink, pnorm = function(x) pnorm(x / sqrt(2)),
                  stop("No other inverse links are available atm."))
  lapply(sm, function(x) {
    out <- matrix(0, nrow = nrow(x), ncol = 4)
    colnames(out) <- c(colnames(x)[c(1, ncol(x)-1)], "lwr", "upr")
    out[, 1] <- x[, 1]
    se <- rev(x)[, 1]
    yy <- rev(x)[, 2]
    out[, 2] <- yy
    out[, 3:4] <- yy + qnorm(0.975) * se %o% c(-1, 1)
    if (PI) out[, 2:4] <- PIfun(out[, 2:4])
    out
  })
}

## NOTE: this function assumes the smooth term: s(DOY, by = "treatment")
diff_smooth <- function(mod, DOY = NULL, n = 100) {
  XZ_sm <- function(XZ) {
    X <- XZ$X
    X[, attr(XZ$X, "type") != "sm"] <- 0
    Z_ <- Matrix::t(XZ$Zt)[, attr(XZ$Zt, "type") == "sm", drop = FALSE]
    Z <- tramME:::nullTMatrix(nrow = nrow(Z_), ncol = length(mod$param$gamma))
    Z[, attr(mod$param$gamma, "type") == "sm"] <- Z_
    list(X = X, Z = Z)
  }

  mf <- model.frame(mod, drop.unused.levels = TRUE)
  if (is.null(DOY)) DOY <- seq(range(mf$DOY)[1], range(mf$DOY)[2], length.out = n)
  nd_ <- expand.grid(DOY = DOY,
    treatment = factor(levels(mf$treatment), levels = levels(mf$treatment)))
  nd <- mf[rep(1, nrow(nd_)), ]
  nd[colnames(nd_)] <- nd_
  XZ <- model.matrix(mod, data = nd, type = c("X", "Zt"), keep_sign = FALSE)
  XZ1 <- XZ_sm(XZ)
  nd$DOY <- nd$DOY + 1
  XZ <- model.matrix(mod, data = nd, type = c("X", "Zt"), keep_sign = FALSE)
  XZ2 <- XZ_sm(XZ)
  XZ <- list(X = XZ2$X - XZ1$X, Z = tramME:::as_dgTMatrix(XZ2$Z - XZ1$Z))
  pr <- predict(mod$tmb_obj, newdata = XZ, scale = "lp")
  split(data.frame(DOY = nd_$DOY, df = pr$pred, se = pr$se), nd$treatment)
}

## ----ecoli-m1, echo=FALSE-----------------------------------------------------
fm1 <- log10(ecoli_MPN) ~ treatment + cattle +
  s(DOY, bs = "cr", by = treatment) +
  (1 | year:stream:pasture) + (1 | stream)
ecoli_m1 <- LmME(fm1, data = ecoli)
ecoli_m1_bc <- BoxCoxME(fm1, data = ecoli)

## ----plot-ecoli-m1, echo=FALSE, fig.width=9, fig.height=3---------------------
par(mfrow = c(1, 3), cex = 0.8, mar = c(4, 4, 2, 1), las = 1)
plot2cis(smooth2ci(diff_smooth(ecoli_m1), PI = TRUE),
         smooth2ci(diff_smooth(ecoli_m1_bc), PI = TRUE),
         mains = paste("treatment =", levels(ecoli$treatment)),
         ylabs = rep("PI", 3))
legend("topright", c("normal", "non-normal"), col = c(1, 2),
       lty = 1, bty = "n", cex = 0.9)

## ----plot-ecoli-trafo, echo=FALSE, fig.width=5.5, fig.height=4, out.width="0.5\\linewidth"----
nd <- model.frame(ecoli_m1_bc)[rep(1, 100), ]
nd[[1]] <- do.call(seq,
  c(as.list(range(log10(ecoli$ecoli_MPN))), length.out = 100))

Y <- model.matrix(ecoli_m1_bc, data = nd, type = "Y")$Ye
b <- coef(ecoli_m1_bc, with_baseline = TRUE)[1:7]
vc <- vcov(ecoli_m1_bc, pargroup = "baseline")
ci <- confint(multcomp::glht(multcomp::parm(b, vc), linfct = Y),
              calpha = multcomp::univariate_calpha())$confint

par(mar = c(4, 4, 1, 1), cex = 0.8, las = 1)
plot(0, type = "n", xlim = range(nd[[1]]), ylim = range(ci),
     xlab = variable.names(ecoli_m1_bc, "response"), ylab = "h(y)",
     panel.first = grid(), xaxs = "i", yaxs = "i")
lines(nd[[1]], ci[, 1], lwd = 2)
polygon(c(nd[[1]], rev(nd[[1]])), c(ci[, 2], rev(ci[, 3])), col = grey(0.2, 0.2),
        border = FALSE)
b2 <- coef(ecoli_m1, with_baseline = TRUE)
abline(b2[1:2], lwd = 2, lty = 2)
legend("topleft", c("Non-normal", "Normal"), lwd = 2, lty = 1:2, bty = "n")

## ----ecoli-cens---------------------------------------------------------------
fm1c <- update(fm1, Surv(log10(ecoli_MPN), event = ecoli_MPN < 2419.6) ~ .)
ecoli_m1_cens <- BoxCoxME(fm1c, data = ecoli)
summary(ecoli_m1_cens)

## ----ecoli-fe-tbl, results="asis", echo=FALSE---------------------------------
cis <- lapply(list(ecoli_m1, ecoli_m1_bc, ecoli_m1_cens), function(x) {
  ci <- confint(x, pargroup = "shift", estimate = TRUE)
  ci <- pnorm(ci[, c(3, 1, 2)] / sqrt(2))
  ci <- formatC(ci, format = "f", digits = 2)
  cbind(ci[, 1],
        apply(ci[, 2:3], 1, function(x) paste(x, collapse = "---")))
})

tbl <- do.call("cbind", cis)
rownames(tbl) <- c("treatment = medium", "treatment = short", "cattle = present")
add <- list()
add$pos <- list(0)
add$command <- paste(c("& \\multicolumn{2}{c}{Normal} &",
                       "\\multicolumn{2}{c}{Non-normal} &",
                       "\\multicolumn{2}{c}{Non-normal, censored} \\\\\n",
                       "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                       "\\cmidrule(lr){6-7}",
                       rep("& PI & 95\\% CI ", 3),
                       "\\\\\n"), collapse = " ")
print(xtable(tbl, align = "lrrrrrr"), include.colnames = FALSE,
      floating = FALSE, add.to.row = add,
      sanitize.text.function = function(x) x, booktabs = TRUE)

## ----plot-ecoli-cens, echo=FALSE, fig.width=9, fig.height=3-------------------
par(mfrow = c(1, 3), cex = 0.8, mar = c(4, 4, 2, 1), las = 1)
plot2cis(smooth2ci(diff_smooth(ecoli_m1), PI = TRUE),
         smooth2ci(diff_smooth(ecoli_m1_cens), PI = TRUE),
         mains = paste("treatment =", levels(ecoli$treatment)),
         ylabs = rep("PI", 3))
legend("topright", c("normal", "non-normal, censored"), col = c(1, 2),
       lty = 1, bty = "n", cex = 0.9)

## ----ecoli-notr---------------------------------------------------------------
f_nontr <- update(fm1, Surv(ecoli_MPN, event = ecoli_MPN < 2419.6) ~ .)
ecoli_nontr <- BoxCoxME(f_nontr, data = ecoli, log_first = TRUE)
summary(ecoli_nontr)

## ----algae-data, echo=FALSE---------------------------------------------------
andrew <- read.csv("andrew.csv",
                   colClasses = c(QUAD = "factor", PATCH = "factor"))
andrew$TREAT <- factor(andrew$TREAT, labels = c("control", "removal", "0.33", "0.66"))
andrew$TREAT <- factor(andrew$TREAT, levels = c("control", "0.33", "0.66", "removal"))
andrew$pALGAE <- andrew$ALGAE / 100
## summary(andrew)
ecdfs <- lapply(split(andrew, andrew$TREAT), function(x) ecdf(x$pALGAE))

## ----plot-algae-treatment, echo=FALSE, fig.width=5.5, fig.height=4.5, out.width=".6\\textwidth"----
x <- seq(0, 1, length.out = 100)
plot(x, ecdfs[[1]](x), type = "s", xlab = "Algae cover proportion", ylab = "ECDF",
     lty = 1, lwd = 2, xlim = c(0, 1), ylim = c(0, 1), las = 1, panel.first = grid())
for (ii in 2:4) {
  lines(x, ecdfs[[ii]](x), type = "s", lty = ii, lwd = 2)
}
legend("bottomright", levels(andrew$TREAT), lty = 1:4, lwd = 2, bty = "n", cex = 0.9)

## ----algae-glmm, eval=FALSE---------------------------------------------------
#  urchin_zib <- glmmTMB(pALGAE ~ TREAT + (1 | PATCH), ziformula = ~ TREAT,
#                        data = andrew, family = beta_family())

## ----algae-tram---------------------------------------------------------------
urchin_tram <- ColrME(
  Surv(pALGAE, pALGAE > 0, type = "left") ~ TREAT + (1 | PATCH),
  bounds = c(-0.1, 1), support = c(-0.1, 1), data = andrew,
  order = 6)
summary(urchin_tram)

## ----pointmass-plot, echo=FALSE, fig.width=7, fig.height=3.5, out.width=".8\\textwidth"----
par(mar = c(4, 4, 1, 1))
layout(mat = matrix(1:3, nrow = 1), widths = c(45, 10, 45))
nd <- model.frame(urchin_tram)[rep(21, 100), ]
nd[[1]] <- seq(-0.1, 1, length.out = 100)
ccdf_e <- predict(urchin_tram, newdata = nd, type = "distribution", ranef = "zero")
plot(nd[[1]], ccdf_e, type = "l", ylim = c(0, 1), xlim = c(-0.1, 1),
     panel.first = grid(), xlab = "y", ylab = "h(y)")
idx <- which(nd[[1]] < 0)
lines(nd[[1]][idx], ccdf_e[idx], col = 2, lwd = 3)
par(mar = c(1, 1, 1, 1))
plot(0:1, 0:1, type = "n", yaxt = "n", xaxt = "n", xlab = "", ylab = "", bty = "n")
arrows(x0 = 0, x1 = 1, y0 = 0.5, length = 0.1, lwd = 3)
nd[[1]] <- seq(0, 1, length.out = 100)
ccdf <- predict(urchin_tram, newdata = nd, type = "distribution", ranef = "zero")
par(mar = c(4, 4, 1, 1))
plot(nd[[1]], ccdf, type = "l", ylim = c(0, 1), xlim = c(-0.1, 1), lwd = 2,
     panel.first = grid(), xlab = "y", ylab = expression(F[Y] * "(y)"))
points(0, 0, pch = 21, cex = 1, lwd = 2)
points(0, ccdf[1], pch = 20, cex = 1, lwd = 2)
segments(x0 = -0.1, x1 = -0.01, y0 = 0, lwd = 2)

## ----algae-mCDF1, echo=FALSE--------------------------------------------------
## Numerical integration to get an estimate of the marginal distribution
marginalize.zib.glmmTMB <- function(model, data,
                                     type = c("distribution", "density"),
                                     lower = -Inf, upper = Inf) {
  type <- match.arg(type)
  ## for a single data point
  joint <- function(re, nd, mod, type) {
    mu <- plogis(predict(mod, newdata = nd, type = "link", re.form = NA) + re)
    mu <- ifelse(mu == 0, .Machine$double.eps, mu)
    mu <- ifelse(mu == 1, 1 - .Machine$double.eps, mu)
    sig <- predict(mod, newdata = nd, type = "disp")
    nu <- predict(mod, newdata = nd, type = "zprob")
    out <- sapply(mu, function(m) {
      switch(type,
             distribution = gamlss.dist::pBEZI(nd[[1]], mu = m, sigma = sig, nu = nu),
             density = gamlss.dist::dBEZI(nd[[1]], mu = m, sigma = sig, nu = nu),
             stop("Unknown function!"))
    })
    out * dnorm(re, mean = 0, sd = sqrt(VarCorr(mod)$cond[[1]][1]))
  }
  ## integrate out the random effects row by row
  res <- parallel::mclapply(split(data, seq(nrow(data))), function(nd) {
    out <- integrate(joint, lower = lower, upper = upper, nd = nd, mod = model,
                     type = type)
    if (out$message == "OK") return(out$value)
    return(NA)
  }, mc.cores = 8)
  unlist(res)
}

marginalize.tramME <- function(model, data, type = "survivor", add = 0,
  lower = -Inf, upper = Inf) {
  ## fun(y|x,g) * density(g)
  joint <- function(re, nd, mod, type) {
    nd <- nd[rep(1, length(re)), ]
    rv <- variable.names(mod, "response")
    nd[[rv]] <- nd[[rv]] + add
    nd[[variable.names(mod, "grouping")]] <- seq(nrow(nd)) ## to take vector-valued REs
    pr <- predict(mod, newdata = nd, ranef = re, type = type) *
      dnorm(re, 0, sd = sqrt(varcov(mod)[[1]][1, 1]))
    c(pr)
  }
  ## integrate out the random effects row by row
  res <- parallel::mclapply(split(data, seq(nrow(data))), function(nd) {
    out <- integrate(joint, lower = lower, upper = upper, nd = nd, mod = model,
                     type = type)
    if (out$message == "OK") return(out$value)
    return(NA)
  }, mc.cores = 8)
  unlist(res)
}


nd <- expand.grid(pALGAE = seq(0, 1, length.out = 100),
                  TREAT = factor(levels(andrew$TREAT),
                                 levels = c("control", "0.33", "0.66", "removal")),
                  PATCH = andrew$PATCH[1],
                  KEEP.OUT.ATTRS = FALSE)
colnames(nd)[1] <- variable.names(urchin_tram, "response")

if (!file.exists("urchin.rda")) {
  ## estimate the glmmTMB model
  urchin_zib <- glmmTMB(pALGAE ~ TREAT + (1 | PATCH), ziformula = ~ TREAT,
                        data = andrew, family = beta_family())
  ## save results in a single df
  urchin_res <- nd
  colnames(urchin_res)[1] <- "pALGAE"
  urchin_res$tramME_shift <- marginalize.tramME(urchin_tram,
                                                data = nd, type = "distribution")
  urchin_res$zib <- marginalize.zib.glmmTMB(urchin_zib, data = nd,
                                            type = "distribution")
  urchin_zib_ll <- data.frame(ll = logLik(urchin_zib), np = length(urchin_zib$obj$par))
} else{
  load("urchin.rda")
}

## ----plot-algae-mCDF1, echo=FALSE, fig.width=9, fig.height=3.5----------------
multiplot <- function(dfs, yn, xn = colnames(dfs[[1]])[1], ecdfs = NULL,
                      cols = 1:(length(yn)+!is.null(ecdfs)),
                      ltys = rep(1, length(yn)+!is.null(ecdfs)),
                      lwds = rep(1, length(yn)+!is.null(ecdfs)),
                      mains = NULL, ...) {
  if (!is.null(ecdfs)) stopifnot(length(dfs) == length(ecdfs))
  for (i in seq_along(dfs)) {
    col <- if (is.null(ecdfs)) cols else cols[-length(cols)]
    lty <- if (is.null(ecdfs)) ltys else ltys[-length(ltys)]
    lwd <- if (is.null(ecdfs)) ltys else lwds[-length(lwds)]
    main <- if (is.null(mains)) names(dfs)[i] else mains[i]
    matplot(dfs[[i]][[xn]], dfs[[i]][yn], type = "l",
            col = col, lty = lty, lwd = lwd, panel.first = grid(), ...,
            main = main)
    if (!is.null(ecdfs)) {
      lines(dfs[[i]][[xn]], ecdfs[[i]](dfs[[i]][[xn]]), type = "s",
            lty = rev(ltys)[1], col = rev(cols)[1], lwd = rev(lwds)[1],
            ...)
    }
  }
}

cols <- c("#E16A86", "#00AD9A")
layout(mat = matrix(c(1:4, rep(5, 4)), nrow = 2, byrow = TRUE),
       heights = c(7, 1))
par(mar = c(4, 4, 3, 1), las = 1)
multiplot(split(urchin_res, urchin_res$TREAT),
          yn = c("tramME_shift", "zib"),
          ecdfs = ecdfs,
          mains = paste("treatment =", levels(urchin_res$TREAT)),
          cols = c(cols, 1), lwds = c(rep(2, 2), 1),
          xlab = "pALGAE", ylab = "prob",
          xlim = c(0, 1))
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       c("Linear transformation model", "Zero-inflated beta", "ECDF"),
       col = c(cols, 1), lty = 1, bty = "n", lwd = c(rep(2, 2), 1), horiz = TRUE)

## ----algae-glmm2, eval=FALSE--------------------------------------------------
#  urchin_zib_disp <- glmmTMB(pALGAE ~ TREAT + (1 | PATCH),
#                             ziformula = ~ TREAT, dispformula = ~ TREAT,
#                             data = andrew, family = beta_family())

## ----algae-tram2--------------------------------------------------------------
urchin_tram_strat <- ColrME(
  Surv(pALGAE, pALGAE > 0, type = "left") | 0 + TREAT ~ 1 + (1 | PATCH),
  bounds = c(-0.1, 1), support = c(-0.1, 1), data = andrew,
  order = 6, control = optim_control(iter.max = 1e3, eval.max = 1e3,
                                     rel.tol = 1e-9))
summary(urchin_tram_strat)

## ----algae-mCDF2, echo=FALSE--------------------------------------------------
if (!("tramME_strat" %in% colnames(urchin_res))) {
  ## fit the glmm
  urchin_zib_disp <- glmmTMB(pALGAE ~ TREAT + (1 | PATCH),
                             ziformula = ~ TREAT, dispformula = ~ TREAT,
                             data = andrew, family = beta_family())
  ## marginalize
  urchin_res$tramME_strat <- marginalize.tramME(urchin_tram_strat,
    data = nd, type = "distribution")
  urchin_res$zib_disp <- marginalize.zib.glmmTMB(urchin_zib_disp,
    data = nd, type = "distribution")
}
if (nrow(urchin_zib_ll) == 1) {
  urchin_zib_ll <- rbind(urchin_zib_ll,
    c(logLik(urchin_zib_disp), length(urchin_zib_disp$obj$par)))
}
if (!file.exists("urchin.rda")) {
  save(urchin_res, urchin_zib_ll, file = "urchin.rda")
}

## ----plot-algae-mCDF2, echo=FALSE, fig.width=9, fig.height=3.5----------------
cols <- c("#E16A86", "#00AD9A")
layout(mat = matrix(c(1:4, rep(5, 4)), nrow = 2, byrow = TRUE),
       heights = c(7, 1))
par(mar = c(4, 4, 3, 1), las = 1)
multiplot(split(urchin_res, urchin_res$TREAT),
          yn = c("tramME_strat", "zib_disp"),
          ecdfs = ecdfs,
          mains = paste("treatment =", levels(urchin_res$TREAT)),
          cols = c(cols, 1), lwds = c(rep(2, 2), 1),
          xlab = "pALGAE", ylab = "prob",
          xlim = c(0, 1))
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       c("Stratified linear transformation model",
         "Zero-inflated beta with dispersion model",
         "ECDF"),
       col = c(cols, 1), lty = 1, bty = "n", lwd = c(rep(2, 2), 1), horiz = TRUE)

## ----algae-loglik, results="asis", echo=FALSE---------------------------------
lls <- data.frame(c(urchin_zib_ll$ll[1], urchin_zib_ll$ll[2],
                    logLik(urchin_tram), logLik(urchin_tram_strat)),
                  c(urchin_zib_ll$np[1], urchin_zib_ll$np[2],
                    length(urchin_tram$tmb_obj$par),
                    length(urchin_tram_strat$tmb_obj$par)))
rownames(lls) <- c("Zero-inflated beta without dispersion model",
                   "Zero-inflated beta with dispersion model",
                   "Linear transformation model",
                   "Stratified linear transformation model")
colnames(lls) <- c("$\\log\\mathcal{L}$", "Number of parameters")
print(xtable(lls, align = "lrr", digits = c(0, 2, 0)), sanitize.text.function = function(x) x,
      booktabs = TRUE, floating = FALSE)

## ----mosquito-data, echo=FALSE------------------------------------------------
AGO <- read.csv("Juarez2021.csv")
factors <- c("Community", "HouseID", "Year", "Income", "Placement")
AGO[factors] <- lapply(AGO[factors], factor)
AGO$AEAfemale <- as.integer(AGO$AEAfemale)

## ----cotram-model, echo=TRUE, eval=FALSE--------------------------------------
#  ## Additive count transformation model
#  ## See ?cotram::cotram for the documentation
#  CotramME <- function(formula, data,
#                       method = c("logit", "cloglog", "loglog", "probit"),
#                       log_first = TRUE, plus_one = log_first, prob = 0.9,
#                       ...) {
#    method <- match.arg(method)
#    rv <- all.vars(formula)[1]
#    stopifnot(is.integer(data[[rv]]), all(data[[rv]] >= 0))
#    data[[rv]] <- data[[rv]] + as.integer(plus_one)
#    sup <- c(-0.5 + log_first, quantile(data[[rv]], prob = prob))
#    bou <- c(-0.9 + log_first, Inf)
#    data[[rv]] <- as.Surv(R(data[[rv]], bounds = bou))
#    fc <- match.call()
#    fc[[1L]] <- switch(method, logit = quote(ColrME), cloglog = quote(CoxphME),
#                       loglog = quote(LehmannME), probit = quote(BoxCoxME))
#    fc$method <- NULL
#    fc$plus_one <- NULL
#    fc$prob <- NULL
#    fc$log_first <- log_first
#    fc$bounds <- bou
#    fc$support <- sup
#    fc$data <- data
#    out <- eval(fc, parent.frame())
#    out$call$data <- match.call()$data
#    class(out) <- c("CotramME", class(out))
#    out
#  }
#  mosquito_tram <- CotramME(AEAfemale ~ Year + Income*Placement
#    + s(Week) + s(CovRate200) + (1|HouseID)
#    + (1|Community), offset = -log(daystrapping), data = AGO,
#    method = "logit", order = 5, log_first = TRUE, prob = 0.9)

## ----mosquito-est, echo=FALSE-------------------------------------------------
if (file.exists("mosquito_models.rda")) {
  load("mosquito_models.rda")
} else {
  ## GAMMs
  mosquito_pois <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement
                         + s(Week) + s(CovRate200), random =~ (1|HouseID)
                         + (1|Community), data = AGO, family=poisson)

  mosquito_nb <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement
                       + s(Week) + s(CovRate200), random =~ (1|HouseID)
                       + (1|Community), data = AGO, family = negative.binomial(1))

  ## additive count transformation model
  CotramME <- function(formula, data,
                       method = c("logit", "cloglog", "loglog", "probit"),
                       log_first = TRUE, plus_one = log_first, prob = 0.9,
                       ...) {
    method <- match.arg(method)
    rv <- all.vars(formula)[1]
    stopifnot(is.integer(data[[rv]]), all(data[[rv]] >= 0))
    data[[rv]] <- data[[rv]] + as.integer(plus_one)
    sup <- c(-0.5 + log_first, quantile(data[[rv]], prob = prob))
    bou <- c(-0.9 + log_first, Inf)
    data[[rv]] <- as.Surv(R(data[[rv]], bounds = bou))
    fc <- match.call()
    fc[[1L]] <- switch(method, logit = quote(ColrME), cloglog = quote(CoxphME),
                       loglog = quote(LehmannME), probit = quote(BoxCoxME))
    fc$method <- NULL
    fc$plus_one <- NULL
    fc$prob <- NULL
    fc$log_first <- log_first
    fc$bounds <- bou
    fc$support <- sup
    fc$data <- data
    out <- eval(fc, parent.frame())
    out$call$data <- match.call()$data
    class(out) <- c("CotramME", class(out))
    out
  }

  mosquito_tram <- CotramME(AEAfemale ~ Year + Income*Placement
                            + s(Week) + s(CovRate200) + (1|HouseID)
                            + (1|Community), offset = -log(daystrapping), data = AGO,
                            method = "logit", order = 5, log_first = TRUE, prob = 0.9)
  stopifnot(mosquito_tram$opt$convergence == 0)
}

## ----mosquito-ll-tbl, echo=FALSE, results="asis"------------------------------
if (!file.exists("mosquito_models.rda")) {
  ll_mosquito <- c("Poisson GAMM" = logLik(mosquito_pois$mer),
                   "Negative binomial GAMM" = logLik(mosquito_nb$mer),
                   "Additive count transformation model" = logLik(mosquito_tram))
  ll_mosquito <- data.frame(ll_mosquito)
  names(ll_mosquito) <- "Log-likelihood"
}
print(xtable(ll_mosquito, align = "lr"), floating = FALSE, booktabs = TRUE)

## ----mosquito_summaries, echo=FALSE-------------------------------------------
if (!file.exists("mosquito_models.rda")) {
  sums_mosquito <- list(
    nb = summary(mosquito_nb$gam),
    tram = summary(mosquito_tram)
  )
}

## ----plot-mosquito-smooth, echo=FALSE, fig.width=9, fig.height=7--------------
if (!file.exists("mosquito_models.rda")) {
  nd <- AGO[rep(1, 100), ]
  nd$Week <- seq(min(AGO$Week), max(AGO$Week), length.out = 100)
  nd$CovRate200 <- seq(min(AGO$CovRate200), max(AGO$CovRate200), length.out = 100)
  pr <- predict(mosquito_nb$gam, type = "terms",
                terms = c("s(Week)", "s(CovRate200)"),
                newdata = nd, se.fit = TRUE)

  sm_mosquito_nb <- lapply(colnames(pr[[1]]), function(n) {
    ci <- pr$fit[, n] + qnorm(0.975) * pr$se.fit[, n] %o% c(-1, 1)
    out <- data.frame(cbind(pr$fit[, n], ci))
    colnames(out) <- c("fit", "lwr", "upr")
    out
  })
  names(sm_mosquito_nb) <- colnames(pr[[1]])
  sm_mosquito_nb[[1]]$x <- nd$Week
  sm_mosquito_nb[[2]]$x <- nd$CovRate200

  sm_mosquito_tram <- smooth_terms(mosquito_tram)
}

plotsm <- function(sm, xlab, ylab) {
  plot(sm$x, sm$fit, type = "l", xlim = range(sm$x),
       ylim = range(sm[, c("lwr", "upr")]),
       xlab = xlab, ylab = ylab, panel.first = grid())
  polygon(c(sm$x, rev(sm$x)), c(sm$lwr, rev(sm$upr)),
          border = NA, col = grey(0.5, 0.25))
}

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), cex = 0.8, las = 1)
plotsm(sm_mosquito_nb[[1]], "Week", "s(Week)")
plotsm(sm_mosquito_nb[[2]], "CovRate200", "s(CovRate200)")
mtext("Negative binomial model", side = 3, line = -2, outer = TRUE)
sm <- sm_mosquito_tram
for (i in seq_along(sm_mosquito_tram)) {
  sm[[i]][, 2] <- -sm[[i]][, 2]
  plot(sm[i], panel.first = grid())
}
mtext("Count transfromation model", side = 3, line = -23, outer = TRUE)

## ----mosquito-fe-tbl, echo=FALSE, results="asis"------------------------------
if (!file.exists("mosquito_models.rda")) {
  formatCI <- function(x, digits = 2) {
    fx <- formatC(x, format = "f", digits = digits)
    fx <- matrix(paste0("$",
                        ifelse(c(x) > 0, paste0("\\phantom{-}", fx), fx),
                        "$"), ncol = 2)
    apply(fx, 1, paste, collapse = " ---")
  }

  b_nb <- sums_mosquito$nb$p.table[-1, 1]
  se_nb <- sums_mosquito$nb$p.table[-1, 2]
  ci_nb <- b_nb + qnorm(0.975) * se_nb %o% c(-1, 1)
  ci_nb <- formatCI(ci_nb)
  b_tr <- -sums_mosquito$tram$coef[, 1]
  se_tr <- sums_mosquito$tram$coef[, 2]
  ci_tr <- b_tr + qnorm(0.975) * se_tr %o% c(-1, 1)
  ci_tr <- formatCI(ci_tr)
  ci_mosquito <- data.frame(paste0("$", formatC(b_nb, format = "f", digits = 2), "$"),
                            ci_nb,
                            paste0("$", formatC(b_tr, format = "f", digits = 2), "$"),
                            ci_tr)
  rownames(ci_mosquito) <- c("Year = 2018", "Income = middle", "Placement = out",
                             "Income = middle \\& Placement = out")

  save(sm_mosquito_nb, sm_mosquito_tram,
       ll_mosquito, ci_mosquito, sums_mosquito,
       file = "mosquito_models.rda")
}
add <- list()
add$pos <- list(0)
add$command <- paste(c("& \\multicolumn{2}{c}{Negative binomial} &",
                       "\\multicolumn{2}{c}{Count transformation} \\\\\n",
                       "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
                       rep("& $\\widehat\\beta$ & 95\\% CI ", 2),
                       "\\\\\n"), collapse = " ")
print(xtable(ci_mosquito, align = c("@{}l", rep("r", 3), "r@{}")),
      include.colnames = FALSE, floating = FALSE, add.to.row = add,
      sanitize.text.function = function(x) x, booktabs = TRUE)

## ----dgp-code-----------------------------------------------------------------
##' @param n Numeric vector, number of observations in each group
##' @param beta Numeric, effect size
##' @param sfun Function with a numeric argument, non-linear smooth shift term
##' @param sigma Numeric, SD of random intercepts
##' @param hinv Function with a numeric argument, inverse transformation function
##' @param link Function with a single numeric argument on [0, 1] link function
##' @param scale Numeric, optional scaling constant on the transformation scale
##' @param seed Seed for the random number generator
##' @param two_sets Logical; generate both an estimation and a test sample
gen_smpl <- function(n, beta, sfun, sigma, hinv, link, scale = 1,
                     seed = NULL, two_sets = TRUE) {
  ## -- setting up the seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  if (two_sets) n <- rep(n, 2)
  x1 <- runif(sum(n))
  x2 <- runif(sum(n))
  gr <- factor(rep(seq_along(n), n))
  re <- rep(rnorm(length(n), mean = 0, sd = sigma), n)
  lp <- x1 * beta + sfun(x2) + re
  y  <- hinv((link(runif(sum(n))) + lp) * scale)
  if (two_sets) {
    n <- sum(n) / 2
    out <- list(est  = data.frame(x1 = x1[1:n], x2 = x2[1:n], gr = gr[1:n],
                                  y = y[1:n]),
                test =  data.frame(x1 = tail(x1, n), x2 = tail(x2, n),
                                   gr = tail(gr, n), y = tail(y, n)))
  } else {
    out <- data.frame(x1 = x1, x2 = x2, gr = gr, y = y)
  }
  out
}

## ----sim-setup, echo=FALSE----------------------------------------------------
## ==== Simulation inputs
b   <- qnorm(0.7) * sqrt(2) ## PI = 0.7, SD of the parametric part is sqrt(1/12) * b = 0.2141 on the LP scale
sf  <- function(x) 2 * b * sqrt(1/12) * sin(x * pi) / sqrt(2 - 16 / pi^2) ## same SD as the parametric terms on the LP scale
sig <- sqrt(1/12) * b ## same RE SD on LP scale
hi1 <- function(x) x
hi2 <- function(x) qchisq(pnorm(x), df = 3)
lnk <- qnorm

## ==== Helper function to extract results
get_res <- function(res, what, which = NULL, simplify = TRUE) {
  out <- lapply(res, function(x) {
    if (!is.null(which)) x <- x[which]
    sapply(x, function(xx) {
      if (length(xx) == 1 && is.na(xx)) NA
      else {
        if (is.function(what)) return(what(xx))
        xx[[what]]
      }
    })
  })
  if (simplify) {
    out <- do.call("rbind", out)
    if (!is.function(what))
      colnames(out) <- paste0(names(res[[1]]), "_", what)
  }
  out
}

## ----intcens-dgp--------------------------------------------------------------
make_intcens <- function(x, length = 1) {
  Surv(floor(x / length) * length, ceiling(x / length) * length,
       type = "interval2")
}

## ----sim1, echo=FALSE, results="hide"-----------------------------------------
if (!file.exists("sim1.rda")) {
  seeds <- 301:800
  system.time({
    sims1 <- parallel::mclapply(seeds, function(s) {
      smpl <- gen_smpl(rep(4, 100), beta = b, sfun = sf, sigma = sig,
                       hinv = hi1, link = lnk, seed = s, two_sets = TRUE)
      m1 <- LmME(y ~ x1 + s(x2) + (1 | gr), data = smpl$est)
      if (m1$opt$convergence == 0) {
        lmm <- list(ll    = logLik(m1),
                    lloos = logLik(m1, newdata = smpl$test, type = "fix_smooth"),
                    beta  = coef(m1), sigma = sqrt(varcov(m1)$gr[1,1]),
                    theta = m1$param$beta[attr(m1$param$beta, "type") == "bl"],
                    beta_se = sqrt(vcov(m1, parm = "x1")[1,1]))
      } else {
        lmm <- NA
      }
      m2 <- BoxCoxME(y ~ x1 + s(x2) + (1 | gr), data = smpl$est)
      if (m2$opt$convergence == 0) {
        bcm <- list(ll    = logLik(m2),
                    lloos = logLik(m2, newdata = smpl$test, type = "fix_smooth"),
                    beta  = coef(m2), sigma = sqrt(varcov(m2)$gr[1,1]),
                    theta = m2$param$beta[attr(m2$param$beta, "type") == "bl"],
                    beta_se = sqrt(vcov(m2, parm = "x1")[1,1]),
                    intercept = mean(sf(smpl$est$x2)))
      } else {
        bcm <- NA
      }
      list(LmME = lmm, BoxCoxME = bcm)
    }, mc.cores = 6)
  })
  save(sims1, file = "sim1.rda")
} else {
  load("sim1.rda")
}

## -- non-convergence
rowSums(sapply(sims1, function(x)
  sapply(x, function(xx) c(length(xx) == 1 && is.na(xx))))
  )

## ----sim1-beta-plot, echo=FALSE, fig.width=4, fig.height=3.5, out.width=".5\\textwidth"----
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
boxplot(pnorm(get_res(sims1, "beta") / sqrt(2)) - 0.7,
        ylab = "PI", names = c("Normal", "Non-normal"),
        boxwex = 0.5)
grid(nx = NA, ny = NULL)
abline(h = 0, lwd = 2, lty = 2, col = 2)

## ----sim1-beta-ci-------------------------------------------------------------
cvi <- get_res(sims1, what = function(x) {
  ci <- x$beta + x$beta_se * qnorm(0.975) * c(-1, 1)
  (ci[1] <= b) && (b < ci[2])
})
(cvr <- colMeans(cvi, na.rm = TRUE)) ## Coverage rate
sqrt(cvr * (1 - cvr) / nrow(cvi)) ## Monte Carlo SE

## ----sim1-pred-plot, echo=FALSE, fig.width=4, fig.height=3.5, out.width=".5\\textwidth"----
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
boxplot(get_res(sims1, "lloos"), ylab = "Out-of-sample log-likelihood",
        names = c("Normal", "Non-normal"), boxwex = 0.5)
grid(nx = NA, ny = NULL)

## ----sim1-baseline-plot, echo=FALSE, fig.width=4, fig.height=3.5, out.width=".5\\textwidth"----
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
th <- get_res(sims1, "theta", "BoxCoxME", simplify = FALSE)
th <- do.call("cbind", th)
## NOTE: to compare with the true transformation function (identity) we need to
## adjust the fitted baseline transformation with the expectation of the smooth
## term
ic <- get_res(sims1, "intercept", "BoxCoxME", simplify = FALSE)
ic <- sapply(ic, `[`, 1)
dat <- gen_smpl(rep(4, 100), beta = b, sfun = sf, sigma = sig,
                hinv = hi1, link = lnk, seed = 1, two_sets = FALSE)
bcm <- BoxCoxME(y ~ x1 + s(x2) + (1 | gr), data = dat, nofit = TRUE)
nd  <- dat[rep(1, 100), ]
nd[[variable.names(bcm, "response")]] <- seq(-1, 2, length.out = 100)
mm  <- model.matrix(bcm$model$ctm$bases$response, data = nd)
tr <- mm %*% th + matrix(rep(ic, each = 100), nrow = 100)
matplot(nd$y, tr, type = "l", col = grey(0.1,0.1), lty = 1,
        ylab = "h(y)", xlab = "y", panel.first = grid())
abline(0,1, lwd = 2, col = 2)

## ----intcens-example----------------------------------------------------------
dat <- gen_smpl(rep(4, 100), beta = b, sfun = sf, sigma = sig,
                hinv = hi2, link = lnk, seed = 1, two_sets = FALSE)
head(dat$y) ## exact values
head(make_intcens(dat$y, length = 2)) ## interval-censored version

## ----sim2, echo=FALSE, results="hide"-----------------------------------------
nd <- data.frame(x2 = seq(0, 1, length.out = 100))
if (!file.exists("sim2.rda")) {
  seeds <- 1001:1500
  system.time({
    sims2 <- parallel::mclapply(seeds, function(s) {
      smpl <- gen_smpl(rep(4, 100), beta = b, sfun = sf, sigma = sig,
                       hinv = hi2, link = lnk, seed = s, two_sets = TRUE)
      m1 <- LmME(y ~ x1 + s(x2) + (1 | gr), data = smpl$est)
      if (m1$opt$convergence == 0) {
        lmm <- list(ll    = logLik(m1),
                    lloos = logLik(m1, newdata = smpl$test, type = "fix_smooth"),
                    beta  = coef(m1), sigma = sqrt(varcov(m1)$gr[1,1]),
                    beta_se = sqrt(vcov(m1, parm = "x1")[1,1]))
      } else {
        lmm <- NA
      }
      m2 <- BoxCoxME(y ~ x1 + s(x2) + (1 | gr), data = smpl$est)
      if (m2$opt$convergence == 0) {
        sm <- smooth_terms(m2, newdata = nd)[[1]][,2] + mean(sf(smpl$est$x2))
        bcm <- list(ll    = logLik(m2),
                    lloos = logLik(m2, newdata = smpl$test, type = "fix_smooth"),
                    beta  = coef(m2), sigma = sqrt(varcov(m2)$gr[1,1]),
                    beta_se = sqrt(vcov(m2, parm = "x1")[1,1]),
                    smooth = sm)
      } else {
        bcm <- NA
      }
      smpl$est$y <- make_intcens(smpl$est$y, length = 2)
      m3 <- BoxCoxME(y ~ x1 + s(x2) + (1 | gr), data = smpl$est)
      if (m3$opt$convergence == 0) {
        sm <- smooth_terms(m3, newdata = nd)[[1]][, 2] + mean(sf(smpl$est$x2))
        icm <- list(ll    = logLik(m3),
                    lloos = logLik(m3, newdata = smpl$test, type = "fix_smooth"),
                    beta  = coef(m3), sigma = sqrt(varcov(m3)$gr[1,1]),
                    beta_se = sqrt(vcov(m3, parm = "x1")[1,1]),
                    smooth = sm)
      } else {
        icm <- NA
      }
      list(LmME = lmm, BoxCoxME = bcm, BoxCoxME_ic = icm)
    }, mc.cores = 8)
  })
  save(sims2, file = "sim2.rda")
} else {
  load("sim2.rda")
}

## -- non-convergence
rowSums(sapply(sims2, function(x)
  sapply(x, function(xx) c(length(xx) == 1 && is.na(xx))))
  )

## ----sim2-beta-ci-------------------------------------------------------------
cvi <- get_res(sims2, what = function(x) {
  ci <- x$beta + x$beta_se * qnorm(0.975) * c(-1, 1)
  (ci[1] <= b) && (b < ci[2])
})
(cvr <- colMeans(cvi, na.rm = TRUE)) ## Coverage rate
sqrt(cvr * (1 - cvr) / nrow(cvi)) ## Monte Carlo SE

## ----sim2-beta-plot, echo=FALSE, fig.width=4.5, fig.height=3.5, out.width=".5\\textwidth"----
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
boxplot(pnorm(get_res(sims2, "beta") / sqrt(2)) - 0.7,
        ylab = "PI", names = c("Normal", "Non-normal", "Non-normal (IC)"),
        boxwex = 0.5)
grid(nx = NA, ny = NULL)
abline(h = 0, lwd = 2, lty = 2, col = 2)

## ----sim2-pred-plot, echo=FALSE, fig.width=4.5, fig.height=3.5, out.width=".5\\textwidth"----
par(cex = 0.8, mar = c(4, 4.2, 2, 1), las = 1)
boxplot(get_res(sims2, "lloos"), ylab = "Out-of-sample log-likelihood",
        names = c("Normal", "Non-normal", "Non-normal (IC)"),
        boxwex = 0.5)
grid(nx = NA, ny = NULL)

## ----sim2-smooth-plot, echo=FALSE, fig.width=6, fig.height=3.5, out.width=".8\\textwidth"----
par(mfrow = c(1, 2), las = 1, mar = c(4, 5, 1, 1), cex = 0.8)
sms <- get_res(sims2, what = "smooth", which = c("BoxCoxME", "BoxCoxME_ic"),
               simplify = FALSE)
sm1 <- sapply(sms, function(x) x[, 1])
matplot(nd$x2, sm1, type = "l", col = grey(0.1, 0.1), lty = 1,
        ylab = expression(hat(f) * "(x)"), xlab = "x",
        panel.first = grid())
lines(nd$x2, sf(nd$x2), col = 2, lwd = 2)
sm2 <- sapply(sms, function(x) x[, 2])
matplot(nd$x2, sm2, type = "l", col = grey(0.1, 0.1), lty = 1,
        ylab = expression(hat(f) * "(x)"), xlab = "x",
        panel.first = grid())
lines(nd$x2, sf(nd$x2), col = 2, lwd = 2)

## ----info---------------------------------------------------------------------
sessionInfo()

