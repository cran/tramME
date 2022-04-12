### R code from vignette source 'RJ-2021-075.Rnw'

###################################################
### code chunk number 1: RJ-2021-075.Rnw:21-26 (eval = FALSE)
###################################################
## Y <- model.matrix(ecoli_m1_bc, data = nd, type = "Y")$Ye
## b <- coef(ecoli_m1_bc, with_baseline = TRUE)[1:7]
## vc <- vcov(ecoli_m1_bc, pargroup = "baseline")
## ci <- confint(multcomp::glht(multcomp::parm(b, vc), linfct = Y),
##               calpha = multcomp::univariate_calpha())$confint


