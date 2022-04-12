## =======================================================================
## Example for using tramME for inidvidual participant data meta-analysis
## =======================================================================

## The dataset is a !!!_random subset_!!! of the 3CIA database
## of COPD patients
##
## Variable definitions
##  follow_up_months: length of follow-up in months
##  alive1_death2: indicator of the status at the end of follow-up
##                 (1: alive, 2: death)
##  age: age in years
##  fev1pp: FEV1 measurement
##  mmrc0_4: dyspnea score (mMRC)
##  ages, fev1pps, mmrcs: rescaled versions of the three covariates
##                        (using their ranges)
##  cohort: study indicator
##  su: Surv object treating survival times as right-censored
##  sui: Surv object taking interval censoring into account

## Reference for the data:
## J. B. Soriano, B. Lamprecht, A. S. Ramirez, P. Martinez-Camblor, B. Kaiser,
## I. Alfageme, P. Almagro, C. Casanova, C. Esteban, J. J. Soler-Cataluna,
## J. P. de Torres, M. Miravitlles, B. R. Celli, J. M. Marin, M. A. Puhan,
## P. Sobradillo, P. Lange, A. L. Sternberg, J. Garcia-Aymerich,
## A. M. Turner, M. K. Han, A. Langhammer, L. Leivseth, P. Bakke,
## A. Johannessen, N. Roche, and D. D. Sin. -- Mortality prediction in chronic
## obstructive pulmonary disease comparing the GOLD 2007 and 2011 staging
## systems: A pooled analysis of individual patient data.
## The Lancet Respiratory Medicine, 3(6):443--450, 2015.
## <doi:10.1016/S2213-2600(15)00157-5>

library("tram")
library("tramME")
library("survival")

load(system.file(file.path("demo-data", "ipd.rda"), package = "tramME"))

sup <- c(1, 150) ## setting the support for basis approximation
opt <- optim_control(iter.max = 1e4, eval.max = 1e4, rel.tol = 1e-8,
                     method = "nlminb") ## optimization controls

## ====== MODEL 1
## - random proportional frailty term for heterogeneity in baseline risk
## - time independent prognostic factor effects
## - random slopes for heterogeneity of prognostic factor effects
fit_m1 <- CoxphME(sui ~ ages + mmrcs + fev1pps + (ages + mmrcs + fev1pps | cohort),
                  data = data3CIA, log_first = TRUE, order = 4, support = sup,
                  control = opt)
summary(fit_m1)


## ====== MODEL 2
## - stratification by study for heterogeneity in baseline risk
## - time independent prognostic factor effects
## - random slopes for heterogeneity of prognostic factor effects
fit_m2 <- CoxphME(sui | 0 + cohort ~ ages + mmrcs + fev1pps + (0 + ages + mmrcs + fev1pps | cohort),
                  data = data3CIA, log_first = TRUE, order = 4, support = sup,
                  control = opt)
summary(fit_m2)


## ====== MODEL 3
## - random proportional frailty term for heterogeneity in baseline risk
## - time dependent prognostic factor effects
## - random slopes for heterogeneity of prognostic factor effects
fit_m3 <- CoxphME(sui | ages + mmrcs + fev1pps ~ 1 + (ages + mmrcs + fev1pps | cohort),
                  data = data3CIA, log_first = TRUE, order = 4, support = sup,
                  control = opt)
summary(fit_m3)


## ====== MODEL 4
## - stratification by study for heterogeneity in baseline risk
## - time dependent prognostic factor effects
## - random slopes for heterogeneity of prognostic factor effects
fit_m4 <- CoxphME(sui | 0 + cohort + ages + mmrcs + fev1pps ~ 1 + (0 + ages + mmrcs + fev1pps | cohort),
                  data = data3CIA, log_first = TRUE, order = 4, support = sup, control = opt)
summary(fit_m4)


## ====== MODEL 5
## - stratification by study for heterogeneity in baseline risk
## - time independent prognostic factor effects
## - fixed effects only
fit_m5 <- Coxph(sui | 0 + cohort ~ ages + mmrcs + fev1pps, data = data3CIA,
                log_first = TRUE, order = 4, support = sup)
summary(fit_m5)


## ====== MODEL 6
## - stratification by study for heterogeneity in baseline risk
## - time dependent prognostic factor effects
## - fixed effects only
fit_m6 <- Coxph(sui | 0 + cohort + ages + mmrcs + fev1pps ~ 1, data = data3CIA,
                log_first = TRUE, order = 4, support = sup)
coef(fit_m6, with_baseline = TRUE)
logLik(fit_m6)
