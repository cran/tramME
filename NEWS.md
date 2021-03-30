# tramME 0.1.0 (2021-03-29)

* Updated internal functions for estimation transformation models with TMB:
    * Fixed effects only models can be estimated
    * Exported functions (`coef`, `logLik`) communicate with the TMB model more smoothly
    * Updated internal structure to help future development (not exported currently)
* Fixing coefficients with the argument `fixed = c(name = value)` 
* New model classes in `SurvregME` using the parameter fixing option
* `AaregME` extends the (parametric) Aalen regression model of `tram` with mixed-effects.
* Out-of-sample log-likelihoods using the `newdata` argument of `logLik`
* Score residuals with the `resid` method. When frequently recalculated (e.g. in boosting), setting `resid = TRUE` in model definition increases efficiency.
* Updating models via the `update` method
* Setting observation weights and offsets efficiently (activated  with `do_update = TRUE`). Currently, only possible through directly manipulating the `tmb_obj` of the model. 
* Parametric bootstrap with the `parboot` method
* New optimization options (e.g. internal scaling of fixed effects design matrix to improve convergence, more sensible initial values) and improved control over the optimization process with `optim_control()`
* Various methods (including `analytical` in the case of fixed effects only models) to calculate the Hessian; trying harder to invert the Hessian in numerically unstable cases. 
* Calculate the linear predictor with `type = "lp"` in `predict`
* Several additional methods to help the user working with `tramME` objects: `model.frame`, `model.matrix`, `fitmod`, `duplicate`.
* Improved unit testing
* Improved documentation
* Demo for IPD meta-analysis 

# tramME 0.0.4 (2021-02-04)

* fixed bug in setting error distributions of 'dummy' ctms for predict and simulate methods 
* updated Figure 6 in vignette, because the bug above affected predict in the case of CorlME
* fixed bug in unit test for simulate that caused error with mlt 1.2-1 
* fixed simulate output structure with `what = "joint"` option 

# tramME 0.0.3 (2020-07-30)

* Added author ORCID
* Fixed CRAN issue with unit test using MKL
* Figure 3 color/legend changed in Vignette

# tramME 0.0.2 (2020-03-30)

* Fixed numerical precision problem in unit tests

# tramME 0.0.1 (2020-03-20)

* First CRAN version

