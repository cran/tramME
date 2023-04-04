# tramME 1.0.4 (2023-04-04)

* Activate the multi-core capabilities of `TMB`
* `Resp()` function to extend the functionality of `Surv` objects
* `trafo` option for `plot.smooth.tramME` to plot monotonic transformations of
  the smooth terms
* Faster and more robust `vcov` and `edf_smooth` calculations.
* `simulate` method (with limited functionality in this version) returned 
* Updated additive `tramME` package vignette

# tramME 1.0.3 (2022-09-05)

* Adapted matrix type conversions to Matrix 1.5.0.
* Fixed typo in vignette `mixed-effects-additive-models`.

# tramME 1.0.2 (2022-07-04)

* Updated vignette `mixed-effects-additive-models`.
* Cleaned up the documentation. 

# tramME 1.0.1 (2022-04-23)

* Updated test utilities to avoid issues with MKL on CRAN. 

# tramME 1.0.0 (2022-04-12)

* Smooth shift terms as defined by `mgcv`. Currently only `s()` smooths are
  allowed.
* Cleaned up the source code and the documentation.
* The official reference is [Tamasi and Hothorn (2021)](https://doi.org/10.32614/RJ-2021-075), 
  vignette was also updated.
* New vignette with examples of mixed-effects transformation models with
  smooth shift terms.
* Some methods have been moved to development branch for proper testing and
  refactoring (`Aareg`, `lpterms`, `simulate`, `parboot`, `trafo`). These
  features may eventually find their way back to the package, probably in a
  slightly changed form.

# tramME 0.1.3 (2022-03-08)

* Fixed `duplicate.tramTMB` method to work with TMB 1.8.0.

# tramME 0.1.2 (2021-08-16)

* Added R-Forge URL

# tramME 0.1.1 (2021-04-05)

* Fixed uninitialized value problem indicated by valgrind

# tramME 0.1.0 (2021-03-29)

* Updated internal functions for estimation transformation models with TMB:
    * Fixed effects only models can be estimated
    * Exported functions (`coef`, `logLik`) communicate with the TMB model more
      smoothly
    * Updated internal structure to help future development (not exported
      currently)
* Fixing coefficients with the argument `fixed = c(name = value)` 
* New model classes in `SurvregME` using the parameter fixing option
* `AaregME` extends the (parametric) Aalen regression model of `tram` with
  mixed-effects.
* Out-of-sample log-likelihoods using the `newdata` argument of `logLik`
* Score residuals with the `resid` method. When frequently recalculated (e.g. in
  boosting), setting `resid = TRUE` in model definition increases efficiency.
* Updating models via the `update` method
* Setting observation weights and offsets efficiently (activated  with
  `do_update = TRUE`). Currently, only possible through directly manipulating
  the `tmb_obj` of the model. 
* Parametric bootstrap with the `parboot` method
* New optimization options (e.g. internal scaling of fixed effects design matrix
  to improve convergence, more sensible initial values) and improved control
  over the optimization process with `optim_control()`
* Various methods (including `analytical` in the case of fixed effects only
  models) to calculate the Hessian; trying harder to invert the Hessian in
  numerically unstable cases. 
* Calculate the linear predictor with `type = "lp"` in `predict`
* Several additional methods to help the user working with `tramME` objects:
  `model.frame`, `model.matrix`, `fitmod`, `duplicate`.
* Improved unit testing
* Improved documentation
* Demo for IPD meta-analysis 

# tramME 0.0.4 (2021-02-04)

* fixed bug in setting error distributions of 'dummy' ctms for predict and
  simulate methods 
* updated Figure 6 in vignette, because the bug above affected predict in the
  case of CorlME
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

