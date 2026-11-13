adrftools News and Updates
======

# *adrftools* (development version)

* Error messages now use *arg*, which is a new dependency.

* With `svyglm` models, unconditional standard errors now account for the survey design, including clusters, strata, finite population corrections, calibration, and replicate weights. Previously, only the survey weights were used. `cluster` is now ignored with a warning for these models.

* `subset`, `by`, `cluster`, and `wts` supplied as a variable name can now refer to variables in the original dataset that are not in the model. Previously, only variables in the model could be used.

* Bootstrapping now refits the model with its own weights multiplied by the bootstrap weights, and, for models fit with a `weightit` object from *WeightIt*, re-estimates the weights in each replicate. Previously, the model was refit with the bootstrap weights multiplied by `wts`, which dropped any weights used to fit the model, and the refit failed when the model was fit inside a function or to data with missing values. Bootstrapping is no longer allowed with `svyglm` models.

* `adrf()` now computes predictions from the model matrix and coefficients only when doing so reproduces the model's own predictions, and otherwise uses `marginaleffects::get_predict()`. Previously, the choice was made by the model's class, which gave incorrect estimates for `lm()` and `glm()` models with offsets and errors for *fixest* models and for `ordinal_weightit()` models with `type = "link"`.

* In `summary.effect_curve()`, `method = "imhof"` and `method = "davies"` now use the saddlepoint approximation for p-values below .0001, where the numerical error of these methods can exceed the p-value itself. `method = "saddlepoint"` is now implemented in *adrftools* and no longer requires *survey*; previously, with infinite `df`, it silently used the Satterthwaite approximation. When *CompQuadForm* is not installed, the default method is now `"saddlepoint"` rather than `"sim"`.

* With `glm_weightit()` models fit with `vcov = "HC0"`, unconditional standard errors in `adrf()` now treat the weights as fixed, consistent with the model's own covariance matrix. Previously, the estimation of the weights was always accounted for.

* *fwb* (>= 1.0.0) is now required for bootstrapping.

* `wts` in `adrf()` is now checked for missing, non-finite, and negative values and for a positive sum in each subgroup.

* Multiway cluster-robust covariance matrices now have negative eigenvalues set to 0, as with `sandwich::vcovCL(fix = TRUE)`. Clustering levels with no units are now dropped, and a clustering variable with a single cluster gives an error. Previously, these gave missing or infinite variances.

* `anova()` for `curve_projection` objects now gives an error when the projections use different transformations.

* Models with offsets now use the faster path in which predictions are computed from the model matrix and coefficients.

* Evaluating bootstrapped and Bayesian effect curves and computing contrasts and reference curves from them is now much faster.

* Estimates on a very small scale are now printed with enough digits to be read rather than as 0.

* With multiply imputed data, the degrees of freedom used for tests and confidence intervals now account for the between-imputation variance, using the method of Reiter (2007), or of Li, Raghunathan, and Rubin (1991) when the complete-data degrees of freedom are infinite. Previously, the smallest complete-data degrees of freedom across the imputed datasets were used, so inference for generalized linear models used the normal distribution regardless of the number of imputations.

* `vcov = "HC3"` can now be used with `cluster` in `adrf()` and is passed to `sandwich::vcovCL()`. Supplying `cluster` with any other string or with a function or formula supplied to `vcov` now gives an error. Previously, these were ignored and the default cluster-robust covariance matrix was used.

* `subset` in `adrf()` now gives an informative error when it is a factor, contains character values that are not row names, or contains row numbers that are out of range or that mix positive and negative values. Previously, factors gave an unrelated error and the other values gave incorrect results or were ignored.

* `transform` in `summary.curve_est()`, `plot.effect_curve()`, `summary.effect_curve()`, and `curve_projection()` now accepts only transformations with a known inverse and derivative: `TRUE`, `FALSE`, a family object, a link object (e.g., from `make.link()`), or one of `log`, `sqrt`, `exp`, `qlogis`, and `qnorm`. Previously, any function could be supplied and was inverted and differentiated numerically, which gave incorrect confidence limits when more than a few estimates were computed or when an estimate was near the boundary of the transformation.

* `amef()` applied to a reference curve now returns a plain AMEF, which is what it computed, since subtracting a reference value does not change the derivative. Previously, the result was labeled, plotted, and tested as a reference curve even though it was not 0 at the reference point. Use `reference_curve()` on the output of `amef()` to compute a reference AMEF.

* `adrf()` now gives informative errors when `subset` selects no units, the two values in `range` are equal, or a variable in `by` has missing values.

* Simultaneous p-values for Bayesian models are now computed by inverting the simultaneous intervals. Previously, they were smaller than the pointwise p-values. The calibrated tail probability of the simultaneous intervals (Montiel Olea and Plagborg-Møller, 2019) is now found exactly rather than by root finding, which also fixes an error that could occur with few posterior draws.

* Fixed a bug in `adrf()` in which a transformed version of the treatment (e.g., a `poly()` basis) could be used in place of the treatment when the original dataset could not be recovered from the model, e.g., when the model was fit inside a function. Now, an error asks for the dataset to be supplied to `data`.

* Fixed a bug in `adrf()` in which an object with the same name as the dataset used to fit the model, found in the environment of the model's formula, could be used in place of that dataset, even when the correct dataset was supplied to `data`. This gave incorrect estimates or unrelated errors when, for example, a formula created outside a function was used to fit a model inside it.

* Fixed a bug in `adrf()` with Bayesian models and `by` in which estimates were assigned to the wrong subgroups and treatment values.

* Fixed a bug in `adrf()` in which unconditional standard errors were incorrect when the target weights were not all equal, including with `svyglm` models and *WeightIt* models fit with sampling weights.

* Fixed a bug in `adrf()` with `svyglm` models in which unconditional standard errors were too small when `subset` was supplied or the treatment was the only predictor.

* Fixed a bug in `adrf()` with `svyglm` models fit to replicate-weight designs in which the survey weights were not used to average the predictions.

* Fixed a bug in `adrf()` with multiply imputed data in which `cluster` and values of `vcov` other than `"unconditional"` and `"none"` were ignored.

* Fixed several bugs in `adrf()` with multiply imputed data, which now uses the same code as a single dataset with the estimates pooled afterward. Previously, `type = "link"`, `subset`, `vcov = "none"`, Bayesian models, and models other than `lm()` and `glm()` fits gave errors.

* Fixed a bug in `adrf()` in which negative indices supplied to `subset` kept the indexed units rather than dropping them.

* Fixed a bug in `curve_projection()` with Bayesian models.

* Fixed a bug in `curve_projection()` in which `model = "flat"` failed with bootstrapped or Bayesian effect curves.

* Fixed a bug in `curve_contrast()` in which contrasts of bootstrapped or Bayesian effect curves failed.

* Fixed a bug in `point_contrast()` in which `summary()` failed when `by` was supplied in the original call to `adrf()`.

* Fixed a bug in `summary.curve_est()`, `plot.effect_curve()`, and `summary.effect_curve()` in which a transformation supplied to `transform` that produced non-finite values for some estimates (e.g., `qlogis` with estimates outside 0 and 1) gave an error or reported standard errors on a different scale from the estimates. Now, the estimates are left untransformed in this case.

* Fixed a bug in which interpolation between the points of an effect curve lost precision when the kernel weights were very unequal. This mostly affected `amef()` at the ends of the treatment range, where derivatives could be off by about 0.1%.

* Fixed a bug in which simultaneous inference failed when the degrees of freedom were not a whole number. They are now rounded down for simultaneous inference.

* Fixed a bug in `adrf()` in which unconditional standard errors were too small for models fit with weights that included zeros. The default `wts` also failed for such models.

* Fixed a bug in `adrf()` in which `vcov = "conditional"` with `cluster` ignored the clusters when `cluster` was a formula with a transformation (e.g., `interaction()`), and gave an error when it was a data frame.

* Fixed a bug in `adrf()` in which `subset` supplied as row names selected the wrong units when the model dropped rows with missing values.

* Fixed a bug in `adrf()` in which a dataset supplied to `data` was not restricted to the rows used to fit the model.

* Fixed a bug in which estimates between the points of an effect curve, `amef()`, and the bandwidth used for interpolation depended on the scale of the treatment. With treatments on very small or very large scales, estimates were taken from the nearest grid point or failed.

* Fixed a bug in `curve_projection()` and `summary.effect_curve()` in which standard errors for projections with `transform` omitted a term for the curvature of the transformation. Quadratic and cubic projections also no longer fail for treatments on a large scale.

* Fixed a bug in which p-values for estimates with no variance (e.g., a reference curve at its reference point) were 2 or 0 with Bayesian or bootstrapped curves, and simultaneous p-values for Bayesian curves failed. These p-values are now `NA`, as with Wald inference. Pointwise p-values for Bayesian curves are now consistent with the pointwise intervals.

* Fixed a bug in `summary.curve_est()` in which summarizing a single estimate from an effect curve with a non-identity link failed.

* Fixed a bug in which `vcov()` and `summary()` of a projection failed when `vcov = "none"` in the original call to `adrf()`.

* Fixed a bug in `summary.effect_curve()` in which `subset` used the wrong bootstrap or posterior draws.

* Fixed a bug in `summary.effect_curve()` in which p-values depended on the scale of the outcome and could be 1, `NaN`, or negative. Numerical tolerances are now relative to the variance of the estimates, and p-values are restricted to be between 0 and 1.

* Fixed a bug in which estimates between grid points were inaccurate when `n` was 2 or 3 in the call to `adrf()`.

* Fixed a bug in which `summary()` and `plot()` failed for models whose family objects store only the name of the link function, such as those from *bartisan*.

* Fixed a bug in `curve_projection()` and `summary.effect_curve()` in which a projection model that produced non-finite values gave an unrelated error instead of an informative one.

* Fixed a bug in `adrf()` in which bootstrapping failed when `cluster` was supplied. Supplying more than one clustering variable with bootstrapping now gives an informative error.

* Fixed a bug in `adrf()` in which `cluster` supplied as a formula failed when the model was fit to data with missing values or when the formula contained a transformation, such as `interaction()`.

* Fixed a bug in `summary.effect_curve()` and `summary.curve_projection()` in which `subset` was silently ignored when `by` was not supplied in the original call to `adrf()` or when the curve was a contrast. Now, an error is given, as with `plot()`.

* Fixed a bug in `point_contrast()` in which supplying the output of `summary()` gave contrasts printed as estimates of the original curve, with the critical value from the original summary.

* Fixed a bug in `adrf()` in which models with aliased (i.e., `NA`) coefficients failed with `vcov = "conditional"` and robust standard errors, and bootstrapping failed when a replicate aliased a coefficient, e.g., when all units in a rare factor level received a weight of 0.

* Fixed a bug in `anova()` for `curve_projection` objects in which the null hypothesis was not printed when a model formula was long. `anova()` also now gives an informative error when the extra terms cannot be tested because the projection has more terms than the effect curve can identify.

# *adrftools* 0.1.0

* Initial CRAN submission.
