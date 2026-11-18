# Changelog

## *adrftools* (development version)

- New functions
  [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md)
  and [`ade()`](https://ngreifer.github.io/adrftools/reference/ade.md).
  [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md)
  estimates the mean outcome when every unit’s treatment is shifted from
  its observed value (by a fixed amount, by a fixed proportion, or by a
  user-supplied function), as a function of the size of the shift, and
  so requires less extrapolation than the ADRF when the treatment is
  strongly related to the covariates.
  [`ade()`](https://ngreifer.github.io/adrftools/reference/ade.md)
  computes the average derivative effect, the derivative of the shift
  curve at 0. Both support the same inference options as
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).
  See the new
  [`vignette("shift-effects")`](https://ngreifer.github.io/adrftools/articles/shift-effects.md)
  for details.

- [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  and
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  have a new argument, `integrate`, which controls how treatment values
  are weighted in the omnibus test statistic and in the projection. The
  default, `"uniform"`, weights all treatment values in the range of the
  effect curve equally, as before; `"empirical"` weights them by a
  kernel density estimate of the treatment, computed by
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  from all units used (within each subgroup, and weighted by `wts`).
  Empirical weighting gives less weight to sparse regions of the
  treatment, where the effect curve is estimated least precisely, and
  had higher power in simulations with an unevenly distributed
  treatment.

- With `svyglm` models, unconditional standard errors now account for
  the survey design, including clusters, strata, finite population
  corrections, calibration, and replicate weights. Previously, only the
  survey weights were used. `cluster` is now ignored with a warning for
  these models.

- `subset`, `by`, `cluster`, and `wts` supplied as a variable name can
  now refer to variables in the original dataset that are not in the
  model. Previously, only variables in the model could be used.

- Bootstrapping now refits the model with its own weights multiplied by
  the bootstrap weights, and, for models fit with a `weightit` object
  from *WeightIt*, re-estimates the weights in each replicate.
  Previously, the model was refit with the bootstrap weights multiplied
  by `wts`, which dropped any weights used to fit the model, and the
  refit failed when the model was fit inside a function or to data with
  missing values. Bootstrapping is no longer allowed with `svyglm`
  models.

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) now
  computes predictions from the model matrix and coefficients only when
  doing so reproduces the model’s own predictions, and otherwise uses
  [`marginaleffects::get_predict()`](https://rdrr.io/pkg/marginaleffects/man/get_predict.html).
  Previously, the choice was made by the model’s class, which gave
  incorrect estimates for [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) models with offsets and
  errors for *fixest* models and for
  [`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.html)
  models with `type = "link"`.

- In
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md),
  `method = "imhof"` and `method = "davies"` now use the saddlepoint
  approximation for p-values below .0001, where the numerical error of
  these methods can exceed the p-value itself. `method = "saddlepoint"`
  is now implemented in *adrftools* and no longer requires *survey*;
  previously, with infinite `df`, it silently used the Satterthwaite
  approximation. When *CompQuadForm* is not installed, the default
  method is now `"saddlepoint"` rather than `"sim"`.

- With
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
  models fit with `vcov = "HC0"`, unconditional standard errors in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) now
  treat the weights as fixed, consistent with the model’s own covariance
  matrix. Previously, the estimation of the weights was always accounted
  for.

- The unconditional variance in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md),
  the default, now applies the HC2 small-sample adjustment of MacKinnon
  and White (1985), dividing each unit’s influence on the model
  coefficients by \\\sqrt{1 - h_i}\\, where \\h_i\\ is its leverage.
  Previously, the variance had the downward bias of the HC0 sandwich
  estimator, which is largest at the extremes of the treatment and made
  the omnibus tests in
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  reject a true null hypothesis too often in small samples. The
  adjustment is not used with `svyglm` models; for models whose
  leverages cannot be found, the HC1 adjustment \\\frac{n}{n-p}\\ is
  used instead with a message.

- With `cluster`, the unconditional variance in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) now
  uses the CR2 adjustment described by Pustejovsky and Tipton (2018),
  the clustered analog of HC2, for
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) models, in place of the
  \\G/(G-1)\\ adjustment. Other models use the CR1 adjustment,
  \\\frac{G}{G-1}\frac{n-1}{n-p}\\, with a message.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which unconditional standard errors for models fit with
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) or
  [`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) were too small
  for families with an estimated dispersion (e.g.,
  [`gaussian()`](https://rdrr.io/r/stats/family.html)), because the
  bread *sandwich* computes for these models leaves out the dispersion
  that its estimating functions divide by. Models with extended
  families, such as
  [`mgcv::nb()`](https://rdrr.io/pkg/mgcv/man/negbin.html), now use the
  conditional variance with a warning.

- Error messages now use *arg*, which is a new dependency.

- *fwb* (\>= 1.0.0) is now required for bootstrapping.

- `wts` in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) is
  now checked for missing, non-finite, and negative values and for a
  positive sum in each subgroup.

- Multiway cluster-robust covariance matrices now have negative
  eigenvalues set to 0, as with `sandwich::vcovCL(fix = TRUE)`.
  Clustering levels with no units are now dropped, and a clustering
  variable with a single cluster gives an error. Previously, these gave
  missing or infinite variances.

- [`anova()`](https://rdrr.io/r/stats/anova.html) for `curve_projection`
  objects now gives an error when the projections use different
  transformations.

- Models with offsets now use the faster path in which predictions are
  computed from the model matrix and coefficients.

- Evaluating bootstrapped and Bayesian effect curves and computing
  contrasts and reference curves from them is now much faster.

- Estimates on a very small scale are now printed with enough digits to
  be read rather than as 0.

- With multiply imputed data, the degrees of freedom used for tests and
  confidence intervals now account for the between-imputation variance,
  using the method of Reiter (2007), or of Li, Raghunathan, and
  Rubin (1991) when the complete-data degrees of freedom are infinite.
  Previously, the smallest complete-data degrees of freedom across the
  imputed datasets were used, so inference for generalized linear models
  used the normal distribution regardless of the number of imputations.

- `vcov` can now be `"HC0"`, `"HC1"`, `"HC2"`, or `"HC3"` when `cluster`
  is supplied to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md), in
  which case it is passed as the type to
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html).
  Supplying `cluster` with any other string besides `"unconditional"`,
  `"conditional"`, `"boot"`, and `"fwb"`, or when `vcov` is a function
  or formula, now gives an error. Previously, these were ignored and the
  default cluster-robust covariance matrix was used.

- `subset` in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) now
  gives an informative error when it is a factor, contains character
  values that are not row names, or contains row numbers that are out of
  range or that mix positive and negative values. Previously, factors
  gave an unrelated error and the other values gave incorrect results or
  were ignored.

- `transform` in
  [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md),
  [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md),
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md),
  and
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  now accepts only transformations with a known inverse and derivative:
  `TRUE`, `FALSE`, a family object, a link object (e.g., from
  [`make.link()`](https://rdrr.io/r/stats/make.link.html)), or one of
  `log`, `sqrt`, `exp`, `qlogis`, and `qnorm`. Previously, any function
  could be supplied and was inverted and differentiated numerically,
  which gave incorrect confidence limits when more than a few estimates
  were computed or when an estimate was near the boundary of the
  transformation.

- [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md)
  applied to a reference curve now returns a plain AMEF, which is what
  it computed, since subtracting a reference value does not change the
  derivative. Previously, the result was labeled, plotted, and tested as
  a reference curve even though it was not 0 at the reference point. Use
  [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  on the output of
  [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md) to
  compute a reference AMEF.

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) now
  gives informative errors when `subset` selects no units, the two
  values in `range` are equal, or a variable in `by` has missing values.

- Evaluating an effect curve at treatment values outside the range of
  its grid points now gives an error, as does a value of `reference`
  outside that range in
  [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md).
  Values outside it by no more than 0.1% of its width (e.g., an end
  point supplied as it was printed, after rounding) are treated as the
  nearest end point. Previously, effect curves were extrapolated beyond
  their grid points, which gave unstable estimates and `NA` estimates
  just past the end points, and `reference` could be up to 1% of the
  width beyond them.

- Informative errors are now given when
  [`anova()`](https://rdrr.io/r/stats/anova.html) is used with
  `vcov = "none"`,
  [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md)
  is used with a single subgroup, a projection model has more
  coefficients than the number of treatment values at which the effect
  curve was computed or has linearly dependent terms, or the treatment
  is named `subset` or shares a name with an object used internally by
  the effect curve. `subset = TRUE` when evaluating or summarizing an
  effect curve now selects all subgroups.

- Simultaneous Wald p-values in
  [`summary()`](https://rdrr.io/r/base/summary.html) are now much
  faster, and simultaneous Wald critical values are faster in
  [`summary()`](https://rdrr.io/r/base/summary.html) and more accurate
  in [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`point_contrast()`](https://ngreifer.github.io/adrftools/reference/point_contrast.md).
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) now compute
  the critical value with the same settings, so their intervals and
  bands agree. Previously,
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`point_contrast()`](https://ngreifer.github.io/adrftools/reference/point_contrast.md)
  used less accurate settings than
  [`summary()`](https://rdrr.io/r/base/summary.html), and a simultaneous
  p-value was computed for every estimate, which took about a second
  each with many estimates.

- `method = "sim"` in
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  is now faster and uses much less memory.

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) is
  faster and uses much less memory. The model matrix is now built for
  one treatment value at a time rather than for all of them at once (on
  the [`?adrf`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  example with `n = 201`, peak memory fell from about 450 MB to about
  170 MB). Influence functions for unconditional standard errors for
  models fit with *WeightIt* are also no longer computed twice. Omnibus
  tests and projections of bootstrapped and Bayesian effect curves are
  faster when no transformation is used.

- Simultaneous p-values for Bayesian models are now computed by
  inverting the simultaneous intervals. Previously, they were smaller
  than the pointwise p-values. The calibrated tail probability of the
  simultaneous intervals (Montiel Olea and Plagborg-Møller, 2019) is now
  found exactly rather than by root finding, which also fixes an error
  that could occur with few posterior draws.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which a transformed version of the treatment (e.g., a
  [`poly()`](https://rdrr.io/r/stats/poly.html) basis) could be used in
  place of the treatment when the original dataset could not be
  recovered from the model, e.g., when the model was fit inside a
  function. Now, an error asks for the dataset to be supplied to `data`.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which an object with the same name as the dataset used to fit the
  model, found in the environment of the model’s formula, could be used
  in place of that dataset, even when the correct dataset was supplied
  to `data`. This gave incorrect estimates or unrelated errors when, for
  example, a formula created outside a function was used to fit a model
  inside it.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with Bayesian models and `by` in which estimates were assigned to the
  wrong subgroups and treatment values.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which unconditional standard errors were incorrect when the target
  weights were not all equal, including with `svyglm` models and
  *WeightIt* models fit with sampling weights.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with `svyglm` models in which unconditional standard errors were too
  small when `subset` was supplied or the treatment was the only
  predictor.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with `svyglm` models fit to replicate-weight designs in which the
  survey weights were not used to average the predictions.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with multiply imputed data in which `cluster` and values of `vcov`
  other than `"unconditional"` and `"none"` were ignored.

- Fixed several bugs in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with multiply imputed data, which now uses the same code as a single
  dataset with the estimates pooled afterward. Previously,
  `type = "link"`, `subset`, `vcov = "none"`, Bayesian models, and
  models other than [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) fits gave errors.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which negative indices supplied to `subset` kept the indexed units
  rather than dropping them.

- Fixed a bug in
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  with Bayesian models.

- Fixed a bug in
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  in which `model = "flat"` failed with bootstrapped or Bayesian effect
  curves.

- Fixed a bug in
  [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md)
  in which contrasts of bootstrapped or Bayesian effect curves failed.

- Fixed a bug in
  [`point_contrast()`](https://ngreifer.github.io/adrftools/reference/point_contrast.md)
  in which [`summary()`](https://rdrr.io/r/base/summary.html) failed
  when `by` was supplied in the original call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- Fixed a bug in
  [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md),
  [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md),
  and
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  in which a transformation that produced non-finite values for some
  estimates, bootstrap or posterior draws, or for the null value (e.g.,
  `qlogis` with values outside 0 and 1, or a null value of 0) gave an
  unrelated error, infinite test statistics, `NaN` p-values, or standard
  errors on a different scale from the estimates. Now, a null value
  outside the domain of the transformation and a transformation supplied
  to `transform` that cannot be applied give informative errors. When
  `transform = TRUE` and the link function of the outcome model cannot
  be applied to some bootstrap or posterior draws (e.g., draws between
  grid points just outside 0 and 1), the estimates are left
  untransformed with a warning.

- Fixed a bug in which interpolation between the points of an effect
  curve lost precision when the kernel weights were very unequal. This
  mostly affected
  [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md) at
  the ends of the treatment range, where derivatives could be off by
  about 0.1%.

- Fixed a bug in which simultaneous inference failed when the degrees of
  freedom were not a whole number. They are now rounded down for
  simultaneous inference.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which unconditional standard errors were too small for models fit with
  weights that included zeros. The default `wts` also failed for such
  models.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which `vcov = "conditional"` with `cluster` ignored the clusters when
  `cluster` was a formula with a transformation (e.g.,
  [`interaction()`](https://rdrr.io/r/base/interaction.html)), and gave
  an error when it was a data frame.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which `subset` supplied as row names selected the wrong units when the
  model dropped rows with missing values.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which a dataset supplied to `data` was not restricted to the rows used
  to fit the model.

- Fixed a bug in which estimates between the points of an effect curve,
  [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md),
  and the bandwidth used for interpolation depended on the scale of the
  treatment. With treatments on very small or very large scales,
  estimates were taken from the nearest grid point or failed.

- Fixed a bug in
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  and
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  in which standard errors for projections with `transform` omitted a
  term for the curvature of the transformation. Quadratic and cubic
  projections also no longer fail for treatments on a large scale.

- Fixed a bug in which p-values for estimates with no variance (e.g., a
  reference curve at its reference point) were 2 or 0 with Bayesian or
  bootstrapped curves, and simultaneous p-values for Bayesian curves
  failed. These p-values are now `NA`, as with Wald inference. Pointwise
  p-values for Bayesian curves are now consistent with the pointwise
  intervals.

- Fixed a bug in
  [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  in which summarizing a single estimate from an effect curve with a
  non-identity link failed.

- Fixed a bug in which [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) of a projection
  failed when `vcov = "none"` in the original call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- Fixed a bug in
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  in which `subset` used the wrong bootstrap or posterior draws.

- Fixed a bug in
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  in which p-values depended on the scale of the outcome and could be 1,
  `NaN`, or negative. Numerical tolerances are now relative to the
  variance of the estimates, and p-values are restricted to be between 0
  and 1.

- Fixed a bug in which estimates between grid points were inaccurate
  when `n` was 2 or 3 in the call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- Fixed a bug in which
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) failed for
  models whose family objects store only the name of the link function,
  such as those from *bartisan*.

- Fixed a bug in
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  and
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  in which a projection model that produced non-finite values gave an
  unrelated error instead of an informative one.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which bootstrapping failed when `cluster` was supplied. Supplying more
  than one clustering variable with bootstrapping now gives an
  informative error.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which `cluster` supplied as a formula failed when the model was fit to
  data with missing values or when the formula contained a
  transformation, such as
  [`interaction()`](https://rdrr.io/r/base/interaction.html).

- Fixed a bug in
  [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  and
  [`summary.curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  in which `subset` was silently ignored when `by` was not supplied in
  the original call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) or
  when the curve was a contrast. Now, an error is given, as with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- Fixed a bug in
  [`point_contrast()`](https://ngreifer.github.io/adrftools/reference/point_contrast.md)
  in which supplying the output of
  [`summary()`](https://rdrr.io/r/base/summary.html) gave contrasts
  printed as estimates of the original curve, with the critical value
  from the original summary.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which models with aliased (i.e., `NA`) coefficients failed with
  `vcov = "conditional"` and robust standard errors, and bootstrapping
  failed when a replicate aliased a coefficient, e.g., when all units in
  a rare factor level received a weight of 0.

- Fixed a bug in [`anova()`](https://rdrr.io/r/stats/anova.html) for
  `curve_projection` objects in which the null hypothesis was not
  printed when a model formula was long.
  [`anova()`](https://rdrr.io/r/stats/anova.html) also now gives an
  informative error when the extra terms cannot be tested because the
  projection has more terms than the effect curve can identify.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which models fit without a `data` argument (e.g., `glm(y ~ x)`) failed
  whenever the environment in which they were fit contained other
  objects.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which a character predictor with values absent from the units selected
  by `subset` gave an error.

- Fixed a bug in
  [`summary.curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  in which `subset` failed with bootstrapped effect curves, and
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) of the output failed or
  returned the coefficients of all subgroups.

- Fixed a bug in
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) in
  which unconditional standard errors were `NA` for models fit with
  `na.action = na.exclude`.

## *adrftools* 0.1.0

CRAN release: 2026-02-17

- Initial CRAN submission.
