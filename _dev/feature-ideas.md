# Feature and improvement ideas (2026-11-12)

Collected while auditing the package and writing the test suite. Ordered roughly by how much I think each would help users, within each section. Multivariate treatments are already in progress (`_dev/adrf_mv.R`) and are not repeated here.

## Inference

**Respecting the `vcov` setting of *WeightIt* fits.** `adrf()` always uses *WeightIt*'s M-estimation estimating functions, even for a `glm_weightit()` fit with `vcov = "HC0"`. Passing `asympt = (x$vcov_type == "asympt")` to `estfun()` would honor the user's choice.

**Barnard–Rubin degrees of freedom for multiply imputed data.** Currently the minimum complete-data df is used. With `df` per estimate, `summary()` would also give sensible t intervals for small samples.

**Bootstrapping with multiply imputed data.** Bootstrap within each imputation and pool, or allow `vcov = "fwb"` with the weights shared across imputations. `adrf.mira()` currently has no `fwb.args` argument.

**Parallel bootstrapping.** Documenting `fwb.args = list(cl = <n>)` would help with large models.

**One-sided and equivalence omnibus tests.** For example, a monotonicity test (H0: the AMEF is nonpositive everywhere) or a test that the ADRF stays within a margin of a reference value. The simultaneous-band machinery already in `summary.curve_est()` is most of what's needed.

**A leverage-adjusted unconditional variance.** The omnibus tests are anti-conservative in moderate samples (7.5% to 9.75% rejection at the .05 level with n = 500; see `_dev/omnibus-pvalue-study/README.md`), because the unconditional variance shares the small-sample bias of HC0 sandwich estimators. Scaling the coefficient influence functions by leverage, as HC2 or HC3 do (e.g., residuals divided by `1 - h_ii` for linear models), would bring the size closer to nominal; conditional HC3 rejected 4.0% in the same simulation. The study in `_dev/omnibus-size-study/` found HC2 accurate (within a percentage point of nominal from n = 250) and HC3 conservative, and, for models without hat values, the average of the current variance and a grouped jackknife variance close to HC2.

**Testing part of the curve.** A `range` argument to `summary.effect_curve()` so the omnibus test can be restricted to a subinterval without re-running `adrf()`.

## Estimands

**Shift and incremental interventions.** `E[Y(A + delta)]` (or `E[Y(A * (1 + delta))]`) averaged over the observed treatment distribution answers "what if everyone's exposure were reduced by delta", which is often the question asked with continuous exposures and needs less extrapolation than the ADRF. The averaging machinery in `.effect_curve_internal()` handles it directly by setting the counterfactual treatment per unit instead of to a constant.

**Ratio contrasts.** `reference_curve(type = "ratio")` and `curve_contrast(type = "ratio")` for risk ratios and rate ratios, with inference on the log scale.

**Thresholds.** The treatment value at which the ADRF crosses a given outcome value, or the minimum effective dose, with an interval by inverting the simultaneous band.

**Average derivative effect.** The AMEF averaged over the treatment distribution (a single number summarizing the curve's slope), with its variance from the existing contrast-matrix approach.

**Multi-category outcomes.** Ordinal and multinomial models (`polr`, `ordinal_weightit`, `multinom_weightit`, `nnet::multinom`, *brms* `categorical`) return one prediction per category, which `adrf()` cannot currently average. An ADRF per category, labeled like `by` groups, would make them work, and with a `type` suited to each model (`"probs"` for `polr`).

**Survival outcomes.** `coxph_weightit` and other survival models with `type = "survival"` at a chosen time, or restricted mean survival time.

**Doubly robust estimation.** The pseudo-outcome kernel-smoothing estimator of Kennedy et al. (2017) for continuous treatments, taking a treatment-density model and an outcome model. It would fit the existing object design: the result is an `effect_curve` with an influence-function variance.

## Data handling

**`subset`, `by`, and `wts` from the full data**, not just model predictors (audit issue 31).

**A target-population argument.** Compute the ADRF in a supplied `newdata` (for transporting to another population) instead of the estimation sample, with the conditional variance by default since the unconditional variance would need the target's sampling design.

**Memory for large samples.** The counterfactual grid has `N * n` rows; for N = 100,000 and the default `n = 51` the model matrix has 5.1 million rows. Processing the grid in chunks of treatment values (as the Bayesian path already does) would cap memory.

**Quantile-spaced grids.** An option to place grid points at treatment quantiles so they are densest where the data are, which also improves the trapezoidal weights used by the omnibus test.

## Output and plotting

**Treatment distribution on plots.** A rug or density strip beneath the curve, so readers can see where the curve is extrapolating. A positivity diagnostic (conditional treatment density at each grid value) would go further.

**`tidy()`, `glance()`, and `as.data.frame()` methods.** *generics* is already in Suggests. `confint()` methods for `curve_est` and `curve_projection` would also be natural.

**Plot data access.** `plot(..., plot = FALSE)` or an exported helper returning the data frame behind the plot, for users who build their own figures.

**Pointwise and simultaneous bands together.** An option to draw both, as is common in the functional data literature.

**Richer print output.** The sample size, number of grid points, and number of subgroups in `print.effect_curve()`.

## Package internals

**Expose the interpolation settings** (`degree`, `kernel`, `bw`) in `adrf()` or as an option, and cap the degree for small grids (audit issue 9).

**Remove unused helpers** (`simultaneous_conf_level()`, `.gradientRich()`, `...mget()`, `all_apply()`, `.ul()`, `len()`), or export `simultaneous_conf_level()` if it is meant for users.

**Plot snapshot tests.** *vdiffr* snapshots would catch visual regressions that the current `layer_data()` checks miss (labels, legends, colors).
