# Feature and improvement ideas (updated 2026-11-24)

First collected on 2026-11-12 while auditing the package and writing the test suite, and checked against the development version on 2026-11-24. Ordered roughly by how much I think each would help users, within each section. Ideas added on 2026-11-24 are marked *(new)*. Multivariate treatments are in progress (`_dev/adrf_mv.R`), and decisions already on the task list (the default `integrate`, the default `hypothesis` for reference curves, the open questions about `shift_curve()`) are in `_dev/TASKS.md`; neither is repeated here. Ideas that have been implemented are listed at the end.

## Inference

**Bootstrapping with multiply imputed data.** Bootstrap within each imputation and pool, or allow `vcov = "fwb"` with the weights shared across imputations. `adrf.mira()` and `shift_curve.mira()` have no `fwb.args` argument, so `"boot"` and `"fwb"` are not available with `mira` objects.

**A small-sample reference distribution for the omnibus test** *(new)*. With HC2, the omnibus test still rejects 6.6% to 8.0% at the .05 level with a flexible logistic model (n up to 1,000) or with errors largest at the extremes of the treatment (n up to 500), because the variance estimate itself is variable (`_dev/omnibus-size-study/README.md`). A degrees-of-freedom correction in the spirit of Bell and McCaffrey's for CR2, or an F-type reference distribution for the quadratic form, would address it. The size of CR2 with `cluster` has not been checked either.

**One-sided and equivalence omnibus tests.** For example, a monotonicity test (H0: the AMEF is nonpositive everywhere) or a test that the ADRF stays within a margin of a reference value. The simultaneous-band machinery already in `summary.curve_est()` is most of what's needed.

**Uncertainty from smoothing-parameter selection in `gam()` fits** *(new)*. The unconditional variance treats the smoothing parameters as fixed, and on `nhanes3lead` it is smaller than a bootstrap that reselects them (0.066 against 0.086 at `logBLL = 1`). An option that adds a correction like the one behind `vcov.gam(unconditional = TRUE)` would close most of the gap; extended families (`nb()`, `tw()`) currently fall back to the conditional variance. Details are in `_dev/TASKS.md`.

**Parallel bootstrapping.** Documenting `fwb.args = list(cl = <n>)` in `?adrf` and the vignette would help with large models; `fwb()` already supports it.

## Estimands

**Doubly robust estimation.** The pseudo-outcome kernel-smoothing estimator of Kennedy et al. (2017), taking the user's outcome model and a treatment-side input (GPS values or stabilized weights). Planned in `_dev/plan-doubly-robust.md`: first an exported constructor that turns curves estimated elsewhere into an `effect_curve` (with the `adrf()` method for *discoverHTE*'s `<erf>` objects living in *discoverHTE*), then possibly Kennedy's estimator. The constructor is useful on its own, since it would give curves from any method the plotting, contrasts, and omnibus tests of *adrftools*.

**Multi-category outcomes.** Ordinal and multinomial models (`polr`, `ordinal_weightit`, `multinom_weightit`, `nnet::multinom`, *brms* `categorical`) return one prediction per category. `adrf()` now gives an informative error for them. Planned (2026-11-24) in `_dev/plan-categorical-outcomes.md`: the curves for all categories are estimated at once, with the category stored as a column of `.by_grid`, so `subset`, `summary()`, and `plot()` select a category as they select a subgroup. Sums across categories (e.g., the probability of the top two categories of an ordinal outcome) would be a natural addition once that exists.

**Thresholds.** The treatment value at which the ADRF crosses a given outcome value, or the minimum effective dose, with an interval from inverting the simultaneous band.

**Doubly robust shift effects** *(new)*. Shift curves currently use only the outcome model. A treatment-density ratio between the shifted and observed treatment gives a weighting or doubly robust estimator of the same estimand, as in the modified treatment policy literature. It would share the treatment-side input of the doubly robust ADRF above, so it is best considered alongside it.

**`ade()` for fitted models** *(new)*. `ade()` takes a shift curve, so a user who wants only the number must first call `shift_curve()`. A method for fitted models that calls `shift_curve()` with a small range internally was left for later in `_dev/plan-shift-effects.md`.

**Survival outcomes.** `coxph_weightit` and other survival models with `type = "survival"` at a chosen time, or restricted mean survival time.

**Effect curves over a continuous moderator** *(new)*. `by` handles discrete moderators; for a continuous one, users must cut it into groups. ADRFs at chosen values of the moderator (e.g., its quartiles), with the moderator set to each value for all units, would show how the curve changes with it. This is a special case of the multivariate-treatment work and could reuse its grid.

## Data handling

**A target-population argument.** Compute the ADRF or shift curve in a supplied `newdata` (for transporting to another population) instead of the estimation sample, with the conditional variance by default, since the unconditional variance would need the target's sampling design.

**A positivity diagnostic** *(new; deferred in `_dev/plan-shift-effects.md`)*. For shift curves, the share of shifted treatment values outside the observed support (or beyond its .01 and .99 quantiles) at each shift, overall and by subgroup; in `nhanes3lead`, a shift of $-0.5$ SD moves 10.8% of units below the observed minimum. For the ADRF, an estimate of the conditional treatment density at each grid value, averaged over units, would show where the curve rests on extrapolation. Either could be printed by `summary()` or drawn by `plot()` (see the next section).

## Output and plotting

**Treatment distribution on plots.** A rug or density strip beneath the curve, so readers can see where the curve is extrapolating. `adrf()` already stores the treatment density used by `integrate = "empirical"` (`.treat_density`), so the strip needs no extra data. For shift curves, the positivity diagnostic above would be the natural equivalent.

**Several curves on one plot** *(new)*. A way to overlay effect curves from different outcome models (e.g., a GLM and a BART fit) or estimators, for checking sensitivity to the model specification. Plotting needs only the estimates and bands; contrasts between the curves would need their joint covariance, which is available for the bootstrap with shared weights but not otherwise.

**`tidy()`, `glance()`, and `as.data.frame()` methods.** These would need *generics*, which was removed from Suggests because nothing used it (audit R44). `confint()` methods for `curve_est` and `curve_projection` would also be natural.

**Pointwise and simultaneous bands together.** An option to draw both, as is common in the functional data literature.

## Package internals

**Parallel predictions for slow models** *(new)*. Predictions for each treatment value are independent, and for some models they dominate the run time: a *bartisan* soft-rule fit on `nhanes3lead` takes about 1.7 s per grid point, so the vignette uses `n = 21`. Running `.data_at()` predictions across treatment values in parallel (with *mirai* or `parallel::mclapply()`, opt in) would help Bayesian and machine-learning models. The fast path for `lm()` and `glm()` models is already quick.

**Expose the interpolation settings** (`degree`, `kernel`, `bw`) in `adrf()` or as an option. (Capping the degree for small grids, audit issue 9, is done.)

**Plot snapshot tests.** *vdiffr* snapshots would catch visual regressions that the current `layer_data()` checks miss (labels, legends, colors).

## Implemented since 2026-11-12

- **Respecting the `vcov` setting of *WeightIt* fits** (audit M16): fits with `vcov = "HC0"` use `estfun(asympt = FALSE)`.
- **Degrees of freedom for multiply imputed data** (audit M12): done differently from the Barnard–Rubin df proposed here, with a single Reiter (2007) df for the whole curve, or that of Li, Raghunathan, and Rubin (1991) when the complete-data df are infinite.
- **A leverage-adjusted unconditional variance** (audit M19, 2026-11-16): HC2 by default (HC1 when the model has no hat values) and CR2 with `cluster` (CR1 for models other than `lm()` and `glm()`). The remaining size error is the small-sample item above.
- **Shift and incremental interventions** (2026-11-17): `shift_curve()`, with additive, multiplicative, and user-supplied shifts, and `vignette("shift-effects")`.
- **Average derivative effect** (2026-11-17): `ade()` computes $E[\partial_a \mu(A, X)]$, the derivative at 0 of a shift curve. The AMEF averaged over the treatment density, the version proposed here, was set aside as hard to interpret; `curve_projection(integrate = "empirical")` gives a density-weighted slope of the ADRF instead.
- **Density-weighted omnibus tests and projections** (2026-11-16): the `integrate` argument of `summary.effect_curve()` and `curve_projection()`, which also covers part of the quantile-grid idea.
- **`subset`, `by`, `wts`, and `cluster` from the full data**, not just model predictors (audit issue 31).
- **Memory for large samples** (audit R38, 2026-11-16): the model matrix is built one treatment value at a time from a model frame of the whole grid, which cut peak memory on the `?adrf` example from 447 MB to 169 MB at `n = 201`.
- **Unused helpers** (audit M13): `simultaneous_conf_level()` moved to `_dev/`; the other helpers and `make_rmvt()` are kept for possible future features, such as multivariate treatments.
- **Plot data access** (2026-11-25): an `as.data.frame()` method for effect curves returns the estimates at the grid points and the confidence bands `plot()` draws, computed by the same code.
- **Richer print output** (2026-11-25): `print()` reports the number of grid points and subgroups, the sample size and number of imputed datasets, and the small-sample adjustment (HC2, HC1, CR2, CR1, or the survey design), the type of bootstrap or named covariance matrix, and the number of replicates or draws.
