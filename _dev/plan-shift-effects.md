# Plan: shift-response curves and the average derivative effect (2026-11-16)

Adds `shift_curve()`, which estimates the mean outcome when every unit's treatment is shifted from its observed value, and `ade()`, which computes the average derivative effect (ADE) from a shift curve. Implements two items from `_dev/feature-ideas.md` ("Shift and incremental interventions" and "Average derivative effect", the latter reinterpreted as below).

**Decided (2026-11-16):** the ADE is computed from a shift curve, not by a function of its own that takes a fitted model; users can supply their own shift function; no positivity diagnostic in this version; the argument choosing the kind of shift is not named `type`, which `adrf()` already passes to the predictions; the curve's argument is named `shift`.

**Implemented (2026-11-16, branch `shift-curve`), as planned except:**

- The curve's argument is `shift`. Instead of a separate `.x_name` attribute, the `.treat` attribute holds `"shift"` for shift curves (it names the curve's argument everywhere it is used) and the treatment's name is stored in `.shift$treat`, which the print and plot labels use. This left the ~55 uses of `.treat` unchanged.
- `.curve_type` stays `"ADRF"`/`"AMEF"`; `.curve_type_label()` gives "shift curve" and "shift curve derivative" for printing and titles.
- The shift function is not stored in the object (it would carry its environment); `.shift` keeps its label, type, the treatment SD, the SD units of the default range, and whether $g(A_i, 0) = A_i$.
- Shifted values at which the model cannot be evaluated give an error both when the model frame drops rows and when the estimates are not finite (`log()` of a negative value gives `NaN` predictions, not dropped rows).
- `ade()` uses the interpolated curve, as planned; the fallback was not needed (see below). Its output is `ade_est`, labeled `"ADE"` through a new `.est_names` attribute read by `.get_curve_est_labels()`.
- **A bug found on the way:** `amef()` returned all `NA` with 101 or more grid points, for ADRFs as well. The local polynomial in `get_locpoly_w()` used a basis orthogonal over the whole grid, which is nearly collinear within a narrow kernel window, so `qr()` dropped columns and `qr.coef()` returned `NA`. Each local fit now uses coordinates centered at the interpolated point and scaled by the window, and takes the intercept. The fit is the same (the bases span the same space); it was introduced earlier in this development cycle (commit a7d492c replaced the normal equations with QR), so it gets no NEWS item.

**Verification** (`_dev/variance-verification/verify-shift.R`, results in `shift-results.csv`): the unconditional variance matches a hand-computed HC2 influence function and *geex* stacked M-estimation (HC0) to every printed digit, including a logistic model with target weights that differ from the model weights and a user-supplied shift within a subgroup; a fractional weighted bootstrap written independently (10,000 replicates) agrees within 1% for the curve and 1.5% for the ADE. The ADE from the interpolated curve differs from an exact central difference of the predictions by less than $10^{-7}$ standard errors for polynomial, spline, and logistic models at 21, 51, and 101 grid points (up to $2 \times 10^{-5}$ for splines at 21), and matches `marginaleffects::avg_slopes()` (estimate and conditional HC0 standard error). Tests are in `tests/testthat/test-shift-curve.R`.

**Not done:** the vignette section; the open questions below are unchanged.

## Estimands

Write $\mu(a, x) = E[Y \mid A = a, X = x]$ and let $g(a, \delta)$ be the treatment a unit with observed treatment $a$ receives under a shift of size $\delta$. The shift-response curve is

$$\psi(\delta) = E[Y(g(A, \delta))] = E[\mu(g(A, \delta), X)],$$

averaged over the joint distribution of $(A, X)$ in the target population (the units in `subset`, weighted by `wts`, within subgroups defined by `by`). The ADRF is the special case $g(a, \delta) = \delta$.

- Additive shift: $g(a, \delta) = a + \delta$.
- Multiplicative shift: $g(a, \delta) = a(1 + \delta)$, so $\delta = 0.1$ is a 10% increase.
- User-supplied: any function `function(a, delta)`, e.g., `function(a, delta) pmin(a, delta)`, which caps every unit's treatment at $\delta$.

When $g(a, 0) = a$, the derivative of the curve at 0 is

$$\psi'(0) = E\left[\dot{g}(A)\, \partial_a \mu(A, X)\right], \qquad \dot{g}(a) = \partial_\delta g(a, \delta)\big|_{\delta = 0}.$$

For the additive shift this is the ADE, $E[\partial_a \mu(A, X)]$, the effect of an infinitesimal shift of everyone's treatment (Hines, Diaz-Ordaz, and Vansteelandt, 2026). For the multiplicative shift it is $E[A\, \partial_a \mu(A, X)]$, the change in the mean outcome per unit proportional increase in the treatment (multiply by .01 for a 1% increase). For other functions it is the conditional slope averaged with weights $\dot{g}(A)$.

The ADE differs from the AMEF: the AMEF at $a$ averages the conditional slope at a fixed treatment value over the distribution of $X$, whereas the ADE averages it at each unit's own treatment value. The two agree only when $A$ is independent of $X$. The ADE requires positivity only near the observed treatment values, which is its main appeal; it still relies on the outcome model being correct near those values.

## The core change

`.adrf_fit_one()` builds counterfactual data by setting every unit's treatment to `values[k]`. It will instead set it to $g(A_i, v_k)$, with $g(a, v) = v$ for the ADRF. Everything downstream already averages predictions at counterfactual treatment values within subgroups, so it applies unchanged:

- the fast path (the model frame evaluated on the whole grid, the model matrix built one value at a time);
- the unconditional variance, with the covariate-sampling term (which now also captures the sampling of $A$, since the counterfactual treatment depends on it), HC2/CR2, and `svyglm` design-based variances;
- the conditional variance through `marginaleffects::get_vcov()`;
- the bootstrap, including the cache of model-matrix columns that are the same at every value, which is found by comparing columns and so stays correct;
- posterior draws, which are computed with `.data_at()`;
- multiply imputed data (each imputation shifts its own treatment values) and Rubin's rules with the Reiter degrees of freedom;
- `by`, `subset`, `wts`, and `cluster`.

Specifically, in `R/adrf.R`:

- `.effect_curve_internal()` gains a `shift_fun` argument (`NULL` for the ADRF). When it is supplied, the grid comes from a new `process_shift()` in place of `process_range()`, `.treat_density` is not computed, and `shift_fun` is passed to `.adrf_fit_one()`.
- `.adrf_fit_one()` computes the counterfactual treatment once as a vector over the grid (unit within value, the order of `.data_grid()`), and `.data_at()`, `.data_grid()`, and `data_check` take their treatment values from it. The vector is the size of one column of the grid, which is already built.
- The `.quick` shortcut (`R/adrf.R:378`), which collapses the data to one row when the treatment is the only predictor, is skipped for shifts, because units no longer share a counterfactual treatment value.
- A shifted value at which the model cannot be evaluated (e.g., a negative value under `log(treat)`) gives an error that names the shift and suggests a narrower `shift`, rather than the current message about the grid.

## `shift_curve()`

```r
shift_curve(x, treat, shift = NULL, shift_fun = "additive", vcov = "unconditional",
            cluster = NULL, type = "response", data = NULL, subset = NULL,
            by = NULL, wts = NULL, n = 51, fwb.args = list(), ...)
```

A generic with `default` and `mira` methods, mirroring `adrf()` with `range` replaced by `shift` and `shift_fun`.

**`shift_fun`.** `"additive"`, `"multiplicative"`, or a function of two arguments, the vector of observed treatment values and a single shift value, returning a numeric vector of the same length. It is checked on the observed treatment for every grid value: wrong length or non-finite values give an error. `"multiplicative"` requires a nonnegative treatment. At creation, `shift_curve()` records whether $g(A_i, 0) = A_i$ for every unit (always true for the two built-in shifts), which `ade()` requires.

Alternative names: `policy` (after modified treatment policies) or `shift_type`. `shift_fun` reads naturally when a function is supplied.

**`shift`.** Two numbers give the bounds of $\delta$ in the units of $\delta$ (treatment units for an additive shift, proportions for a multiplicative one, whatever the function takes otherwise). A single number $s$ gives the symmetric range $[-s, s]$, in standard deviations of the treatment for an additive shift and in the units of $\delta$ otherwise. This parallels `range` in `adrf()`, where a single number is also on a relative scale.

**Default shift size.** `shift = NULL` means 0.5 for an additive shift, i.e., $\pm 0.5$ standard deviations of the treatment, and gives an error asking for `shift` with any other `shift_fun`.

- The standard deviation is of the treatment as it enters the model (e.g., `logBLL`, not BLL), computed over the units in `subset` and weighted by `wts`, pooled over subgroups so every subgroup has the same grid and their curves can be contrasted, and over imputations, as `process_range()` does with quantiles.
- Standard deviation units are scale-free and familiar ("a half-SD increase"). The IQR would be more robust to skew but is less familiar, and for most treatments the two are close (for `logBLL` in `nhanes3lead`, SD = 0.726 and IQR/1.349 = 0.733).
- A half SD keeps most shifted values within the observed support. In `nhanes3lead`, $+0.5$ SD moves 0.4% of units above the observed maximum and 3.4% above the .99 quantile; $+1$ SD moves 1.2% and 8.4%. Downward shifts are worse there, because `logBLL` has many units near its minimum: $-0.5$ SD moves 10.8% below the observed minimum ($-0.25$ SD 9.1%, $-1$ SD 18.9%). This is the kind of thing the deferred positivity diagnostic would report.
- With confounding, the treatment varies less within covariate strata than marginally, so a shift of a marginal SD is larger in conditional terms. In `nhanes3lead`, $R^2$ of `logBLL` on the covariates is 0.19, so 0.5 marginal SD is 0.55 conditional SD; with $R^2 = .5$ it would be 0.71. That argues for a modest default.
- The ADE does not depend on `shift` beyond the accuracy of the derivative (see `ade()`), so the default matters mostly for the plotted curve and the omnibus tests.
- With the default odd `n`, a symmetric range puts 0 on the grid.

**Output.** An `effect_curve` with `.curve_type = "ADRF"` (it is a curve of mean outcomes on the outcome scale, so `transform`, the default `null`, and `.is_pure_adrf()` behave as for an ADRF), an added class `shift_curve`, and a new attribute `.shift`: a list with the label of the shift function, the treatment SD, the identity-at-zero flag, and the range in both units. `.make_fun()` gains `.shift` as a copied attribute and adds the `shift_curve` class whenever it is non-`NULL`, so `amef()`, `reference_curve()`, `curve_contrast()`, and `curve_projection()` keep it.

**The curve's argument.** `.make_fun()` names the callable's argument after `.treat` (`rlang::fn_fmls_names(fn)[1L] <- .treat` in `R/effect_curve.R`), which would make `sc(logBLL = .2)` read as a treatment value. A new attribute `.x_name` (default `.treat`; `"shift"` for shift curves) will name the argument, the column of `curve_est` objects, and the plot's x variable, while `.treat` stays in labels (x-axis "Shift in logBLL", or "Proportional change in logBLL", or "Shift (`shift_fun`)" for a function). There are about 55 uses of `.treat` in `R/`; each needs sorting into "names the curve's argument" or "describes the treatment".

**Interaction with existing functions.**

- `plot()`, `print()`, `summary.curve_est()`, `point_contrast()`: work, with the labels above. `print()` reports the range in treatment units and, for the default, in SD units.
- `summary.effect_curve()`: omnibus tests work; `hypothesis = "flat"` tests that no shift in the range changes the mean outcome. `integrate = "empirical"` gives an error, because there is no treatment density over shifts.
- `reference_curve(sc, 0)`: the shift effects $\psi(\delta) - \psi(0)$. Consider making 0 the default reference for shift curves.
- `amef(sc)`: the derivative over $\delta$.
- `curve_contrast()`: contrasts of shift curves between subgroups.
- `curve_projection()`: works with `integrate = "uniform"`.

## `ade()`

```r
ade(x, subset = NULL, eps = 1e-5)
```

Takes a shift curve (including a reference curve or curve contrast made from one) and returns the derivative at $\delta = 0$ for each subgroup or contrast, i.e., `amef(x)(shift = 0)`, relabeled. Applied to a curve contrast, it gives differences in ADEs between subgroups, with their joint covariance. The variance, bootstrap replicates, and posterior draws come from the contrast-matrix machinery `amef()` already uses.

Errors:

- `x` is not a shift curve: the message points to `shift_curve()`, and to `curve_projection(integrate = "empirical")` for a one-number summary of an ADRF.
- 0 is outside the range of the curve.
- The shift function does not satisfy $g(A_i, 0) = A_i$ (e.g., the cap above).
- `x` is an AMEF or a projection.

**Output.** A `curve_est` subclass, `ade_est`, without the shift column, so `summary()`, `coef()`, `vcov()`, and `marginaleffects::hypotheses()` work. `print()` labels the estimate "Average derivative effect" for an additive shift and names the shift otherwise.

**Accuracy.** The derivative comes from the local polynomial interpolation of the grid, not from an exact finite difference of the predictions. With the default range and `n`, the grid spacing is 0.02 SD, so for smooth models the error should be far below the standard error; this needs to be checked (see Tests). For models whose predictions are step functions of the treatment (*dbarts*, *stan4bart*, and *bartisan* fits), an exact finite difference with a tiny step is zero or erratic, and the grid-based derivative is a smoothed version over the grid spacing, which is an advantage worth documenting. If the check fails for smooth models, the fallback is to have `shift_curve()` also estimate $\psi(\pm \text{eps})$ in the same pass and store those two estimates (with their replicates) in `.shift`, for `ade()` to use as an exact central difference.

## What this offers beyond `marginaleffects::avg_slopes()`

For the documentation and vignette, stated carefully:

- The unconditional variance accounts for sampling of the covariates and of the treatment, which the ADE averages over. As far as I know, `avg_slopes()` uses the model's covariance by default, which treats them as fixed.
- HC2/CR2 small-sample adjustments, design-based variances for `svyglm` fits, and, for `glm_weightit()` fits, M-estimation that accounts for estimating the weights.
- The whole shift curve, with simultaneous bands and omnibus tests, not only its slope at 0.
- Subgroup ADEs with simultaneous inference and their contrasts.
- A derivative that is well defined for tree-based models.
- The ADE alongside `curve_projection(integrate = "empirical")` in the same object system; the two agree when the treatment is unconfounded, so a gap between them indicates that the conditional slopes vary with the covariates in a way related to the treatment.

Not a distinguishing feature: bootstrap inference, which `marginaleffects::inferences()` also provides.

## Documentation

Use the `/r-doc-style` skill for all of it.

- `?shift_curve` and `?ade`: the estimands, the default shift, positivity (shifts need overlap only near the observed values, but large shifts give that up), and the difference between the ADE and the AMEF.
- `?amef`: its `@seealso` currently calls `avg_slopes()` "similar to the AMEF"; `avg_slopes()` estimates the ADE, so point to `ade()` instead and say how the two differ.
- `?effect_curve`, `?adrf`, `?summary.effect_curve` (`integrate = "empirical"` not available for shift curves), `?reference_curve`.
- Vignette: a section on shift curves and the ADE with `nhanes3lead`, including the many units near the minimum of `logBLL`.
- References: Díaz Muñoz and van der Laan (2012, *Biometrics*, doi:10.1111/j.1541-0420.2011.01685.x); Haneuse and Rotnitzky (2013, *Statistics in Medicine*, doi:10.1002/sim.5907); Hines, Diaz-Ordaz, and Vansteelandt (2026, *Biometrika*, doi:10.1093/biomet/asag012); Lundberg and Brand (2026, *Sociological Methodology*, doi:10.1177/00811750261459364).
- NEWS.

## Tests

`tests/testthat/test-shift-curve.R`, using the `/r-patrick` skill for the cases over shift functions and variance types.

- **ADRF unchanged.** `adrf()` results identical to before the change, with and without the `.quick` path.
- **Point estimates.** $\psi(\delta)$ equals the weighted mean of `predict()` on data with the treatment shifted, for additive, multiplicative, and a user function; with `subset`, `by`, and `wts`; for `lm()`, `glm()`, a model on the fallback (non-fast) path, and a treatment-only model (checking `.quick` is off).
- **Conditional variance.** Matches `marginaleffects::avg_predictions()` on the shifted data with the same `vcov`.
- **Unconditional variance.** Matches a hand-computed influence function for a small `lm()`; checked against the bootstrap and against *geex* or *stdReg* in `_dev/variance-verification/`, as for the ADRF.
- **ADE.** `ade(sc)` matches `avg_slopes()` in point estimate (relative tolerance of about 1e-6) and conditional standard error with the same `vcov`; `ade(sc)` equals `amef(sc)(shift = 0)`; the error of the interpolated derivative relative to an exact central difference, for polynomial, spline, and logistic models at the default `n`, is below 1% of the standard error and shrinks as `n` grows.
- **Multiplicative.** `ade()` equals the weighted mean of $A_i$ times the unit-level slopes.
- **Contrasts.** `ade(curve_contrast(sc))` equals the difference of the subgroup ADEs, with the covariance from the joint estimates.
- **Defaults and grid.** The default range is $\pm 0.5$ weighted SD, pooled over subgroups and imputations, with 0 on the grid; a single number for a multiplicative shift is a proportion.
- **Errors.** `ade()` on an ADRF, on a curve whose range excludes 0, and on a shift function that is not the identity at 0; a shift function returning the wrong length or non-finite values; a multiplicative shift with a negative treatment; a shift outside the domain of a transformation in the formula; `integrate = "empirical"` on a shift curve; `shift = NULL` with a non-additive `shift_fun`.
- **Other inference.** Bootstrap, multiply imputed data, and a Bayesian model (`skip_if_not_installed()`), each checking that the estimates and the shape of the replicates or draws are right.

## Open questions

1. The default shift: $\pm 0.5$ SD as proposed, or $\pm 0.25$ SD to stay further within the support.
2. The name `shift_fun`, or `policy`.
3. Whether a single number for `shift` should mean SD units only for additive shifts (proposed) or never, with SD units used only for the default.
4. Whether to give `reference_curve()` a default reference of 0 for shift curves.
5. A positivity diagnostic (the share of shifted values outside the observed support, overall and by subgroup) is deferred, not rejected.
6. A method of `ade()` for fitted models, which would call `shift_curve()` with a small range, could be added later for users who want only the number.
