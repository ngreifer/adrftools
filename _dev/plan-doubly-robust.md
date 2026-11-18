# Plan: doubly robust ADRFs (2026-11-18)

Investigates adding the doubly robust (DR) estimator of Kennedy et al. (2017) while keeping *adrftools* a package for visualizing and testing effect curves rather than for fitting models. Item "Doubly robust estimation" in `_dev/feature-ideas.md`. No package code has been changed.

**Recommendation.** Do the first option now: an exported constructor that turns estimates computed elsewhere into an `effect_curve`, with the `adrf()` method for *discoverHTE*'s objects living in *discoverHTE*. Then, if a native DR estimator is wanted, implement Kennedy's estimator from the user's outcome model and the user's estimate of the generalized propensity score (GPS) at each unit's observed treatment, which balancing weights can supply. That avoids most of the density problem: Kennedy's estimator evaluates the conditional density of the treatment only at the observed treatment values, never on a grid. TMLE and the estimator of Colangelo and Lee (2026) need the conditional density at every grid value for every unit, so they are harder here, not easier. Two later additions reuse the same inputs and suit the package's focus on testing, and neither needs a bandwidth: Westling's (2022) DR test of a flat ADRF, and DR projections of the ADRF onto parametric models.

A correction to the request: *discoverHTE*'s objects have class `<erf>` (from `estimate_erf()`), not `erf_estimate`.

## Notation and Kennedy's estimator

Write $\mu(a, x) = E[Y \mid A = a, X = x]$, $\pi(a \mid x)$ for the conditional density of the treatment (the GPS), $\varpi(a) = E[\pi(a \mid X)]$ for its marginal density, and $m(a) = E[\mu(a, X)]$, which is the ADRF $\theta(a)$ under the usual assumptions. `adrf()` currently estimates $\hat m(a) = P_n \hat\mu(a, X)$.

Kennedy et al. (2017, Section 3; read from the LaTeX source, arXiv 1507.00747) regress the pseudo-outcome

$$\hat\xi = \frac{Y - \hat\mu(A, X)}{\hat\pi(A \mid X)}\,\hat\varpi(A) + \hat m(A)$$

on $A$ by local linear kernel regression at each anchor $a$, with bandwidth $h$. $E[\xi \mid A = a] = \theta(a)$ if either $\mu$ or $\pi$ is correct. Their sandwich variance for the local linear coefficients at $a$ uses the estimated influence function

$$\hat D_{ha}^{-1}\Big[g_{ha}(A) K_{ha}(A)\big\{\hat\xi - g_{ha}(A)^\top \hat\beta_h(a)\big\} + \int g_{ha}(t) K_{ha}(t)\big\{\hat\mu(t, X) - \hat m(t)\big\}\hat\varpi(t)\,dt\Big],$$

with $g_{ha}(t) = (1, (t - a)/h)^\top$, $K_{ha}(t) = K\{(t - a)/h\}/h$, $\hat D_{ha} = P_n\{g_{ha}(A) K_{ha}(A) g_{ha}(A)^\top\}$, and $\hat\varpi(t) = P_n \hat\pi(t \mid X)$. The intervals cover the smoothed curve $\theta_h(a) = \theta(a) + O(h^2)$ rather than $\theta(a)$ unless one undersmooths; Kennedy et al. say so and center their intervals at $\theta_h$. The bandwidth is chosen by leave-one-out cross-validation of the pseudo-outcome regression. *discoverHTE* implements all of this (`estimate_pseudo_Y()`, `make_erf()`, `get_bw()` in `R/estimate_erf.R`), with the anchor covariance computed from the cross-products of each unit's influence values across anchors.

The conditional density enters only through $\hat\pi(A_i \mid X_i)$, at each unit's observed treatment. Everything else is the marginal density $\hat\varpi$ (in the pseudo-outcome and in the variance's integral term) or predictions from the outcome model. Kennedy et al. estimate $\hat\varpi(t)$ as the average of $\hat\pi(t \mid X_i)$, which does need the density everywhere, but any consistent estimate of the marginal density should serve: when $\pi$ is correct and $\mu$ is not, an error in $\hat\varpi$ enters $E[\hat\xi \mid A]$ multiplied by the outcome model's error, and when $\mu$ is correct the residual has mean zero whatever multiplies it. A kernel density estimate of the treatment is even consistent when the GPS model is wrong. *adrftools* already computes one, weighted and by subgroup, for `integrate = "empirical"` (`.get_treat_density()`). So the treatment side reduces to one number per unit: the GPS, or equivalently the stabilized weight $\hat\varpi(A_i)/\hat\pi(A_i \mid X_i)$. (This argument is mine, not Kennedy et al.'s; the simulations below should confirm it.)

## Inputs Each Method Needs

| Method | Outcome model | Treatment side | Bandwidth | Rate |
|---|---|---|---|---|
| g-computation (current `adrf()`) | $\hat\mu$ on the grid | none | none | $\sqrt n$ under a correct parametric model |
| Kennedy et al. (2017), local linear | $\hat\mu$ at the observed $(A_i, X_i)$ and on a grid spanning the support of $A$ | GPS at $A_i$ only | yes | $\sqrt{nh}$ |
| Kennedy et al. (2017), TMLE (their appendix) | same | clever covariate at $(X_i, t)$ for every grid value $t$, so the GPS everywhere; one fluctuation fit per anchor | yes | $\sqrt{nh}$ |
| Colangelo and Lee (2026) | $\hat\mu$ on the grid | GPS at every anchor $t$ for every unit | yes | $\sqrt{nh}$ |
| DR projection onto a parametric model in $A$ | as Kennedy | GPS at $A_i$ only | none | $\sqrt n$ |
| Westling (2022) test of a flat ADRF | $\hat\mu(A_j, X_i)$ for all pairs (a grid would do) | $\hat\pi(A_i \mid X_i)/\hat\varpi(A_i)$ only | none | $\sqrt n$ |
| DR shift curve (Haneuse and Rotnitzky, 2013) | $\hat\mu$ at $A_i$ and at the shifted value | ratio of the shifted to the natural treatment density at $A_i$ | none | $\sqrt n$ |
| Variance-weighted ADE (Hines et al., 2026) | $E[Y \mid X]$, not $\mu(a, x)$ | $E[A \mid X]$ predictions | none | $\sqrt n$ |

The treatment-side inputs fall into two tiers. **Tier 1** is one number per unit, the GPS at the observed treatment or a stabilized weight; it suffices for Kennedy's estimator, DR projections, and Westling's test. **Tier 2** is a conditional density that can be evaluated at any treatment value; TMLE, Colangelo and Lee, and DR shift curves need it. Tier 1 is what users of *WeightIt* already have. Tier 2 is what makes the density "not straightforward for all models."

## Option 1: Curves Estimated Elsewhere

### The Constructor

An exported function (name open; `as_adrf()` for now) builds an ADRF `effect_curve` from estimates on a grid:

```r
as_adrf(estimate, values, treat, vcov = NULL, influence = NULL,
        cluster = NULL, by = NULL, treat_values = NULL, weights = NULL,
        family = NULL, response = NULL, df = Inf, label = NULL)
```

- `estimate`, `values`, `treat`: the estimates at equally spaced grid values and the treatment's name. With `by`, a data frame with one row per subgroup, and the estimates ordered by subgroup, then value (the package's internal order).
- `vcov` or `influence`: the covariance of the estimates, or an $N \times k$ matrix of per-unit influence values from which *adrftools* computes it. Influence values make `cluster` possible (summed within clusters, with the $G/(G-1)$ factor); a covariance alone does not.
- `treat_values` and `weights`: the observed treatment (and sampling weights), used only to estimate the treatment density for `integrate = "empirical"`. Without them, `integrate = "empirical"` errors.
- `family`: e.g., `binomial()` for a probability, so `summary(transform = TRUE)` computes intervals on the logit scale. A local linear DR estimate can leave $[0, 1]$ (Kennedy et al. note this); `.check_transform_values()` already falls back to the identity with a warning.
- `label`: printed as the inference type, e.g., "doubly robust (Kennedy et al., 2017)". `.vcov_type` can hold any string: the only branches on it test for `"none"`, `"bootstrap"`, and `"posterior"` (checked by grep), and `.get_vcov_type_name()` prints the string.

It calls `.make_fun()` with `.curve_type = "ADRF"`, so `.make_fun()` needs no change. Checks: 2 to 1000 increasing, equally spaced values (the interpolation in `get_locpoly_w()` uses a bandwidth of half the mean spacing); finite estimates; a symmetric, positive semidefinite covariance of the right dimension; influence values with one column per estimate.

What works on the result: printing, plotting, evaluation between grid points, `summary()` of estimates (pointwise and simultaneous), omnibus tests, `amef()`, `reference_curve()`, `point_contrast()`, `curve_projection()` and `anova()`, and `curve_contrast()` when `by` is supplied. What doesn't: anything that re-estimates (`subset`, `by`, `wts`, `range`, `n`, the bootstrap, a different `vcov`), and `shift_curve()` and `ade()`, which need a model.

Allowing posterior draws (the `"posterior"` path needs only a draws matrix) would cover Bayesian estimates made elsewhere; bootstrap replicates would need wrapping in an *fwb*-like object for percentile intervals. Both are open questions.

### The Method for *discoverHTE*'s `<erf>` Objects

*discoverHTE* is not on CRAN (checked against CRAN on 2026-11-18), and its own `_dev/TASKS.md` notes its GitHub URL does not resolve publicly yet, so *adrftools* cannot list it in Suggests. The method should live in *discoverHTE*, registered with `@exportS3Method adrftools::adrf` (delayed registration; *discoverHTE* already Suggests *adrftools*), and call the constructor.

Everything it needs is in documented components of `<erf>` (`?estimate_erf`): `x$erf(x$values)` returns the estimates and covariance exactly at the anchors (no interpolation when the requested values are the anchors); `x$treat` carries the treatment's name; `x$weights`, `x$formula`, `x$bw`, and `x$cluster` are stored. Because the method lives with the nuisance fits, it can go further than conversion: `range` and `n` by re-smoothing at new anchors without refitting the nuisance models (as `update(bw = )` does), and later `subset` and `by`, with pseudo-outcomes recomputed within subgroups as the trees do. `make_erf()` computes the per-unit influence matrix and discards it; keeping it (about 1 MB for `nhanes3lead` at 51 anchors) would let *adrftools* apply `cluster`.

Details for the method:

- *discoverHTE* defaults to `range = 1`, the full exposure range; *adrftools* defaults to .95 because the variance at the extremes dominates the omnibus tests. Re-smoothing on .95 by default would match `adrf()`.
- ERFs of the derivative (`d = 1`) are AMEFs; map them to `.curve_type = "AMEF"` or reject them.
- Node ERFs extracted from an `<erf_tree>` carry estimates but no nuisance fits, so they can only be converted. Leaves of an honest tree are estimated on disjoint units, so their covariance is zero under Kennedy's variance, and `curve_contrast()` between two leaves would follow from stacking them with `by`.

### The Limits of This Route

Other producers rarely supply a joint covariance. *npcausal*'s `ctseff()` (Kennedy's own implementation) returns pointwise standard errors in its `res` data frame (source checked on GitHub, 2026-11-18), and *CausalGPS* 0.5.1's `estimate_erf()` returns predictions without standard errors. Pointwise standard errors support pointwise intervals only; simultaneous bands, omnibus tests, projections, AMEFs, and contrasts all need the joint covariance. The constructor should require a full covariance (or influence values), and its documentation should say why. The point raised in the request stands for `vcov`, which comes entirely from the producer, but not for `cluster`, which works whenever the producer supplies influence values.

## Option 2: Kennedy's Estimator From the User's Inputs

### Interface

A treatment-side argument on `adrf()` (name open: `gps`, `treat_model`, or `dr`), plus `bw` and `trim`, which error unless it is supplied, the way `adrf()` and `shift_curve()` reject each other's arguments. A separate function is the alternative, but the estimand, the outputs, and nearly all arguments are the same as `adrf()`'s. The argument would accept:

- a numeric vector of GPS values $\hat\pi(A_i \mid X_i)$, which *adrftools* turns into ratios with its own kernel estimate of $\varpi$;
- stabilized weights, as a numeric vector or a `weightit` object, used directly as $\hat\varpi(A_i)/\hat\pi(A_i \mid X_i)$. Here the numerator matters: with a correct GPS and a wrong outcome model, $E[\hat\xi \mid A = a]$ is off by the factor $\tilde\varpi(a)/\varpi(a)$ for a numerator $\tilde\varpi$. *WeightIt*'s numerator under the default `density = "dnorm"` is a normal density, so it is fine only for a roughly normal treatment; `density = "kernel"` uses a kernel estimate. *WeightIt* 2.1.0 does not store the GPS for continuous treatments (checked); having `weightit()` store it would make a `weightit` object a clean Tier-1 input. For a `glm_weightit()` fit, the weights are already in the object, so `dr = TRUE` could use them;
- a fitted treatment model whose conditional density *adrftools* can compute: `lm()`, and `glm()` with gaussian, Gamma, or inverse.gaussian families, possibly *betareg* and *gamlss*, through a small internal generic. This makes the input Tier 2 and permits Kennedy's own $\hat\varpi(t) = P_n \hat\pi(t \mid X)$. As an opt-in fallback for any model with predictions, the density of the standardized residuals can be estimated by a kernel (see Kernel Density Estimation below).

### Computation

1. The observed outcomes $Y_i$, which *adrftools* never uses now (`insight::get_response()`; binomial models with two-column responses need care), and the fitted values $\hat\mu(A_i, X_i)$.
2. $\hat\mu(t, X_i)$ on a nuisance grid spanning the whole observed treatment, not just `range`, because units outside `range` enter the local regressions at anchors near its ends. The fast path in `.adrf_fit_one()` computes predictions on any grid; a second grid of about 100 points roughly doubles the prediction cost.
3. $\hat m(t)$ and $\hat\varpi(t)$ within each subgroup with the target weights (`wts`), interpolated to each $A_i$.
4. Pseudo-outcomes, with optional trimming of the ratio at a quantile (*discoverHTE* trims at .99 by default).
5. The bandwidth by Kennedy's leave-one-out criterion (*discoverHTE*'s `get_LOOCV_obj()` has the closed form for local linear regression) or a supplied value.
6. Local linear fits at the anchors and the $N \times k$ influence matrix, which gives the covariance, and with `cluster` the clustered covariance. The HC2 and CR2 adjustments do not carry over, since they are built on the leverages of the outcome model.
7. Multiply imputed data: one curve per imputation, pooled with `pool_est()` and the existing degrees of freedom. Bayesian outcome models and the bootstrap are out of scope for a first version.

In code, a new `.adrf_dr_one()` beside `.adrf_fit_one()`, called from `.effect_curve_internal()`, which already handles `subset`, `by`, `wts`, `cluster`, and multiple imputation. The statistical code would be ported from *discoverHTE*, since *discoverHTE* Suggests *adrftools* and not the reverse. Longer term, *discoverHTE* could Import *adrftools* for this engine and keep its nuisance methods and trees; that is not needed now.

The result should carry its pseudo-outcomes and the integral-term matrix, so that DR projections and Westling's test (below) can be added later without recomputing them.

### Statistical Points for the Documentation

- Intervals and bands cover $\theta_h$. Takatsu and Westling (2025) give a debiased local linear estimator whose intervals cover $\theta$ itself without undersmoothing, with uniform bands; a possible later refinement (from their abstract; not yet read in full).
- **Tests of a flat or linear ADRF have no smoothing bias under the null.** Local linear regression reproduces linear functions exactly: if $\theta(t) = c_0 + c_1 t$, then $\theta(t) = g_{ha}(t)^\top (c_0 + c_1 a,\, c_1 h)^\top$, so $\theta_h = \theta$. The flat test in `summary.effect_curve()` and linear projections keep their size at any fixed bandwidth; smoothing costs only power. Nulls with curvature (a quadratic projection, for example) are biased by $O(h^2)$, though mostly by a constant in the interior. For a package centered on testing, this matters more than the pointwise coverage problem.
- The cross-validated bandwidth is treated as fixed, as in Kennedy et al.
- There is no cross-fitting, since users supply models fit to the full sample. Kennedy et al.'s Theorem 2 instead requires the nuisance estimators to lie in classes with finite uniform entropy integrals, which covers parametric models and smooth GAMs but not BART or random forests. Users who want cross-fitting can use *discoverHTE*.
- The influence-function variance assumes both nuisance estimators converge fast enough (their condition (f)); with one misspecified, the estimate stays consistent but the intervals are not guaranteed.
- Near-violations of positivity produce extreme ratios; trimming trades a little bias for stability.

### Verification

- Against *discoverHTE*: fit the models `method_glm()` fits (a normal linear model for the treatment; a GLM for the outcome with `splines::ns(A, df = 4)`), pass them to *adrftools* as a Tier-2 treatment model so $\hat\varpi$ is computed the same way, fix the bandwidth, and set `trim = 1`. The estimates and covariance should match `estimate_erf()` up to its integration grid and its smoothing-spline interpolation of $\hat m$ and $\hat\varpi$.
- Calibration under known truth, as the global instructions require: coverage of $\theta_h$, and the size of the flat and linear tests, with (a) both models correct, (b) the outcome model wrong, and (c) the treatment model wrong, comparing power with g-computation. Also check that the kernel estimate of $\varpi$ keeps double robustness, which the argument above claims. Run through *pueue* with the live-progress widget.

## Alternatives Raised in the Request

### TMLE

The appendix of Kennedy et al. gives a TMLE: a logistic fluctuation of $\hat\mu$ with clever covariate $g_{ha}(A) K_{ha}(A)\,\hat\varpi(A)/\hat\pi(A \mid X)$, one fit per anchor, after which $\hat m^*_{ha}(t) = P_n \hat\mu^*_{ha}(X, t)$ needs the clever covariate at $(X_i, t)$, and so the GPS on the whole grid (Tier 2). Its advantage is staying within the bounds of a binary outcome. *adrftools* already respects bounds in intervals with `transform`, and the same appendix suggests simpler fixes (fitting the local regression on the logit scale, or projecting the estimate onto the outcome's range). Not easier; skip.

### Double/Debiased Machine Learning

DML is a recipe, orthogonal scores plus cross-fitting, rather than a separate estimator. Kennedy's pseudo-outcome is already the orthogonal score: it is the uncentered efficient influence function of $\psi = \int m(a)\varpi(a)\,da$ (their Theorem 1). What DML adds is cross-fitting, which means refitting the nuisance models on folds. A package that takes the user's fitted models cannot do that, and parametric or GAM nuisance models don't need it. The R packages for DML target other estimands: *ddml* 0.3.1 (installed here) has `ddml_ate()`, `ddml_att()`, `ddml_late()`, `ddml_plm()`, `ddml_pliv()`, and `ddml_fpliv()`, for binary treatments and the partially linear coefficient, not dose-response curves.

### Colangelo and Lee

Their estimator at an anchor $t$ (Step 2 of their procedure; read from the LaTeX source, arXiv 2004.03036, which may differ from the published version) is the g-computation estimate plus a kernel-weighted residual correction,

$$\hat\beta_t = \frac{1}{n}\sum_{i=1}^n \left\{\hat\mu(t, X_i) + \frac{K_h(A_i - t)}{\hat\pi(t \mid X_i)}\big\{Y_i - \hat\mu(t, X_i)\big\}\right\},$$

with cross-fitting (with a single fold it is the DR estimator without sample splitting). It fits *adrftools* neatly: the first term is what `adrf()` computes now, and the influence function adds one term to the existing unconditional variance. But the correction evaluates the GPS at the anchor $t$, not at the unit's own treatment, so it needs $\hat\pi(t \mid X_i)$ for every unit at every anchor (Tier 2), and its local-constant kernel has the usual boundary bias. Their own GPS estimator (ReGPS) inverts conditional distribution functions estimated by a separate regression for each value of $t$, which is modeling the package should not take on. It could be a second DR option later for users who supply a Tier-2 treatment model; it is not easier than Kennedy's.

### Kernel Density Estimation

A full conditional kernel density estimate (e.g., with *np*) suffers the curse of dimensionality with more than a few covariates. A location(-scale) estimate works better: a kernel density of the standardized residuals from a model for the treatment's mean (and optionally its scale) gives a Tier-2 density from any model with predictions. *WeightIt*'s `density = "kernel"` and *discoverHTE*'s `method_superlearner()` work this way. It assumes the residuals' shape does not depend on the covariates. It is a reasonable opt-in fallback and a small amount of estimation (one kernel density estimate) for the package to own.

## Other Ideas

### Weighted g-Computation and Double Robustness

Users can already fit the outcome model with stabilized GPS weights (`glm_weightit()`, as in the vignette), and that is DR in a limited sense. With weights $\tilde\varpi(A)/\pi(A \mid X)$ for any density $\tilde\varpi$ of $A$ alone, $A$ and $X$ are independent in the weighted population, with $X$ distributed as before. For a GLM with canonical link, the score equations for the intercept and the columns $b(A)$ that depend on the treatment alone then give

$$\int b(a)\big\{\theta(a) - \tilde\theta(a)\big\}\tilde\varpi(a)\,da = 0,$$

where $\tilde\theta$ is the limit of the g-computation estimate. For a linear outcome model, $\tilde\theta$ lies in the span of $1$ and $b$ whatever the covariate terms (interactions $b(A) \otimes c(X)$ included), so $\tilde\theta$ is the $\tilde\varpi$-weighted least-squares projection of $\theta$ onto that span. With a correct GPS model, the g-computation ADRF is therefore consistent for any ADRF in that span (e.g., the span of `splines::ns(logBLL, df = 5)`), however wrong the covariate part of the outcome model, and the numerator only changes the projection's weighting. For other links, only the orthogonality holds. Kennedy's estimator removes the restriction to a span at the cost of a bandwidth. This derivation is mine and was not checked against a paper; confirm it by simulation before it goes in the vignette, where it would help users who weight already.

### Westling's Test of a Flat ADRF

Westling (2022; read from the LaTeX source, arXiv 2001.05344) tests that $\theta$ is constant on the support of $A$ through its primitive $\Omega(a_0) = \int_{-\infty}^{a_0}\{\theta(u) - \gamma\}\,dF(u)$, with $\gamma = \int \theta\,dF$. $\Omega$ is pathwise differentiable, so its one-step estimator is $\sqrt n$-consistent and DR with no bandwidth. The statistic is the supremum (recommended) or an $L_p$ norm of $\sqrt n\,\hat\Omega$, with critical values simulated from a Gaussian process with the estimated influence-function covariance. Its inputs are Tier 1: $\hat\mu(A_j, X_i)$ for all pairs and $\hat\pi(A_i \mid X_i)/\hat\varpi(A_i)$. The machinery for simulated sup-type critical values exists (`.sim_p_value()`, `.mvt_crit()`). It answers the package's central question, whether the treatment matters at all, with a DR test that needs no bandwidth. Westling cross-fits; the non-cross-fitted one-step estimator in the same paper is valid when the nuisance estimators fall in small enough function classes, as parametric fits do. Separately, the same primitive computed from any curve on its grid (a linear map of the estimates) would give a sup-type flatness test for g-computation curves, worth adding to the test-choice study (`_dev/omnibus-test-choice/`).

### DR Projections

Regressing Kennedy's pseudo-outcome on a parametric basis $b(A)$ by least squares estimates the $\varpi$-weighted projection of $\theta$ onto the span of $b$, a $\sqrt n$-estimable parameter needing no bandwidth. Neugebauer and van der Laan (2007) define such projections as nonparametric marginal structural models, and Semenova and Chernozhukov (2021) develop inference for best linear predictors of an orthogonal signal (Colangelo and Lee note that their example uses this pseudo-outcome). This is what `curve_projection(integrate = "empirical")` targets, up to the kernel density estimate and `range`, so a DR `curve_projection()` would give DR $\sqrt n$ versions of the linearity tests (`anova()` of nested projections). Its influence function should be the projection's residual term plus an integral term like Kennedy's with $b(t)$ in place of $g_{ha}(t) K_{ha}(t)$; that is my extrapolation and needs deriving properly.

### DR Shift Curves

`shift_curve()` estimates $\psi(\delta) = E[\mu(A + \delta, X)]$, a regular parameter. Its DR estimator adds to the current plug-in the mean of the residual $Y_i - \hat\mu(A_i, X_i)$ times the ratio of the shifted treatment's density to the natural treatment's density at $A_i$ (Haneuse and Rotnitzky, 2013; check the exact form against the paper before implementing). The variance extends the current unconditional variance, which already has the covariate-sampling term. It needs Tier 2, a density evaluable at $A_i - \delta$; a normal linear treatment model gives the ratio in closed form. Multiplicative and user-supplied shifts would also need the shift to be invertible, with its Jacobian, which `shift_fun` does not require now. *lmtp* and *txshift* implement DR estimators, TMLE among them, for such interventions.

### A Variance-Weighted ADE

Hines et al. (2026; estimators read from the LaTeX source of the preprint, arXiv 2109.13124) give a weighted ADE whose efficient estimator needs only the regressions $E[A \mid X]$ and $E[Y \mid X]$, with no density: $\hat\Psi = \sum_i (a_i - \hat\pi_i)(y_i - \hat\mu_i)/\sum_i (a_i - \hat\pi_i)^2$, where $\hat\pi_i$ and $\hat\mu_i$ are those regressions' predictions. It is the partially linear coefficient that `ddml_plm()` estimates. It needs $E[Y \mid X]$ rather than the outcome model of $Y$ on $A$ and $X$ that users give `adrf()`, so it would be a separate function with a different input pattern; a possible companion to `ade()`, but not part of this plan.

## Open Questions

- The constructor's name (`as_adrf()`, `make_adrf()`, or an `adrf()` method for numeric vectors), and whether it should accept posterior draws or bootstrap replicates.
- Whether the `<erf>` method should only convert, or re-smooth to honor `range`, `n`, and later `by` and `subset`.
- The name of the treatment-side argument for Option 2, which inputs it accepts, and whether `weightit()` should store the GPS for continuous treatments.
- Whether the DR curve uses Kennedy's $\hat\varpi(t) = P_n \hat\pi(t \mid X)$ when a Tier-2 model is supplied (matching *discoverHTE*) or always the kernel estimate (which Tier 1 needs anyway).
- The default bandwidth (Kennedy's cross-validation, as in *discoverHTE*, or a rule of thumb) and the default trimming.
- Whether to pursue Option 2 at all, or stop at Option 1 and point users who want a native DR estimator to *discoverHTE*. Option 2 makes *adrftools* own a bandwidth, a marginal density estimate, and trimming, which is some estimation, though not of models.

## References

All are in the Zotero library; DOIs from the library's metadata.

Colangelo, K., & Lee, Y.-Y. (2026). Double debiased machine learning nonparametric inference with continuous treatments. *Journal of Business & Economic Statistics*, 44(1), 67–79. https://doi.org/10.1080/07350015.2025.2505487

Haneuse, S., & Rotnitzky, A. (2013). Estimation of the effect of interventions that modify the received treatment. *Statistics in Medicine*, 32(30), 5260–5277. https://doi.org/10.1002/sim.5907

Hines, O. J., Diaz-Ordaz, K., & Vansteelandt, S. (2026). Parameterizing the effect of a continuous treatment using average derivative effects. *Biometrika*, 113(2), asag012. https://doi.org/10.1093/biomet/asag012

Kennedy, E. H., Ma, Z., McHugh, M. D., & Small, D. S. (2017). Non-parametric methods for doubly robust estimation of continuous treatment effects. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(4), 1229–1245. https://doi.org/10.1111/rssb.12212

Neugebauer, R., & van der Laan, M. (2007). Nonparametric causal effects based on marginal structural models. *Journal of Statistical Planning and Inference*, 137(2), 419–434. https://doi.org/10.1016/j.jspi.2005.12.008

Semenova, V., & Chernozhukov, V. (2021). Debiased machine learning of conditional average treatment effects and other causal functions. *The Econometrics Journal*, 24(2), 264–289. https://doi.org/10.1093/ectj/utaa027

Takatsu, K., & Westling, T. (2025). Debiased inference for a covariate-adjusted regression function. *Journal of the Royal Statistical Society Series B: Statistical Methodology*, 87(1), 33–55. https://doi.org/10.1093/jrsssb/qkae041

Westling, T. (2022). Nonparametric tests of the causal null with nondiscrete exposures. *Journal of the American Statistical Association*, 117(539), 1551–1562. https://doi.org/10.1080/01621459.2020.1865168
