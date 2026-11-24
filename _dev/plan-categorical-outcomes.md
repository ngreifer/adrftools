# Plan: categorical outcomes (2026-11-24)

Adds support for models of multi-category outcomes (`nnet::multinom()`, `MASS::polr()`, `WeightIt::multinom_weightit()`, `WeightIt::ordinal_weightit()` on the probability scale, *brms* `categorical()` and `cumulative()`), which `adrf()` now rejects with an error raised in `.get_predict()` (commit c24e3e2). Implements the "Multi-category outcomes" item in `_dev/feature-ideas.md`.

**Recommendation:** estimate the curves for all categories at once, with the outcome category stored as an extra column of `.by_grid`. A separate argument to pick one category is unnecessary once this exists, because `subset` on the effect curve, `summary()`, and `plot()` already selects it.

**Implemented (2026-11-24), approach A, with these changes from the plan:**

- Categories are selected with a new `outcome` argument, not `subset`: in `adrf()` and `shift_curve()` (only those categories are computed), in the effect curve, `summary.effect_curve()`, `plot.effect_curve()`, `summary.curve_projection()`, and `ade()`. `subset` refers only to the `by` variables. `process_subset_by_grid()` handles both and returns one logical per block (subgroup, contrast, or category).
- The category column of `.by_grid` is named `outcome`, or `.outcome` (more dots as needed) when a `by` variable already has that name; the `.outcome` attribute holds the name. Its values are a factor in the order of the categories.
- `type` now defaults to `NULL`: `"response"` when `get_predict()` accepts it, and otherwise the type `marginaleffects::predictions()` uses by default (`"probs"` for `polr()`), found by `.default_type()` on one row of the data. `type = "mean"` and `"link"` for `ordinal_weightit()` give a single curve as before.
- Curve contrasts of categorical models keep the category of each contrast as the only column of `.by_grid`; contrast curves of other models now have `.by_grid = NULL` (it was the parent's grid, unused). Contrast labels in plots read `lo: [Male = 1] - [Male = 0]`.
- The y-axis of plots reads `P[Y = k]` for categorical models. With `by` (or for a curve contrast), each category gets its own panel (`facet_wrap()`, with a shared y-axis) and colors distinguish the subgroups or contrasts, so a subgroup has the same color in every panel; without `by`, colors distinguish the categories in one panel.
- **A bug found on the way (2026-11-24):** the unconditional variance of `polr()` curves was 20 to 30 times too large, because `marginaleffects::get_coef()` puts the thresholds first and `sandwich::estfun()` puts them last, so the gradient and the influence of the coefficients were paired in the wrong order. `.align_to_coef()` now reorders the influence columns by name to match the coefficients (and those of the CR2 influence). Among `betareg`, `gam`, *fixest*, `glm_weightit()`, and `glm()` fits, the orders already agree (`gam()` influence is unnamed and left as is). `polr()` could not be used in 0.1.0, so there is no NEWS item. The hand-computed variance test now covers `polr()` as well as `multinom_weightit()`, and fails with the fix removed.
- `.adrf_fit_one()`'s grid argument is now `grid_values`, because a `values` argument passed through `...` (needed by `predict.ordinal_weightit()` with `type = "mean"` and non-numeric categories) collided with it. The collision was introduced in this development cycle (0.1.0 passed `values` through), so there is no NEWS item.

**Verification:** estimates equal `marginaleffects::avg_predictions(by = "group")` exactly and conditional SEs to 1e-7 for `polr()`, `multinom_weightit()`, and `ordinal_weightit()`; the unconditional variance with `by` equals a hand-built influence function (the package's `.coef_influence()` with an independent finite-difference gradient and covariate-sampling term) to a relative 3e-6; *brms* `categorical()` posterior draws equal `posterior_epred()` per category. Tests are in `tests/testthat/test-categorical-outcomes.R` (63 expectations; the *brms* test is opt-in with `ADRFTOOLS_TEST_BRMS=true`); the full suite passes.

**Not done:** the second stage (linear combinations of categories); an analytic gradient for multinomial models. `nnet::multinom()` works but is not tested, because *nnet* is not in Suggests. `rstanarm::stan_polr()` is not supported: `marginaleffects::get_predict()` returns a single dummy group with estimates of 0 for it.

## Estimand

For an outcome with categories $1, \dots, K$, the ADRF for category $k$ is

$$\psi_k(a) = E[P(Y = k \mid A = a, X)],$$

averaged over $X$ in the target population (within subgroups defined by `by`). At every $a$, $\sum_k \psi_k(a) = 1$. For shift curves, $a$ is replaced by $g(A, \delta)$ as in `_dev/plan-shift-effects.md`; `shift_curve()` gets this for free through `.effect_curve_internal()`.

## What the models give now

Checked on `nhanes3lead` with a 3-category version of `Math` (`$TMPDIR/probe.R`, 2026-11-24):

- `marginaleffects::get_predict()` returns a long data frame with a `group` column (one block of rows per category) for all four models; its default `type` is `"probs"` for each. `nnet::multinom()` also accepts `type = "response"`; `polr()` does not.
- `marginaleffects::get_coef()` returns a single named vector (8 coefficients for the multinomial models, 5 for the ordinal ones), and `set_coef()` works, so the numeric gradient in `.adrf_fit_one()` applies unchanged.
- `.coef_influence()` already returns a 2521 x 8 (multinomial) or 2521 x 5 (ordinal) matrix, so the unconditional variance needs nothing new on the coefficient side. HC2 for the WeightIt models is in `_dev/TASKS.md` (2026-11-18).
- `process_family()` returns a logit link for all four.
- `.get_model_type()` falls back to `"default"` for these models (the prediction vector is $K$ times longer than the model-matrix one), so the fast path never applies.

Timing: for `multinom(M3 ~ ns(logBLL, 3) * (PIR + Male + Age + Race))` (56 coefficients) on the default grid (2521 units x 51 values = 128,571 rows), one `get_predict()` call takes 0.18 s (0.05 s for `predict.multinom()` directly). The central-difference gradient needs $2p = 112$ calls, about 20 s, for either the conditional or unconditional variance.

## Approach A: all categories, as a layer of `by`

### Why it is cheap

Every effect curve stores its estimates as one vector ordered by row of `.by_grid`, then grid value. `.make_fun()`, `add_est_labels()`, `process_subset_by_grid()`, `summary.effect_curve()`, `summary.curve_est()`, `plot.effect_curve()`, `print.effect_curve()`, `amef()`, `reference_curve()`, `point_contrast()`, and `curve_projection()` all work block by block over `fnrow(.by_grid)` (`get_n_by()`). If the category is a column of `.by_grid`, all of them work without change: estimates are labeled by category, `f(logBLL = 1, subset = category == "hi")` gives one category, `plot()` draws one line per category, `amef()` gives each category's derivative curve, and `summary()` tests each category's curve.

The joint covariance across categories is singular (the probabilities sum to 1). The omnibus tests use the eigenvalues of the weighted covariance and never invert it, and the correlation matrix for simultaneous bands is already prepared for `mvtnorm` when it is rank deficient, so neither breaks; both still need a test.

### Changes

In `R/get_estimator.r`:

- `.get_predict()` returns the predictions as a matrix with one column per category (rows in the order of `newdata`, matched by `rowid` rather than assumed from the row order of `get_predict()`), with the category labels as column names, in place of the error. Posterior draws get the same layout (rows by category, then unit). Single-prediction models still return a vector.

In `R/adrf.R`, `.adrf_fit_one()`:

- `est_fun()` returns a vector of length $K \cdot n_i \cdot n$, category-major. `cross_id` is built over (category, subgroup, value) and `wi` is repeated $K$ times, so `.fmean()` gives estimates in the order of the expanded `.by_grid`.
- Conditional variance: `gr` from `.gradient()` gains rows automatically; `quad_mult(gr, V)` is unchanged.
- Unconditional variance: `S <- tcrossprod(infl, gr)` is unchanged. The covariate-sampling loop (`for (g in seq_len(n_by))`) runs over (category, subgroup) pairs, with `P` an $n_i \times (n K)$ matrix. Each category uses the same units as its subgroup, since categories do not partition units.
- Bootstrap: `wi_boot` is repeated $K$ times.
- Bayesian: each grid value's draws have $K \cdot n_i$ columns and are averaged within (category, subgroup).
- $K$ is found once from the first prediction call; a refit (bootstrap replicate) that drops a category gives an error.

In `R/adrf.R`, `.effect_curve_internal()`:

- `by_grid` is expanded to (category x subgroups), with the category as the first column, so that the block for one category is laid out like an ordinary `by` ADRF.
- `.treat_density` is repeated for each category (the treatment distribution does not depend on the outcome category).
- With multiply imputed data, every imputation must give the same categories.
- A new attribute, `.category`, holds the name of the category column, for `curve_contrast()` and labels.
- `.family` is `binomial()` for category probabilities, whatever the model's link (a `polr()` fit with a probit or cloglog link still gives probabilities), so `transform = TRUE` uses the logit.
- When `type` is not supplied, use `"probs"` for models whose `get_predict()` rejects `"response"` (or, more simply, try `"response"` and retry with `"probs"`).

In `R/curve_contrast.R`:

- Contrasts are between `by` subgroups within each category, never between categories; without `by`, the error says that categories are not contrasted. The pooled treatment density uses the two subgroups' densities, as now.

Elsewhere:

- `plot.effect_curve()`: the legend title is "Category" when the only grouping is the category; with `by` as well, colors stay one per (category, subgroup), and faceting by category is worth considering.
- Messages that say `subset` requires `by` (`process_subset_by_grid()`) also allow the category.
- `adrf()` documentation: a section on categorical outcomes; NEWS item replacing the error item.

### Tests

- `nnet::multinom()`, `MASS::polr()`, `multinom_weightit()`, `ordinal_weightit()`: estimates and conditional SEs equal to `marginaleffects::avg_predictions(fit, variables = list(logBLL = values))` with `by = "group"`.
- Probabilities sum to 1 at every grid value, within every `by` subgroup.
- Unconditional SEs match a hand-computed influence function and a fractional weighted bootstrap (in `_dev/variance-verification/`, as for shift curves).
- `by` x category, `subset` on the curve, `summary()` (per-category tests and simultaneous bands across all categories), `amef()`, `curve_contrast()` within category, `plot()`.
- Multiply imputed data (*mice*), bootstrap, and *brms* `categorical()` if installed.
- `ordinal_weightit()` with `type = "link"` is unchanged.

## Approach B: one category chosen by the user

A `category` argument to `adrf()` and `shift_curve()`, passed to `.get_predict()`, which keeps the rows of `get_predict()` for that category (and the matching rows of the posterior draws). Nothing downstream changes. About 20 lines in `R/`, plus argument plumbing through `adrf.default()`, `adrf.mira()`, `shift_curve()`, and `.effect_curve_internal()`, documentation, and tests.

Its limits:

- $K$ calls for $K$ curves, each recomputing the gradient (about 20 s each in the timing above) or the bootstrap.
- No joint covariance across categories, so no sums or comparisons of categories (the second stage below is impossible).
- Separate bootstrap calls use different resamples.
- Once A exists, the `category` argument duplicates `subset` and becomes API to deprecate.

## Optional second stage: linear combinations of categories

For ordinal outcomes, the expected score $\sum_k s_k \psi_k(a)$ and the cumulative probabilities $P(Y \ge k)$ are linear in the category curves. A modifier (name to decide; e.g., `category_combine()`) can form them with `.make_fun(x, .contrast_mat = M, .by_grid = <by_grid without the category>)`, which already transforms estimates, covariance, bootstrap replicates, and posterior draws. `.family` would be the identity for scores and `binomial()` for cumulative probabilities. About 60 lines; needs A.

## Assessment

A is worth doing. The core change is local to `.get_predict()`, `.adrf_fit_one()`, and the end of `.effect_curve_internal()`, plus `curve_contrast()`; everything else comes from treating the category as a subgroup. Most of the work is in testing and the variance verification. B is cheap but is a strict subset of A, so it is worth doing only as a stopgap if A is put off.

The main risk is the one place where categories are not like subgroups: they share units. That matters for the covariate-sampling term, the treatment density, and `curve_contrast()`, all listed above; a grep for code that assumes subgroups partition units should precede the implementation.

Speed is the other cost. The numeric gradient dominates for models with many coefficients; an analytic gradient for multinomial models ($\partial p_k / \partial \beta_j = p_k(\mathbb{1}[k = j] - p_j)\, x$) and for `polr()` would remove it, but can wait.

## Open questions

- The name of the category column: `category` (clashes with a `by` variable of that name), the response name (awkward when the response is an expression such as `ordered(cut(Math, 4))`), or a user-supplied name.
- Whether `plot()` facets by category by default when `by` is also supplied.
- Whether to do the second stage, and its name.
- Whether the default `type` should be chosen automatically for `polr()` or left to the user with a clearer error.
