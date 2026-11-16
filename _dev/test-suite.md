# Test suite report (2026-11-12, updated 2026-11-16)

## Results

Line coverage is **95.7%** (`covr::package_coverage()` with `NOT_CRAN=true` on 2026-11-16; it was 91.8% on the first pass and 96.1% after the first audit). All 14 test files pass with 254 tests (each *patrick* case counts as a test) and 1,082 expectations: 0 failures, 0 warnings, and 1 skip (the opt-in *brms* test). The full suite runs in about a minute.

| File | Coverage |
|---|---|
| `R/aux_functions.R` | 90.4% |
| `R/get_estimator.r` | 93.0% |
| `R/print.effect_curve.R` | 96.4% |
| `R/utils.R` | 97.2% |
| `R/posterior_inference.R` | 97.3% |
| `R/summary.effect_curve.R` | 97.6% |
| `R/adrf.R` | 97.8% |
| `R/reference_curve.R` | 97.9% |
| `R/curve_projection.R` | 98.9% |
| `R/amef.R`, `R/curve_contrast.R`, `R/effect_curve.R`, `R/gradient.R`, `R/plot.effect_curve.R`, `R/point_contrast.R`, `R/summary.curve_est.R` | 100% |

Every issue in `_dev/audit.md` (issue 12 was retracted) and `_dev/audit-2.md` is fixed, and each fix has a regression test in the relevant file, so `test-known-issues.R` was removed. The remaining uncovered lines are mostly error branches for unusual inputs and model types, and the `interactive()` hint in `print.effect_curve()`. `R/simultaneous_ci_level.R`, which held an unused helper, moved to `_dev/`.

On CRAN, the Bayesian and BART files skip (`skip_on_cran()`), so CRAN coverage is lower; that matters only for the coverage number, not for correctness.

## Files

- `helper-adrftools.R`: loads `nhanes3lead` and defines the independent reference computations and accessors. `manual_adrf()` averages `predict()` over counterfactual data; `manual_uncond_vcov_lm()` builds the unconditional influence function by hand, with the HC2 adjustment by default or HC0 or HC1 on request; `ribbon_data()` extracts a plot's confidence band.
- `test-adrf.R`: point estimates for `lm`, `glm`, negative binomial, `betareg`, *fixest*, splines, offsets, matrix responses, models fit without `data`, and aliased coefficients (for every variance, and on the *marginaleffects* path); conditional variance against `marginaleffects::avg_predictions()`; the unconditional variance with HC2 against the hand-built influence function, the leverages against `hatvalues()` (including zero weights), the HC1 fallback, and `na.exclude` against `na.omit`; cluster-robust variance with CR2 (against a hand computation with eigendecompositions, against `sandwich::vcovCL()` for linear models, and equal to HC2 with singleton clusters), two-way clustering, and the CR1 fallback; routing to the fast path (including character predictors with `subset`); the model matrix built one treatment value at a time, against predictions on the whole grid for formulas that depend on other rows or read the treatment without naming it; `range`, `n`, `subset` (factors, row names, row numbers, empty subsets), `by` (including missing values), `wts` (validation, zero weights), reserved treatment names, the treatment-only shortcut, and input validation.
- `test-effect_curve.R`: evaluating curves at and between grid points (exact for cubic ADRFs), values just outside the grid (moved to the end points) and farther outside (errors), invariance to the scale of the treatment, `coef()`/`vcov()` labels, subgroup `subset`, and printing.
- `test-summary.curve_est.R`: Wald intervals and p-values against hand computations, simultaneous against pointwise, transformations (equivalent specifications of the logit and log links, `sqrt` and `exp`, and errors for other functions), errors for a null value or an explicit transformation outside the transformation's domain, null values, zero-variance estimates, and non-integer degrees of freedom.
- `test-curve-modifiers.R`: `amef()` against analytic derivatives (and of a reference curve), `reference_curve()`, `curve_contrast()` (including the error for one subgroup), `point_contrast()` (including of a `summary()`), and `subset = TRUE`, each against hand-built contrast matrices.
- `test-curve_projection.R`: projection coefficients and standard errors against a numerical Jacobian, named and formula models, transformed projections, large-scale treatments, subgroup subsets, `anova()` (F equals the squared t for one df; errors for different transformations, `vcov = "none"`, and extra terms that cannot be identified; long formulas), and errors for projections with too many or linearly dependent terms.
- `test-summary.effect_curve.R`: the omnibus test across all six p-value methods (agreement with Imhof), a direct Imhof computation, hypotheses by name, formula, and number, subgroups (and errors for `subset` without `by`), and transformations.
- `test-plot.R`: plotted values and bands against `summary()` (the same simultaneous critical value under the same seed), null lines, subgroups, contrasts, projections, and transformations, using `ggplot2::layer_data()`.
- `test-bootstrap.R`: *fwb* bootstrap estimates, percentile, Wald, and BC intervals, every curve modifier with bootstrap draws, replicates that alias a coefficient (against refitting with `predict()`), projections of draws (against per-draw fits, with and without a transformation), projection summaries with `subset`, draws outside a transformation's domain (a warning for the model's link, an error for an explicit transformation), and zero-variance p-values.
- `test-integration-survey.R`: *svyglm* with survey weights, unconditional variance against the equivalent weighted `lm()` (with the latter's HC2 adjustment switched off, since survey designs are not adjusted), conditional variance against `vcov(fit)`, and logistic models with `by`.
- `test-integration-weightit.R`: `glm_weightit()` and `lm_weightit()` with estimated weights and with sampling weights, `vcov = "HC0"` fits, whose unconditional variance treats the weights as fixed and matches `lm(weights = )`, and `ordinal_weightit()` on the link scale, which has no hat values and so uses HC1 with a message.
- `test-integration-mi.R`: *mice* `mira` objects checked against Rubin's-rules pooling of single-dataset ADRFs by hand (`pool_by_hand()`), covering `glm` on the response and link scales, `subset` and `by`, `vcov = "none"`, custom `vcov` and `cluster`, *betareg*, `wts`, variables outside the model, the Reiter degrees of freedom, stacked posterior draws for *rstanarm*, and *MatchThem* `mimira` objects.
- `test-integration-bayes.R`: *rstanarm* draws checked against `posterior_epred()` draw by draw, percentile and Wald inference (with p-values that invert the intervals), curve modifiers, and plots. A *brms* test runs when `ADRFTOOLS_TEST_BRMS=true` (it compiles a Stan model; it passed on 2026-11-12).
- `test-integration-bart.R`: *bartisan* (Gaussian and binomial), *dbarts*, and *stan4bart*, each checked against the package's own posterior predictions.
- `test-internals.R`: kernel and local-polynomial weights, bandwidths, numerical gradients, simulated p-values for the omnibus test (against exact chi-squared and F p-values), the exact sup-t bands of Montiel Olea and Plagborg-Møller (against bisection) and posterior p-values that invert them, simultaneous Wald critical values and p-values (against the Sidak formula), Rubin's rules, formatting, and small utilities.

## Parameterized tests

Eleven tests use `patrick::with_parameters_test_that()`, each nested in an outer `test_that()` that calls `skip_if_not_installed("patrick")` after any plain checks: the kernel types in `test-internals.R`, the grid sizes in `test-effect_curve.R` (twice), the equivalent `transform` specifications in `test-summary.curve_est.R` (three, for the logit, the log, and `sqrt` and `exp`), the p-value methods in `test-summary.effect_curve.R`, the named projection models in `test-curve_projection.R`, and the `vcov` abbreviations, `subset` forms, and formulas that challenge building the model matrix by treatment value in `test-adrf.R`. Eight of them came from the conversion of existing tests on 2026-11-12, described next; the others were written that way.

The conversion was checked as follows. Old and new versions of each file gave identical reporter totals, except `test-adrf.R` (104 vs. 103 passing), where the row-name `subset` case now also checks the covariance. Four planted bugs (a family passed to `transform` treated as identity, unnormalized triangular kernel weights, a shifted Liu p-value, and row-name subsets dropping a unit) produced the same failures in both versions, plus one more in the new `test-adrf.R` from that added check. With *patrick* reported as missing, only the eight parameterized tests skipped. No comparison lost its exactness.

## Conventions followed

- Numeric comparisons of separately computed quantities use `expect_equal()` with an explicit tolerance (1e-8 to 1e-10 for exact identities, looser only for Monte Carlo or numerical-derivative checks).
- Every correctness test compares against an independent computation (`predict()`, *marginaleffects*, *sandwich*, the posterior predictions of the fitting package, or algebra), not against stored adrftools output.
- Tests that need a suggested package call `skip_if_not_installed()`. Bayesian and BART tests use `skip_on_cran()` and small numbers of draws.

## Running

```r
devtools::test(filter = "adrf")                 # one file
devtools::test()                                # everything except Bayesian/BART on CRAN settings
Sys.setenv(NOT_CRAN = "true"); devtools::test() # everything
```

