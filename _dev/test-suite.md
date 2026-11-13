# Test suite report (2026-11-12)

## Results

Line coverage is **96.1%** (`covr::package_coverage()` with `NOT_CRAN=true`; it was 91.8% on the first pass). All 14 test files pass with 220 tests (each *patrick* case counts as a test) and 921 expectations: 0 failures, 0 warnings, and 1 skip (the opt-in *brms* test). The full suite runs in about 45 seconds.

| File | Coverage |
|---|---|
| `R/print.effect_curve.R` | 88.7% |
| `R/aux_functions.R` | 90.8% |
| `R/reference_curve.R` | 95.7% |
| `R/get_estimator.r` | 96.3% |
| `R/simultaneous_ci_level.R` | 96.6% |
| `R/adrf.R` | 98.1% |
| `R/summary.effect_curve.R` | 98.3% |
| `R/utils.R` | 98.5% |
| `R/point_contrast.R` | 98.8% |
| `R/plot.effect_curve.R`, `R/summary.curve_est.R`, `R/curve_projection.R` | 99.1% |
| `R/amef.R`, `R/curve_contrast.R`, `R/effect_curve.R`, `R/gradient.R`, `R/posterior_inference.R` | 100% |

Every issue in `_dev/audit.md` is now fixed (issue 12 was retracted), and each fix has a regression test in the relevant file, so `test-known-issues.R` was removed. The remaining uncovered lines are mostly error branches for unusual inputs and the `interactive()` hint in `print.effect_curve()`.

On CRAN, the Bayesian and BART files skip (`skip_on_cran()`), so CRAN coverage is lower; that matters only for the coverage number, not for correctness.

## Files

- `helper-adrftools.R`: loads `nhanes3lead` and defines the independent reference computations. `manual_adrf()` averages `predict()` over counterfactual data; `manual_uncond_vcov_lm()` builds the unconditional influence function by hand.
- `test-adrf.R`: point estimates for `lm`, `glm`, negative binomial, `betareg`, splines, and aliased coefficients; conditional variance against `marginaleffects::avg_predictions()`; unconditional variance against the hand-built influence function; cluster-robust variance (one- and two-way); `range`, `n`, `subset` (including errors for factors, unknown row names, and invalid row numbers), `by`, `wts`, the treatment-only shortcut, missing data, and input validation.
- `test-effect_curve.R`: evaluating curves at and between grid points (exact for cubic ADRFs), `coef()`/`vcov()` labels, subgroup `subset`, and printing.
- `test-summary.curve_est.R`: Wald intervals and p-values against hand computations, simultaneous against pointwise, transformations (logit, log, sqrt, exp, user functions, numeric inverse, and the fallback to the identity when a transformation cannot be reversed), null values, and zero-variance estimates.
- `test-curve-modifiers.R`: `amef()` against analytic derivatives, `reference_curve()`, `curve_contrast()`, and `point_contrast()`, each against hand-built contrast matrices.
- `test-curve_projection.R`: projection coefficients and variance, named and formula models, transformed projections, subgroup subsets, and `anova()` (F equals the squared t for one df).
- `test-summary.effect_curve.R`: the omnibus test across all six p-value methods (agreement with Imhof), a direct Imhof computation, hypotheses by name, formula, and number, subgroups, and transformations.
- `test-plot.R`: plotted values and bands against `summary()`, null lines, subgroups, contrasts, projections, and transformations, using `ggplot2::layer_data()`.
- `test-bootstrap.R`: *fwb* bootstrap estimates, percentile, Wald, and BC intervals, and every curve modifier with bootstrap draws.
- `test-integration-survey.R`: *svyglm* with survey weights, unconditional variance against the equivalent weighted `lm()`, conditional variance against `vcov(fit)`, and logistic models with `by`.
- `test-integration-weightit.R`: `glm_weightit()` and `lm_weightit()` with estimated weights and with sampling weights, and `vcov = "HC0"` fits, whose unconditional variance treats the weights as fixed and matches `lm(weights = )`.
- `test-integration-mi.R`: *mice* `mira` objects checked against Rubin's-rules pooling of single-dataset ADRFs by hand (`pool_by_hand()`), covering `glm` on the response and link scales, `subset` and `by`, `vcov = "none"`, custom `vcov` and `cluster`, *betareg*, `wts`, variables outside the model, stacked posterior draws for *rstanarm*, and *MatchThem* `mimira` objects.
- `test-integration-bayes.R`: *rstanarm* draws checked against `posterior_epred()` draw by draw, percentile and Wald inference, curve modifiers, and plots. A *brms* test runs when `ADRFTOOLS_TEST_BRMS=true` (it compiles a Stan model; it passed on 2026-11-12).
- `test-integration-bart.R`: *bartisan* (Gaussian and binomial), *dbarts*, and *stan4bart*, each checked against the package's own posterior predictions.
- `test-internals.R`: kernel and local-polynomial weights, bandwidths, numerical gradients, `make_inverse()`, simulated p-values, `simultaneous_conf_level()`, Rubin's rules, formatting, and small utilities.

## Parameterized tests

Eight tests use `patrick::with_parameters_test_that()`, each nested in an outer `test_that()` that calls `skip_if_not_installed("patrick")` after any plain checks: the kernel types in `test-internals.R`, the grid sizes in `test-effect_curve.R` (twice), the equivalent `transform` specifications in `test-summary.curve_est.R`, the p-value methods in `test-summary.effect_curve.R`, the named projection models in `test-curve_projection.R`, and the `vcov` abbreviations and `subset` forms in `test-adrf.R`.

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

