# Omnibus test p-values: study and fix for audit issue 11 (2026-11-12)

The omnibus test in `summary.effect_curve()` gave p-values that depended on the scale of the outcome, and some methods gave impossible values. This folder holds the studies used to choose the fix. The scripts run from the package root.

## Problems found

**Scale dependence.** With the outcome multiplied by `1e-8`, every test returned `p = 1`; at `1e-4`, the linearity test did. Two absolute tolerances caused this: `stat / n < 1e-10` declared the hypothesis exactly satisfied, and eigenvalues of the covariance below `1e-10` were dropped (with all of them dropped, the result was `NaN`).

**Impossible values from Imhof's method.** `CompQuadForm::imhof()` controls the absolute error of the probability (default `1e-6`), so below about `p = 1e-4` it returned negative values in 9% of the accuracy scenarios. These were silently floored at `1e-12`. Tightening its tolerances only reduced this to 7%.

**The `"saddlepoint"` method was Satterthwaite's with infinite `df`.** `survey::pFsum()` passes infinite `df` to `survey::pchisqsum()` without its `method` argument, so the Satterthwaite approximation was used. Its errors were identical to `"satterthwaite"` (up to a factor of 6 for unequal eigenvalues).

**Liu's method with finite `df`.** It is not designed for the negative weight that finite `df` introduces, and returned `p = 1.06` when the true value was .001.

## Accuracy of the methods (`accuracy.R`)

Upper-tail probabilities from each method were compared with an exact reference across 10 eigenvalue patterns (equal eigenvalues, one dominant eigenvalue, geometric decay, and eigenvalues from an actual test), exact p-values from .9 to `1e-14`, and `df` of `Inf`, 2500, and 30. The reference is Ruben's chi-squared mixture series in upper-tail form, whose terms are all positive, so it stays exact in the tails; with equal eigenvalues it reduces to exact chi-squared and F probabilities, which check the implementation.

Largest error, in log10 units (.05 is about 12% relative error), and share of impossible values:

| Method | p > .001 | .001 to 1e-5 | 1e-5 to 1e-9 | < 1e-9 | Impossible |
|---|---|---|---|---|---|
| Imhof | .004 | .126 | 1.0 | 6.8 | 8.8% |
| Davies | .001 | .031 | .88 | | 0.6% |
| Saddlepoint (fixed) | .030 | .035 | .046 | .054 | 0% |
| Satterthwaite | .40 | .77 | 2.7 | 6.3 | 0% |
| Liu | 3 | 4 | 8 | 14 | 0% |
| Imhof, saddlepoint below 1e-4 | .004 | .060 | .046 | .054 | 0% |
| Davies, saddlepoint below 1e-4 | .001 | .032 | .046 | .054 | 0% |

## Eigenvalues from actual tests (`realistic.R`)

For untruncated eigenvalues from flat, linear, and AMEF tests at grids of 21 and 51 points, compared with a Monte Carlo reference of `1e7` draws, Imhof and Davies agreed with the reference within its standard error (relative error at most 1.5% at `p = .001`); the saddlepoint approximation was within 6.5%. These tests have only two or three effective eigenvalues, because the ADRF of a cubic-polynomial model varies in few directions.

## The fix

- `.omnibus_p_value()` judges every tolerance relative to the largest eigenvalue of the weighted covariance of the estimates themselves, so results do not depend on the units of the outcome. When the tested quantity has no variance at that scale, the p-value is 1 if the estimates satisfy the hypothesis to numerical precision and 0 otherwise. P-values are restricted to [0, 1].
- `.quad_form_tail()` uses Imhof's or Davies's method where they are accurate and switches to the saddlepoint approximation for p-values below `1e-4`, invalid values, or Davies error codes.
- `.saddlepoint_tail()` implements the saddlepoint approximation of Kuonen (1999) in *adrftools*, in the Barndorff-Nielsen form used by `survey:::saddle()`, with the finite-`df` term entered once with its multiplicity. It no longer requires *survey*, and the default method without *CompQuadForm* is now `"saddlepoint"` rather than `"sim"`.

## Calibration under known truth (`calibration.R`)

Data with a known ADRF were generated (n = 500, outcome model `poly(A, 3) * (X1 + X2)`), so each tested null hypothesis was true: a flat ADRF tested for flatness (also with the outcome multiplied by `1e-6`), a linear ADRF tested for linearity, a flat ADRF for a binary outcome with and without the logit transformation, a reference curve tested against 0, and an AMEF tested against 0. There were 400 replicates per scenario.

- No method gave an impossible value in any of the 11,200 p-values.
- The p-values for the outcome multiplied by `1e-6` matched the unscaled ones to within `3e-12` in all 1,600 pairs.
- The four methods agreed with each other in every scenario (rejection rates within .005).
- **The tests were anti-conservative**, rejecting 7.5% to 9.75% of the time at the .05 level (Monte Carlo SE about .011). Because all methods agreed, this does not come from the p-value computation. A follow-up with the flat scenario (300 replicates) traced it to the variance estimate: 9.0% with the unconditional variance and with conditional HC0, 4.0% with conditional HC3, and 7.7% with the unconditional variance at n = 2000. The unconditional variance shares the small-sample bias of HC0 sandwich estimators.
- **Correction (2026-11-12):** the outcome model in `calibration.R` omits the `X1:X2` term of the data-generating process, so the ADRF it estimates is not exactly flat, and part of the excess (in particular the 7.7% at n = 2000) came from a false null. With the term included, the unconditional variance rejects 9.0%, 7.0%, 6.0%, and 5.2% at n = 250, 500, 1000, and 2000. See `_dev/omnibus-size-study/README.md`, which studies the remaining small-sample bias and its fixes.
- **Update (2026-11-16):** the remaining excess was the downward bias of the HC0-type variance for the units with high leverage at the extremes of the treatment. The unconditional variance now uses the HC2 adjustment by default (CR2 with `cluster`), which brings this scenario to within about a point of the nominal level from n = 250; with a flexible logistic model or with errors largest at the treatment extremes, 6.6% to 8.0% remain in moderate samples (`_dev/omnibus-size-study/README.md`).

## Reference

Kuonen, D. (1999). Saddlepoint approximations for distributions of quadratic forms in normal variables. *Biometrika*, 86(4), 929–935. https://doi.org/10.1093/biomet/86.4.929
