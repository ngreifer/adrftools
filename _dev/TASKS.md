# Tasks

## To Do

- M19 follow-ups: the omnibus test still rejects 6.6% to 8.0% at .05 with a flexible logistic model (n <= 1000) or errors largest at the treatment extremes (n <= 500), from the variability of the HC2 estimate; consider a degrees-of-freedom correction or another reference distribution. The size of CR2 has not been checked.
- Decide whether the default `hypothesis` for reference curves should be `"flat"` (identical to the flat ADRF test for any reference) rather than 0, and correct the docs that call the reference-curve and AMEF tests equivalent to the flat test (`?reference_curve` example, vignette); see `_dev/omnibus-size-study/README.md`, equivalent hypotheses.
- Revisit the default of `integrate` (now `"uniform"`); the test-choice study favors `"empirical"` (`_dev/omnibus-test-choice/README.md`).
- `shift_curve()` and `ade()` (branch `shift-curve`): settle the open questions in `_dev/plan-shift-effects.md` (default shift size, the name `shift_fun`, a default reference of 0 for shift curves).
- Decide which items in `_dev/feature-ideas.md` to pursue.
- Multivariate treatments (`_dev/adrf_mv.R`), to resume later.

## Done

- 2026-11-12: Audit of `R/` (`_dev/audit.md`); 31 confirmed issues plus minor notes.
- 2026-11-12: Test suite in `tests/testthat/` (15 files, 141 tests, 91.8% line coverage); report in `_dev/test-suite.md`.
- 2026-11-12: Feature and improvement list (`_dev/feature-ideas.md`).
- 2026-11-12: Fixed audit issues 1, 2, 3, 6, 8, 10, 19, and 24, with NEWS items and regression tests.
- 2026-11-12: Added *patrick* to Suggests and parameterized eight tests.
- 2026-11-12: Fixed issue 4 and verified the unconditional variance against *stdReg*, *geex*, *survey* replicates, and the bootstrap (`_dev/variance-theory.md`, `_dev/variance-verification/`); retracted issue 12.
- 2026-11-12: Fixed issues 5, 7, 13, and 26, plus replicate-weight `svyglm` target weights; struck through completed audit items.
- 2026-11-12: Fixed issue 27 and two related `cluster` bugs (formulas with transformations, data with missing values); confirmed the subset behavior of omnibus and simultaneous tests.
- 2026-11-12: Fixed issue 31; revised the bootstrap `cluster` error to suggest `fwb::vcovFWB()` for multiway clustering; numbered the minor audit issues M1–M18.
- 2026-11-12: Unified the single-dataset and multiply imputed paths in `adrf()`, fixing issues 14–17 and minor issues M4 and M11.
- 2026-11-12: Fixed issues 18, 20, 21, 22, and 23.
- 2026-11-12: Fixed issues 9, 11, 25, 28, 29, and 30 (all audit issues now resolved); routing by predictions fixed silently wrong estimates for models with offsets; omnibus p-value study in `_dev/omnibus-pvalue-study/`.
- 2026-11-12: Fixed minor issues M1, M5, M7, M8, and M16 (M2 fixed by the maintainer), with regression tests and NEWS items.
- 2026-11-12: Fixed minor issues M3, M6 (an ill-conditioned local polynomial solve), and M12 (Reiter degrees of freedom), plus non-integer `df` in simultaneous inference.
- 2026-11-12: Investigated M19 (`_dev/omnibus-size-study/`): part of it was a misspecified model in the calibration study; the rest is HC0's leverage bias, which HC2 removes. Found a separate bug in `.data_matches_model()` (row-count mismatch accepted as a match).
- 2026-11-12: Fixed M21 (a same-named dataset in the formula's environment used in place of the model's data).
- 2026-11-13: Fixed M15 (`wts` docs) and M17 (exact implementation of Montiel Olea and Plagborg-Møller's Algorithm 2); marked M13, M14, M18, and M20 resolved.
- 2026-11-13: Second audit (`_dev/audit-2.md`): 46 findings (R1–R46), 16 silently wrong results, 15 errors, 8 performance, 7 documentation.
- 2026-11-13: Fixed R1–R3, R5, R6, R8, R10–R13, R16, R20, R23, R27, R29, R32, R35, and R45 from the second audit.
- 2026-11-13: Fixed R7, R9, R15, R17–R19, R22, and R30 from the second audit; resolved R4 and R26 by restricting `transform` to known transformations.
- 2026-11-13: Fixed R21, R28, R31, R36, R37, R40, and R41 from the second audit; `anova()` errors when the extra terms cannot be identified.
- 2026-11-13: Fixed R14: effect curves and references error outside the grid, with values within 0.1% of its width moved to the end points.
- 2026-11-13: Fixed R24 and R25: a null value or an explicit transformation outside the domain gives an error; the outcome model's link falls back to the identity with a warning.
- 2026-11-13: Fixed R33 and R34: shared helpers for simultaneous Wald critical values and p-values.
- 2026-11-13: Fixed R39, R42–R44, R46, and most of R38 from the second audit.
- 2026-11-14: Declared the packages used in tests in Suggests; R CMD check is clean (0 errors, 0 warnings, 0 notes).
- 2026-11-16: M19 study of the source of the size error and of a Wald test (`_dev/omnibus-size-study/run-wald.R`, README).
- 2026-11-16: Fixed M19: HC2-adjusted unconditional variance by default; fixed `NA` unconditional SEs with `na.exclude`.
- 2026-11-16: CR2 for the clustered unconditional variance of `lm()` and `glm()` fits; note on `hatvalues()` methods written to `WeightIt/_dev/hatvalues-ordinal-multinom.md`.
- 2026-11-16: Models without hat values fall back to HC1, and models without CR2 to CR1, each with a message.
- 2026-11-16: Checked the HC2 size with a logistic model and three heteroskedastic patterns (`_dev/omnibus-size-study/run-glm-het.R`, README).
- 2026-11-16: Updated the `_dev` notes for the variance changes: a section on HC2, HC1, CR2, and CR1 in `variance-theory.md`, current results and file descriptions in `test-suite.md` (254 tests, 1,082 expectations, 95.7% line coverage), and implemented ideas struck through in `feature-ideas.md`.
- 2026-11-16: Finished R38: `adrf()` evaluates the model frame on the whole grid and builds the model matrix one treatment value at a time from it (peak memory on the `?adrf` example 118 MB to 53 MB at `n = 51`, 447 MB to 169 MB at `n = 201`); the bootstrap stores the columns that are identical at every treatment value once, found by comparing them. Adversarial formulas (row-dependent variables, a variable that reads the treatment without naming it, and others) are tested against predictions on the whole grid.
- 2026-11-17: Compared omnibus tests of equivalent hypotheses (`_dev/omnibus-size-study/run-equiv.R`, README): all hold their size, but their power differs, and the test of a reference curve against 0 depends on the reference.
- 2026-11-16: Test-choice study (`_dev/omnibus-test-choice/`): with a fixed (oracle) covariance, every statistic holds its size; testing the ADRF itself with empirical treatment weighting had the best or near-best mean and worst-case power for both flatness and linearity.
- 2026-11-16: Added `integrate = c("uniform", "empirical")` to `summary.effect_curve()` and `curve_projection()`; `adrf()` stores a kernel density estimate of the treatment for each subgroup (`_dev/plan-density-weights.md`).
- 2026-11-16: Plan for shift-response curves and the average derivative effect (`_dev/plan-shift-effects.md`): `shift_curve()` generalizes the counterfactual treatment in `.adrf_fit_one()` to a function of the observed treatment, and `ade()` is the derivative of a shift curve at 0.
- 2026-11-16: Implemented `shift_curve()` and `ade()` on branch `shift-curve` (`_dev/plan-shift-effects.md`), with tests (`test-shift-curve.R`) and an independent check of the variance (`_dev/variance-verification/verify-shift.R`); fixed `amef()` returning `NA` with 101 or more grid points.
- 2026-11-17: New vignette `vignettes/shift-effects.Rmd` on shift interventions and the ADE; collapse rewrites of grouped operations in `.adrf_fit_one()` (covariate-sampling term of `S`), `.get_treat_density()`, `.zero_var()`, `.draw_zeta()`, and `.mvt_p_values()`, checked identical to the previous code end to end; reference values rounded in printed notes and plot labels; `shift_curve()` and `ade()` added to `_pkgdown.yml`.
