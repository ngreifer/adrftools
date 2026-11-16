# Tasks

## To Do

- R38 (rest): build the model matrix one treatment value at a time to cut peak memory (`_dev/audit-2.md`).
- Implement the HC2 adjustment to the unconditional variance (M19; see `_dev/omnibus-size-study/README.md`).
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
