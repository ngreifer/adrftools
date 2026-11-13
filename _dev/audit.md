# adrftools audit (2026-11-12)

Audit of the package as of commit `4220e6d` plus the uncommitted working-tree changes present on 2026-11-12. Every numbered issue below was reproduced in R, not only found by reading. At the time of the audit, each issue had a skipped test in `tests/testthat/test-known-issues.R` asserting the correct behavior; every one of them failed on the code as audited.

None of the source files in `R/` were changed during the audit. All issues were fixed afterward on the same day, except issue 12, which was found not to be a bug. Their headings are struck through, each has a note saying what changed, and each fixed issue's test moved from `test-known-issues.R` to the relevant test file as a regression test. With every issue resolved, `test-known-issues.R` was removed. Minor issues are numbered M1–M20 at the end.

Line numbers refer to the files as they stand in the working tree.

## Silent wrong results

These return plausible numbers with no error or warning.

### ~~1. The treatment can be read from a transformed model-frame column~~

**Fixed.** `process_model_data()` now rebuilds the model frame from the recovered data and checks it against the stored one. If they differ, the user's `data` is used when supplied; otherwise an error asks for it. The raw treatment cannot be recovered from a `poly()` basis, so an error is the correct outcome when `data` is missing.

`process_model_data()` (`R/aux_functions.R:203`) uses `insight::get_predictors()`. When *insight* cannot find the original data (the model was fit inside a function, or with a formula created elsewhere), it falls back to the model frame. There, `poly(logBLL, 2)` is a two-column matrix, and `get_predictors()` returns its first column under the name `logBLL`.

Reproduction: `fitter <- function(f, d) lm(f, data = d)`, fit `Math ~ poly(logBLL, 2) + Male + PIR`, call `adrf()`. The printed range is `-0.0364 to 0.0399` (the orthogonal polynomial basis) instead of `-0.36 to 2.42`, and every estimate is wrong. With *rstanarm* in the same setting the column comes back as a matrix and `process_range()` fails with a length error.

Fix ideas: prefer `model[["data"]]` when it is a data frame, then `insight::get_data(source = "environment")`; check that `process_treat()` returns a plain numeric vector (no `dim`); and as a guard, rebuild the model matrix from the extracted data and compare it with `model.matrix(x)` before trusting it. Data modified after fitting has the same failure mode and the same guard catches it.

### ~~2. Bayesian ADRFs with `by` are ordered wrongly~~

**Fixed.** The draws are reordered with `.value_major_to_by_major()` in both the single-dataset and multiply imputed paths.

`R/adrf.R:261-270` (and `582-590` for MI) builds `.draws` by looping over treatment values and `rbind`ing the per-subgroup means. The columns come out value-major (`v1 g1, v1 g2, v2 g1, ...`), but everything downstream assumes subgroup-major (`g1 v1, g1 v2, ...`). The labels then pair estimates with the wrong subgroup and treatment value.

Reproduction: `stan_glm(Math ~ logBLL * Male + PIR)` with `by = ~Male`, `n = 3` gives `9.07 8.79 8.19 7.60 7.31 6.41`; the matching `lm()` gives `9.07 8.19 7.32 8.79 7.60 6.42`.

Fix: reorder the columns of `.draws` (e.g., build them in the `cross_id` order) before computing `.est`.

### ~~3. *svyglm* unconditional variance is scaled by the subset size~~

**Fixed.** The influence function is scaled by its own row count (the full sample size) in both paths.

`R/adrf.R:346-348` multiplies the stored or recomputed influence function by `ni`, the number of rows *after* subsetting. It should use the full sample size. Two consequences:

- With the treatment-only shortcut (`ni` set to 1, `R/adrf.R:204-211`), SEs are too small by a factor of about N: `svyglm(Math ~ logBLL)` gave SEs of `7e-05` against `0.18` from `vcov(fit)`.
- With `subset`, SEs shrink by the subset fraction: `0.114` against `0.220` for the equivalent weighted `lm()` with `subset = Male == 1`.

Fix: use `fnrow(S)` (full sample) or `length(x$y)`.

### ~~4. Weighted unconditional variance omits the unit weights~~

**Fixed.** The term is now `wts[b] * (phat - est) * N / sum(wts[b])` in both paths. Verified against *stdReg*, *geex*, and the bootstrap; see `_dev/variance-theory.md`.

`R/adrf.R:357-368` (and `628-639`) adds `(phat_i - est) * sum(wts) / sum(wts[b])` to the influence function. For a weighted mean, the term is `w_i * (phat_i - est) * N / sum(wts[b])`. The two agree only when all weights are 1. This affects every analysis with non-uniform `wts`, including *svyglm* and *WeightIt* models with sampling weights.

Reproduction: `lm(Math ~ logBLL + PIR)` with `wts = rexp(n)`; the SEs from `adrf()` (`0.1390 0.0855 0.0645`) match the formula in the code but not the correct influence function (`0.1398 0.0870 0.0667`), which a 2,000-replicate fractional weighted bootstrap also supports (`0.1413 0.0873 0.0657`). The gap grows with the variance of the weights; with `rexp(n)^2` weights the known-issue test fails clearly.

Fix: `S[b, k] <- S[b, k] + wts[b] * (phat_split[[k]] - .est[k]) * N / sum(wts[b])`, with `N = fnrow(S)`.

### ~~5. Bootstrapping refits the model with `update()`, which drops weights and needs the original environment~~

**Fixed.** `.make_boot_refit()` refits through the model's own call, on the rows used in the fit, with the model's weights times the bootstrap weights; predictions are averaged with `wts` times the bootstrap weights. Models fit with a `weightit` object have their weights re-estimated in each replicate (matching *WeightIt*'s own `vcov = "FWB"` exactly with the same seed). Confirmed beforehand: the weighted `lm()` case, `wts` used as model weights (bootstrap estimates `9.00 7.87 6.74` vs. `8.91 7.87 6.82`), and fitting inside a function all reproduced. Also found: bootstrapped `svyglm` SEs were about twice the design-based ones, so bootstrapping `svyglm` models is now an error pointing to the design-based variance.

`bootfun()` (`R/adrf.R:288-301`) calls `update(x, weights = wts * w)`. This has three problems.

- It replaces the model's own weights. For `lm(..., weights = w)`, the bootstrap point estimates are those of an unweighted fit (`8.92 8.40 7.89` vs. `9.01 8.45 7.90`). The same holds for *WeightIt* models and anything else fit with weights.
- `wts` are the *target-population* weights, not model weights. When the user supplies `wts`, the model is refit with them, which changes the outcome model.
- `update()` re-evaluates the original call in the caller's frame, so a model fit inside a function fails with `object 'd' not found`.

Fix: multiply the bootstrap weights into the existing model weights (e.g., `update(x, weights = model_weights * w, data = model_frame_with_vars)`), and use the bootstrap weights times `wts` only in the averaging step. Evaluating the refit with the model's own data (`insight::get_data()` result passed explicitly) removes the environment dependence.

### ~~6. Simultaneous posterior p-values go the wrong direction~~

**Fixed.** `posterior_ci()` and `posterior_p_value()` now use the single-pass sup-t method from *fwb* (`sup_t_entry_levels()`/`sup_t_level()`, adapted). Each simultaneous p-value is the proportion of draws whose entry level exceeds one minus the pointwise p-value. In a check with 4,000 draws, these p-values matched a brute-force inversion of the simultaneous band to within 0.0015. In 1,000 simulated flat-prior posteriors with every null true, the familywise error rate was 4.6% (pointwise: 20.3%). This also fixes issue 24.

`posterior_p_value()` (`R/posterior_inference.R:28-31`) passes each pointwise p-value through `.pointwise_p_to_simul_p()`. That function maps a *simultaneous* alpha to the *pointwise* alpha that achieves it, which is smaller. A simultaneous p-value needs the inverse map and should be at least as large as the pointwise one.

Reproduction: four correlated posterior draws with pointwise p-values `0.118 0.073 0.049 0.035` get "simultaneous" p-values `0.040 0.024 0.015 0.011`.

Fix: invert the map, i.e., find the smallest simultaneous alpha whose band excludes the null (`uniroot()` over alpha of `.pointwise_p_to_simul_p(alpha) - p`).

### ~~7. `summary.effect_curve(subset = )` uses the wrong bootstrap or posterior columns~~

**Fixed.** `.boot$t` and `.draws` are now subset along with the estimates. Confirmed beforehand: with `subset = Male == 1`, the bootstrap test gave `0.115` against `0.119` for the matching row of the unsubsetted test.

**On multiplicity.** `summary.effect_curve()` runs a separate omnibus test in each subgroup with no adjustment across subgroups, so restricting it to one subgroup should leave that subgroup's p-value unchanged; after the fix it does, for analytic and bootstrap inference. The pattern before the fix also rules out multiplicity: subsetting to the first subgroup (`Male == 0`) matched the full test exactly, and only the second did not, which is what misaligned columns produce. Where multiplicity does apply is `summary()` of a `curve_est` with `simultaneous = TRUE`: there, subsetting with `a(..., subset = )` leaves pointwise p-values and intervals unchanged but gives smaller simultaneous p-values and narrower bands, because fewer estimates are covered jointly. A test now checks both behaviors.

`R/summary.effect_curve.R:186-200` subsets `.est0` and `.vcov0` but not `.boot$t` or `.draws`. The column indices built later refer to the subset, so they pick out the first subgroup's draws.

Reproduction: bootstrap ADRF `by = ~Male`, `hypothesis = "linear"`; `subset = Male == 1` gives `0.123`, but the second row of the unsubsetted test is `0.155`.

Fix: subset `.boot$t` and `.draws` by `.s` alongside `.est0`.

### ~~8. Multiply imputed analyses silently ignore a supplied `vcov` or `cluster`~~

**Fixed.** The list element is now named `vcov`.

The final `lapply()` of `process_vcov_and_cluster_mi()` (`R/aux_functions.R:875-878`) returns `list(vcov[[.i]])` without the name `vcov`. `.effect_curve_internal_mi()` then reads `cl_list[["vcov"]]` as `NULL` and calls `get_vcov(xi, vcov = NULL)`. `vcov = "HC3"`, `vcov = "conditional"` with a cluster, and a user-supplied matrix all give the default model covariance.

Fix: `c(list(vcov = vcov[[.i]]), ...)`.

### ~~9. Interpolation is wrong when `n` is 2 or 3~~

**Fixed.** The local polynomial's degree is `min(p, length(pos) - 1)` per point, and the orthogonal basis's degree is capped by the number of unique points (its columns are nested, so the first `p_i + 1` span degree `p_i`). Interpolation is exact for a linear ADRF at every grid size; with two grid points it is linear interpolation.

`get_locpoly_w()` fits a local cubic (`p = 3`) whatever the grid size. With fewer than 4 grid points the fit is singular and `.to_psd()` returns a distorted solution. For a linear ADRF, interpolation errors are as large as 1.02 (`n = 2`) and 0.65 (`n = 3`); from `n = 4` up the error is below `1e-11`. `process_range()` allows `n = 2`.

Fix: use `p <- min(p, length(pos) - 1L)` per point, or require `n >= 4`.

### ~~10. Negative `subset` indices keep the units they should drop~~

**Fixed.** The index is now `s <- seq_row(data) %in% seq_row(data)[subset]`.

`process_subset()` (`R/aux_functions.R:61-63`) sets `s <- TRUE` and then `s[subset] <- FALSE`. With negative indices, `s[subset]` refers to every *other* unit, so `subset = -which(Male == 0)` keeps only `Male == 0`.

Fix: `s <- seq_row(data) %in% seq_row(data)[subset]`.

### ~~11. The omnibus test returns `p = 1` for outcomes on a small scale~~

**Fixed.** See `_dev/omnibus-pvalue-study/README.md`. Tolerances are now relative to the variance of the estimates (p-values with the outcome scaled by `1e-6` match the unscaled ones to `3e-12`), p-values are restricted to [0, 1], Imhof's and Davies's methods hand off to a saddlepoint approximation below `1e-4` (Imhof returned negative values there in 9% of accuracy scenarios), and the saddlepoint method is implemented in *adrftools*. The study also found that `method = "saddlepoint"` had been Satterthwaite's approximation with infinite `df` because of how `survey::pFsum()` passes arguments, and that the omnibus tests are anti-conservative at n = 500 because of the small-sample bias of HC0-type variances (see M19).

`R/summary.effect_curve.R:322` returns 1 when `stat / n < 1e-10`, an absolute threshold. Multiplying the outcome by `1e-6` changes the p-value from `< 1e-8` to `1`.

Fix: compare `stat` with `sum(diag(vcov.b))` or with the scale of the estimates, or drop the shortcut and let the zero-eigenvalue branch handle exact fits.

### ~~12. *WeightIt* models: unconditional variance ignores estimation of the weights~~

**Retracted; not a bug.** `glm_weightit` objects have *WeightIt* `estfun()` and `bread()` methods whose estimating functions include the weight model's contribution (`asympt = TRUE` by default), so `adrf()` already applies the M-estimation correction: `adrf(glm_weightit)` equals the delta method on `vcov(glm_weightit)`. The identical results described below arise because `estfun()` reads the attached `weightit` object regardless of the fit's `vcov` setting. A remaining design question is whether `adrf()` should pass `asympt = FALSE` when the fit used `vcov = "HC0"`. See `_dev/variance-theory.md`.

The unconditional variance uses `sandwich::estfun()` and `sandwich::bread()`, which for `glm_weightit` dispatch to the `glm` methods. The M-estimation correction that `glm_weightit()` applies to account for estimating the weights is lost. A `glm_weightit()` fit with the default M-estimation covariance and one with `vcov = "HC0"` give identical `adrf()` SEs, although `vcov()` differs between them. `vcov = "conditional"` uses the M-estimation covariance correctly.

Fix: use the stacked estimating functions *WeightIt* stores (the outcome model's `psi` and the weight model's components) to form the coefficient influence function.

### ~~13. *svyglm*: unconditional variance ignores the design~~

**Fixed.** `.svyglm_vcov()` divides the influence values by the sampling weights and totals them as design variables with `survey::svytotal()`, so *survey* applies the variance estimator for the design (clusters, strata, finite population corrections, calibration, PPS, or replicate weights). `adrf(svyglm)` now matches the hand-computed Binder linearization exactly and the design-based jackknife within 1% for stratified cluster, weights-only, replicate-weight, calibrated, and domain designs, and within subgroups. `cluster` is ignored with a warning for these models. Fixing this also exposed a related bug: `insight::get_weights()` returns `NULL` for `svyglm` models fit to replicate-weight designs, so their predictions were averaged without the survey weights; the weights now come from `prior.weights`.

`R/adrf.R:370-377` computes `crossprod(S / N)`, an i.i.d. sandwich, unless `cluster` is given. PSUs, strata, and finite-population corrections in the design are ignored. A one-stage cluster design with 100 PSUs gives the same SEs as `lm()` with no clustering.

Fix: pass the influence functions to `survey::svyrecvar()` with the design's `cluster`, `strata`, and `fpc` (or `survey::svytotal()` on them as a design variable). Replicate-weight designs need their own handling.

## Errors on supported inputs

### ~~14. Multiply imputed `glm` with `type = "link"` fails~~

**Fixed** by unifying the paths. The single-dataset and multiply imputed paths are now one function: `.effect_curve_internal()` builds a list of models (one, or one per imputed dataset) and processes data, weights, `vcov`, and `subset` per model with the single-dataset functions, forms the `by` subgroups across all models after subsetting, and computes each model's estimates with `.adrf_fit_one()` before pooling. The family and model type are now determined with `type` for every model. `adrf()` on a `mira` object now matches Rubin's-rules pooling of `adrf()` applied to each completed dataset by hand.

`R/adrf.R:464-466` compares `process_family(mod)` without `type` against the first model's family with `type`, so the check fails with "All models must have the same family". Separately, `.get_model_type()` at `R/adrf.R:471` is called without `type`, so even with the check fixed a `glm` would use the response-scale estimator for `type = "link"`.

### ~~15. Multiply imputed `subset` fails~~

**Fixed** by unifying the paths: `subset` is applied per model before the subgroups are formed. Matches pooling by hand, with and without `by`.

`R/adrf.R:531` takes `by_id.complete[.imp]` without applying the subset `s`, so `by_id` and the subsetted data have different lengths ("argument lengths differ").

### ~~16. Multiply imputed `vcov = "none"` and Bayesian models fail~~

**Fixed** by unifying the paths: covariance specifications are processed per model by `process_vcov_and_cluster()`, and posterior draws are stacked across imputed datasets. Bayesian results match the per-imputation draws stacked by hand.

`rep_with(list(vcov = "none"), models)` (`R/aux_functions.R:767`, also `791`) flattens the list: each element becomes `"none"` rather than `list(vcov = "none")`, and `cl_list.list[[1L]][["vcov"]]` errors with "subscript out of bounds". Once that is fixed, the Bayesian branch still fails at `R/adrf.R:673`, which assigns `est.list[[.i]] <- est` although `est` is never set in that branch.

Fix: `rep_with(list(list(vcov = ...)), models)` (or `rep(list(list(...)), length(models))`), and skip the `est.list` assignment for Bayesian models.

### ~~17. Multiply imputed models outside the `lm`/`glm` fast path fail~~

**Fixed** by unifying the paths: `.adrf_fit_one()` defines `beta0` for every non-Bayesian model. *betareg* fits match pooling by hand for unconditional and conditional variances.

`beta0` and `aliased` are defined only inside `if (!is_bayes && is_null(pred_fun))` (`R/adrf.R:562-574`), but the unconditional and conditional branches use them for every model. `betareg`, `glm` with `type = "link"`, and anything else using `marginaleffects::get_predict()` fail with "object 'beta0' not found". The single-dataset path defines them earlier and works.

### ~~18. `curve_contrast()` fails for bootstrap and posterior curves~~

**Fixed.** `.make_fun()` subsets the draws to `seq_len(ncol(.contrast_mat))`, the number of original estimates. The contrast draws equal the differences of the subgroup draws, for bootstrap and posterior curves.

`.make_fun()` (`R/effect_curve.R:89-109`) overwrites `.est` with the contrasts before subsetting the draws with `seq_along(.est)`. When the contrast matrix is not square (always, for `curve_contrast()`), the draws are cut to the wrong width and `tcrossprod()` fails with "non-conformable arguments". `amef()` and `reference_curve()` use square matrices, which is why they work.

Fix: record `n_est <- length(.est)` before applying the contrast.

### ~~19. `curve_projection()` fails for posterior curves~~

**Fixed.**

`R/curve_projection.R:210` calls `stats::vcov(.draws)` on a matrix; it should be `stats::cov(.draws)`.

### ~~20. `curve_projection(x, "flat")` fails for bootstrap and posterior curves~~

**Fixed.** The projected coefficients are coerced to a matrix with one column per coefficient. Each replicate of a flat projection equals the trapezoidal average of that replicate's curve, for bootstrap and posterior curves.

With one coefficient, `dapply(MARGIN = 1L, return = "matrix")` (`R/curve_projection.R:165-172` and `196-203`) returns a vector, so `vcov()`/`cov()` fail. Wrap the result in `matrix(..., ncol = length(proj_coefficients))`.

### ~~21. `vcov()` and projection summaries fail when `vcov = "none"`~~

**Fixed.** `vcov.curve_est()` and `vcov.curve_projection()` return `NULL` when there is no covariance matrix, and `summary.curve_projection()` then reports only the coefficients.

`vcov.curve_est()` and `vcov.curve_projection()` set `dimnames` on `NULL`. `summary.curve_projection()` calls `vcov()` unconditionally (`R/curve_projection.R:255`), so summarizing a projection of a `vcov = "none"` curve fails.

### ~~22. `point_contrast()` loses its attributes with subgroups~~

**Fixed.** `point_contrast()` saves its attributes before labeling and restores them afterward. Estimates, covariance, and labels match hand-built contrasts within each subgroup, with analytic, bootstrap, and posterior inference.

`add_est_labels()` merges with the `by` grid (`R/aux_functions.R:1173-1186`), and `merge()` drops the attributes that `point_contrast()` set just before (`R/point_contrast.R:72-110`). `summary()` then fails because `.vcov_type` is `NULL`. Set the attributes after `add_est_labels()`, as the `effect_curve` call does.

### ~~23. `summary()` of a single estimate on a non-identity link fails~~

**Fixed.** `diag()` is called with `nrow = length(x)` here and at the other three places with the same pattern (`plot()` and twice in `summary.effect_curve()`). A single estimate's interval equals its pointwise interval in a larger summary.

`diag(d_transform(est0))` (`R/summary.curve_est.R:133`) builds an identity matrix of size `est0` when `est0` has length 1. A logistic ADRF evaluated at one treatment value cannot be summarized. Use `diag(x, nrow = length(x))`. The same pattern appears at `R/plot.effect_curve.R:197` and `R/summary.effect_curve.R:294`, though those see longer vectors in practice.

### ~~24. Simultaneous posterior intervals fail with few draws~~

**Fixed** by the change for issue 6: the single-pass method needs no root finding.

`.pointwise_p_to_simul_p()` calls `uniroot()` on `c(p / (2.5 * k), max(p / 2, level / 2))`. When even the lower end does not reach the target coverage, the bracket has no sign change and `uniroot()` errors. This happened for an *rstanarm* fit with 250 draws and 11 grid points (the default `plot()`), and for *bartisan* with 100 draws. Starting the bracket at 0 (where the band spans all draws and coverage is 1) avoids it.

### ~~25. Families without link functions break `summary()` and `plot()`~~

**Fixed.** `process_family()` rebuilds `linkfun`, `linkinv`, and `mu.eta` with `stats::make.link()` from the link name, or falls back to `quasi("identity")` for an unknown link.

`process_family()` keeps only the `link`, `linkfun`, `linkinv`, and `mu.eta` elements. *bartisan*'s family objects carry only `link`, so `transform` becomes `NULL` and `make_inverse()` fails with "could not find function `fun`". Rebuild missing functions with `stats::make.link(.family$link)` when the link name is known, otherwise fall back to `quasi("identity")`.

### ~~26. Bootstrapping fails when the data have missing values~~

**Fixed** by the change for issue 5: the refit uses only the rows used in the original fit.

The bootstrap weights have one entry per complete case, while `update()` refits on the original data with incomplete rows ("variable lengths differ"). Fixed by the refit approach in issue 5.

### ~~27. Bootstrapping with `cluster` fails~~

**Fixed.** `cl_list[["cluster"]][[1L]]` is passed to *fwb*. No clusters are ignored: with two or more clustering variables (a formula like `~a + b` or a multi-column data frame), `adrf()` stops before bootstrapping with an error that suggests `vcov = "unconditional"` with `cluster`, or a multiway cluster bootstrap of the coefficients supplied as `vcov = fwb::vcovFWB(fit, cluster = ~a + b)` (a conditional variance with Wald inference). Clustering on `~interaction(a, b)` is not suggested because it is one-way clustering on the combinations, not the multiway variance. Checking this exposed two related bugs in `process_cluster()`, both present in the release: a cluster formula with a transformation (e.g., `~interaction(Race, Male)`) failed when the variables were not in the model, and any cluster formula failed for models fit to data with missing values, because `expand.model.frame()` was called with `na.expand = FALSE`, which keeps the dropped rows. It now retrieves the named variables with `na.expand = TRUE` (exactly the rows of `model.frame(model)`) and then evaluates the formula.

`fwb.args[["cluster"]] <- cl_list[["cluster"]]` (`R/adrf.R:314`) passes a list of factors; *fwb* wants a vector. Use `cl_list[["cluster"]][[1L]]`.

### ~~28. Projection models with non-finite values give an unrelated error~~

**Fixed.** The model frames use `na.action = na.pass`, so the finiteness check is reached. `summary.effect_curve()` with a formula hypothesis builds the same model frame and now has the same check.

`curve_projection(a, ~ log(logBLL))` should stop with the "non-finite values" error at `R/curve_projection.R:134`. `model.frame()` drops the `NaN` rows first, so the check passes and `.lm.fit()` fails on mismatched dimensions. Use `model.frame(..., na.action = na.pass)`. `model.matrix.curve_projection()` has the same problem.

### ~~29. *fixest* models fail~~

**Fixed** by the routing change for issue 30: *fixest* models now take the `marginaleffects::get_predict()` path, with unconditional and conditional variances.

`.get_model_type()` lists `"fixest"`, sending `feols()`/`feglm()` down the `lm`/`glm` path, but `get_tt()` and `model.matrix()` do not include the fixed effects, so the coefficients and model matrix do not conform. Either route *fixest* to the `marginaleffects::get_predict()` path or drop it from the list.

### ~~30. `ordinal_weightit` with `type = "link"` fails~~

**Fixed.** `.get_model_type()` no longer relies on the class. It computes the fast-path predictions (model matrix times coefficients, through the inverse link) on up to 20 rows of the counterfactual data and compares them with `marginaleffects::get_predict()`; the fast path (`"lm"` for an identity link, `"glm"` otherwise) is used only if they agree. `lm`, `glm`, `glm.nb`, `glm_weightit`, `lm_weightit`, `svyglm`, and `lm_robust` keep the fast path, and `glm` models with `type = "link"` now take the identity fast path instead of numerical gradients. Models whose predictions involve more than the model matrix and coefficients take the general path: offsets, *fixest* fixed effects, *betareg* precision parameters, and ordinal thresholds. **This also fixed a silent bug present in the release:** `lm()` and `glm()` models with offsets were sent to the fast path, which dropped the offset, giving estimates off by 1.6 and 2.0 in the examples tried. Offset variables supplied through an `offset` argument are now added to the prediction data (`add_offset_vars()`), and `type` is ignored for `lm` objects as documented. `ordinal_weightit` and other multi-category models on the response scale still fail, because their predictions have one value per category; that is the multi-category feature in `_dev/feature-ideas.md`.

`.get_model_type()` special-cases it as `"lm"`, but the model matrix and coefficients (which include thresholds) do not line up. The response scale also fails (see the multi-category item in `_dev/feature-ideas.md`).

### ~~31. `wts`, `subset`, and `by` can refer only to model predictors~~

**Fixed.** Not fixed by issue 1, which only checks the predictors. `process_eval_data()` now recovers the full original dataset for the rows used in the fit (from the survey design for `svyglm` models, then a user-supplied `data`, then the model call's `data` argument evaluated in the formula's environment), checks that it reproduces the model frame, and uses it to evaluate `subset`, `by`, `wts` (as a variable name), and `cluster`; `process_eval_data_mi()` does the same for multiply imputed data. When the dataset cannot be recovered (e.g., a model fit inside a function), `subset` and `by` now give an error that suggests supplying `data`.

The documentation says `subset` "will be evaluated in the environment of the original dataset". `process_model_data()` returns only predictors, so `subset = Race == "White"` and `by = ~Race` fail when `Race` is not in the model, and `wts = "w"` fails unless `w` was used in the fit (`process_wts()` uses `insight::get_data()`, which returns model variables only). Use the full data (`insight::get_data(x, source = "environment")` or `model$data`) for `subset`, `by`, and `wts`.

## Minor issues

Numbered M1, M2, ... for reference. Struck-through items are resolved.

- **~~M1. Typos in messages and docs.~~** ~~`"curve_contast"` (`R/aux_functions.R:86`); `'bs'` for `'boot'` in the MI vcov error (removed with that function when the paths were unified); `object1` for `object` (`R/curve_projection.R:534`); "are produces" and "`ceof()`" (`R/effect_curve.R:6`, `20`); "This computed from ... or from a fitting outcome model directly" (`R/amef.R:4`, and `amef()` does not accept a model); "simple a transformation" (`R/summary.curve_est.R:21`).~~ Resolved: all fixed, along with "an effect curves" and "compliment of the cumulative density function" in `?summary.effect_curve`, "a effect curve" in `?summary.curve_est`, and "ARDF", "for for", and "the the" in the vignette.
- **~~M2. `.gradientFD()` initializes the Jacobian when `j == 1L`~~** ~~(`R/gradient.R:47`). It should test `jj == 1L`; it only works because `.parm` is always `NULL`.~~ Resolved by the maintainer.
- **~~M3. `vcov = "HC3"` (or any non-HC0/1/2 string) with `cluster` silently ignores the string.~~** ~~`process_vcov_and_cluster()` passes `cluster %or% vcov` to `get_vcov()`, so the cluster wins: HC3, HC1, and `"conditional"` with `cluster = ~Race` give identical SEs. Warn or error instead.~~ Resolved: with `cluster`, `vcov` can be `"unconditional"`, `"conditional"` (`vcovCL()` default type), or `"HC0"` to `"HC3"` (passed to `vcovCL(type = )`, using the processed clustering variables, so formulas with transformations and data with missing values work); other strings, functions, and formulas give an error. Checked against `vcovCL()` for all four types, one- and two-way. Tested in `test-adrf.R`.
- **~~M4. Leaked dimnames.~~** ~~`.vcov` carries row and column names like `"1.1"`.~~ Resolved when the paths were unified: the covariance matrix is unnamed before `.make_fun()`.
- **~~M5. Mismatched estimate and SE when an inverse transform fails.~~** ~~In `summary.curve_est()` and `plot()`, `do_transform` is `FALSE` when `inv_transform(est) != est0`, but `est` stays on the transformed scale while `se` is on the original scale. If any estimate is `NaN` after transforming, `do_transform` is `NA` and `if (do_transform)` errors.~~ Resolved: `process_transform()` now checks that the transformation gives finite values for every estimate and that its inverse reverses it (to a tolerance relative to the scale of the estimates), and otherwise uses the identity. This fixes `summary.curve_est()`, `plot()`, `summary.effect_curve()`, and `curve_projection()` at once. Tested in `test-summary.curve_est.R`.
- **~~M6. AMEF at the grid endpoints~~** ~~uses one-sided differences through the local-polynomial interpolation and is about `1e-4` relative off for a quadratic ADRF, against `1e-6` in the interior. Acceptable, but worth a note in `?amef`.~~ Resolved, with a different cause than stated: the finite-difference formula was not the problem (forward, second-order one-sided, and central differences past the end all had errors of about 1e-3 at the ends). `get_locpoly_w()` solved the normal equations of a global orthogonal polynomial basis restricted to a window with kernel weights spanning eight orders of magnitude, losing about seven digits, so the local cubic did not reproduce a quadratic exactly (errors of about 5e-7 on values of about 8), and the derivative amplified this. It now solves the weighted least squares problem by QR, which reproduces polynomials to about 1e-11; endpoint AMEF errors fell to 1e-6 to 1e-4, like the interior. The one-sided stencil was kept. Tested in `test-curve-modifiers.R`.
- **~~M7. `subset` given as a factor~~** ~~fails with "object 's' not found" (`process_subset()` has no branch for it). Character `subset` values that are not row names silently extend the logical vector.~~ Resolved: `process_subset()` gives informative errors for factors, character values that are not row names, non-whole or out-of-range row numbers, and row numbers of mixed sign. Tested in `test-adrf.R`.
- **~~M8. `try(..., silent = FALSE)` in `process_wts()`~~** ~~prints the error before the informative one is thrown. The instance in the string branch was fixed with issue 31; one remains where the number of observations is found (`R/aux_functions.R:539`).~~ Resolved: now `silent = TRUE`.
- **~~M9. Bootstrapping a binomial `glm`~~** ~~emits one "non-integer #successes" warning per replicate.~~ Resolved with issue 5: the refit muffles that warning.
- **~~M10. `posterior_ci()` with `parm` missing~~** ~~fails on `length(parm)`.~~ Resolved with issue 6: it now uses the number of columns of the draws.
- **~~M11. `process_wts_mi()`~~** ~~builds the default `wts` with `unlist()` over a mix of `"(s.weights)"` and logicals.~~ Resolved: the function was removed when the paths were unified, and `process_wts()` now handles each imputed dataset.
- **~~M12. MI degrees of freedom~~** ~~use the minimum complete-data df (`process_df()`); Barnard–Rubin df would be the standard choice.~~ Resolved: `.mi_df()` computes the Reiter (2007) degrees of freedom (Li, Raghunathan, and Rubin, 1991, when the complete-data df are infinite) from the within- and between-imputation covariances, computing the average relative increase in variance in the span of the rank-deficient within-imputation covariance; when Reiter's formula does not apply (k(m − 1) ≤ 4), the large-sample df are combined with the observed-data df as in Barnard and Rubin (1999). A single df is stored, which is exact under the proportional fraction of missing information these methods assume. Matched `mitml:::.D1()` exactly for m = 5 (where, for m = 2 and k = 3, `mitml` returns a negative df). Non-integer df then exposed that `mvtnorm` requires integer df, so df are rounded down for simultaneous inference. Tested in `test-integration-mi.R` and `test-summary.curve_est.R`.
- **~~M13. Unused code.~~** ~~`simultaneous_conf_level()`, `.gradientRich()`, `...mget()`, `all_apply()`, `.ul()`, and `len()` are never called. `plot.effect_curve()` sets `attr(conf_level, "crit")`, which nothing reads. The tests cover these, but removing them is simpler.~~ Resolved by the maintainer; the remaining unused functions are kept because they may be useful later (e.g., for multivariate treatments).
- **~~M14. `print.effect_curve()` hint line~~** ~~is shown only when `interactive()`, so it is untestable without mocking base; no action needed.~~ Resolved: no action needed.
- **~~M15. `wts` documentation.~~** ~~`?adrf` describes `wts` only as a numeric vector, but `process_wts()` also accepts `TRUE` (the model's weights), `FALSE`, and the name of a variable in the original dataset, and the defaults for `svyglm` and *WeightIt* models are described only in Details.~~ Resolved: `@param wts` lists the allowed forms briefly, and a new `wts` section in Details describes each form, the required length of a numeric vector, and the defaults for *survey*, *WeightIt*, and other models.
- **~~M16. *WeightIt* fits with `vcov = "HC0"`.~~** ~~`adrf()` always uses *WeightIt*'s M-estimation estimating functions (`estfun(asympt = TRUE)`), even when the `glm_weightit()` fit was made with `vcov = "HC0"` to treat the weights as fixed. Passing `asympt = (x$vcov_type == "asympt")` would respect that choice (see issue 12).~~ Resolved: for `glm_weightit` objects with `vcov_type == "HC0"`, the influence function uses `estfun(asympt = FALSE)`, so the weights are treated as fixed (the model is then a weighted regression; the unconditional variance matches the same model fit with `lm(weights = )`). Other `vcov_type` values (`"asympt"`, `"BS"`, `"FWB"`, `"none"`) keep the M-estimation estimating functions. This is conditional on the weights, not on the covariates; `adrf()`'s own `vcov` still decides whether the covariate sampling term is added. Tested in `test-integration-weightit.R`.
- **~~M17. Coverage of simultaneous posterior bands.~~** ~~A 95% simultaneous band contained 94.9% of 4,000 draws, because the bands use interpolated quantiles while the level is set from ranks with denominator `R + 1`. *fwb*'s bands behave the same way; acceptable.~~ Resolved by implementing Algorithm 2 of Montiel Olea and Plagborg-Møller (2019) exactly. The band is the product of equal-tailed type 7 percentile intervals with the largest tail probability ζ for which the band contains at least the requested proportion of the draws. A draw of rank r is inside its column's interval exactly when ζ ≤ min(r - 1, R - r)/(R - 1), so ζ is an order statistic of each draw's smallest such bound, which matches bisection over [0, α/2] to rounding. The intermediate code chose ζ from ranks over R + 1 but computed the band with type 7 quantiles, so the band contained fewer draws than requested (92.7% at R = 500, 21 estimates). The paper confines the search to [α/(2k), α/2], but the Bonferroni bound behind that does not hold for discrete quantiles when 2ζR < 1, and the exact solution keeps the "at least 1 - α" definition there. Simultaneous p-values invert the band exactly. The paper's band covers the posterior at the requested level only up to Monte Carlo error, which the paper notes: in simulations, a 95% band covered 94.2% of the posterior with 1000 draws and 21 estimates and 90.2% with 500 draws and 51 estimates. A leave-one-out choice of the rank removed this shortfall (95% to 96% wherever attainable) but was not adopted because it is not described in the literature. Each package now uses one quantile definition consistently in both the choice of the tail probability and the band: *adrftools* uses type 7 quantiles (`fquantile()`, the band reaching the order statistic of rank r at ζ = (r - 1)/(R - 1)), and *fwb* uses its normal interpolation (`norm_inter()`, which reaches the order statistic of rank i exactly at p = i/(R + 1), kept for parity with *boot*). *fwb*'s `sup_t_entry_levels()` and `sup_t_level()` already solve Algorithm 2 exactly under its definition (its band contained 100%, 100%, and 95.0% of the draws for a 95% band with 21 estimates at R = 200, 500, and 4000) and need no change. The intermediate *adrftools* code was wrong only because it copied *fwb*'s R + 1 rank formula but built the band with type 7 quantiles. Tested in `test-internals.R`.
- **~~M18. `fwb::vcovFWB()` cannot evaluate transformed cluster formulas.~~** ~~`vcovFWB(fit, cluster = ~interaction(a, b))` fails with "object 'a' not found", the same problem fixed in *adrftools* with issue 27. This is in *fwb*, not *adrftools*, but it affects the multiway alternative suggested in the bootstrap `cluster` error (a plain formula like `~a + b` works).~~ Resolved: no action needed in *adrftools*; the bootstrap `cluster` error suggests `~a + b`, and `interaction()` should not be recommended for multiway clustering.
- **M19. Small-sample size of the omnibus test.** With the unconditional variance (and conditional HC0), the omnibus tests rejected 7.5% to 9.75% of the time at the .05 level in calibration simulations with n = 500 and a 12-parameter outcome model; conditional HC3 rejected 4.0%. The unconditional variance shares the small-sample bias of HC0 sandwich estimators. A leverage-adjusted (HC2- or HC3-type) unconditional variance would address it (see `_dev/omnibus-pvalue-study/README.md`). **Investigated (2026-11-12, `_dev/omnibus-size-study/README.md`):** part of the excess was an artifact of the calibration study's outcome model, which omitted the `X1:X2` term of the data-generating process, so its ADRF was not exactly flat. With the correct model, the unconditional variance rejects 9.0%, 7.0%, 6.0%, and 5.2% at n = 250, 500, 1000, and 2000. The remaining error is the bias of HC0, concentrated in the high-leverage units at the treatment extremes, not its variability. HC2-adjusted influence values give 5.6%, 5.4%, and 5.0% at n = 250, 500, and 1000; HC3 and the grouped jackknife overcorrect (about 3% at n = 250); the bootstrap and influence-function degrees of freedom do not fix it. For models without hat values, the average of the current variance and a grouped jackknife variance approximates HC2 (4.7% at n = 250). Not yet implemented.
- **~~M20. Liu's method with finite `df`.~~** ~~Liu's approximation is not designed for the negative weight that finite `df` introduces and gave `p = 1.06` for a true value of .001. P-values are now restricted to [0, 1], and `?summary.effect_curve` notes that `"liu"` can be inaccurate when `df` is not `Inf`.~~ Resolved.
- **~~M21. A dataset of the wrong size found by name is used in place of `data`.~~** ~~`process_model_data()` uses the dataset recovered by `insight::get_predictors()` when `.data_matches_model()` accepts it, and `.data_matches_model()` returns `TRUE` when the recovered dataset and the model frame have different numbers of rows (`R/aux_functions.R:369`), presumably for data with missing values. If the model's formula was created in an environment holding a different object with the name of the model's `data` argument (e.g., a formula object defined at the top level and a model fit to a local `d` inside a function, with another `d` at the top level), that object is used even when the correct dataset is supplied to `data`. Found in the M19 study: the unconditional variance failed with "subscript out of bounds", and a bootstrapped curve's estimates did not match the model's. A fix would drop the rows in `na.action` before comparing and treat any remaining row mismatch as a failure to match.~~ Resolved: `.data_matches_model()` now returns `FALSE` when the numbers of rows differ, since every candidate is already restricted to the rows used in the fit (*insight* returns only those rows, handling missing values and `subset`, and `process_eval_data()` aligns candidates by row name first). The fix also covers bootstrap refits, which use the same check. The bug was in 0.1.0, which used the recovered dataset without any check, so it has a NEWS item. Tested in `test-adrf.R`.
