# Weights in the unconditional variance of the ADRF

Written 2026-11-12 alongside the fix for audit issue 4. It derives the unconditional variance `adrf()` computes, shows where each kind of weight enters, and covers three general cases: target and model weights (issue 4), complex survey designs (issue 13), and estimated balancing weights (issue 12). Each result is checked numerically against independent software in `_dev/variance-verification/verify-variance.R`; the results are in `_dev/variance-verification/results.csv` and summarized in the section on verification.

## Summary

- **Issue 4 is fixed.** The covariate-sampling term of the influence function is $N w_i (m_i - \mu) / W_b$. The previous code used $(m_i - \mu) W / W_b$, which agrees only when every weight is 1. After the fix, `adrf()` matches *stdReg* (Sjölander, 2016) and *geex* (Saul & Hudgens, 2020) to every printed digit, with and without subgroups, clustering, and model weights that differ from the target weights.
- **Issue 13 is fixed.** For `svyglm` models, `adrf()` computed an i.i.d. variance that ignored PSUs, strata, and finite-population corrections; in a stratified cluster design it differed from a design-based jackknife and a hand-computed Binder (1983) linearization by up to 16%. It now passes its influence values to the design's own variance estimator and matches the linearization exactly.
- **Issue 12 was a mistake in the audit and is retracted.** `glm_weightit` objects have `estfun()` and `bread()` methods (from *WeightIt*) whose estimating functions already include the contribution of the weight model, so `adrf()` inherits the M-estimation correction. The audit's evidence (identical results for fits with `vcov = "asympt"` and `vcov = "HC0"`) reflected that `estfun()` reads the attached `weightit` object regardless of the fit's `vcov` setting, not that the correction was missing.

## The estimator

Let units $i = 1, \dots, N$ make up the estimation sample. An outcome model with parameters $\beta$ is fit by solving

$$\sum_{i=1}^N \psi_i(\beta) = 0,$$

where $\psi_i$ is the estimating function of unit $i$ (for weighted least squares, $\psi_i(\beta) = v_i x_i (y_i - x_i^\top \beta)$, with model weights $v_i$). Write $m_i(a; \beta) = E[Y \mid A = a, X = x_i; \beta]$ for the prediction for unit $i$ with its treatment set to $a$.

The ADRF at $a$ in subgroup $b$ is estimated as a weighted mean of predictions:

$$\hat\mu_b(a) = \frac{\sum_i b_i w_i \, m_i(a; \hat\beta)}{\sum_i b_i w_i},$$

where $b_i$ indicates membership in the subgroup (all 1 without `by` or `subset`), and $w_i$ are the **target weights** supplied to `wts`, which describe the population the curve is averaged over. Write $W_b = \sum_i b_i w_i$.

Two kinds of weights appear, and they play different roles:

- **Model weights** $v_i$ enter $\psi_i$ and affect $\hat\beta$ (balancing weights, survey weights, or any `weights` argument to the fitting function).
- **Target weights** $w_i$ enter only the averaging step and define the population.

They coincide for survey analyses, where the survey weights are both, and differ otherwise. For example, a model fit with inverse probability weights uses $v_i$ to balance confounders while $w_i = 1$ targets the sample. By default, `adrf()` sets $w_i$ to the survey weights for `svyglm` models, to the sampling weights for *WeightIt* models, and to 1 otherwise.

## The unconditional variance by stacked M-estimation

Treat $\theta = (\beta, \mu_b(a))$ as the solution to the stacked estimating equations

$$\sum_i \begin{pmatrix} \psi_i(\beta) \\ b_i w_i \{ m_i(a; \beta) - \mu \} \end{pmatrix} = 0.$$

Stacking estimating functions this way is the standard device of M-estimation (Stefanski & Boos, 2002). Because the second equation involves $\beta$ but the first does not involve $\mu$, the Jacobian is block lower-triangular. With

$$A_{\beta\beta} = -\frac{1}{N}\sum_i \frac{\partial \psi_i}{\partial \beta^\top}, \qquad \bar g_b(a) = \frac{1}{W_b} \sum_i b_i w_i \frac{\partial m_i(a;\beta)}{\partial \beta},$$

inverting the block triangle gives the influence function of $\hat\mu_b(a)$ for unit $i$:

$$\mathrm{IF}_i\{\mu_b(a)\} = \underbrace{\frac{N}{W_b} \, b_i w_i \, \{ m_i(a; \hat\beta) - \hat\mu_b(a) \}}_{\text{sampling of the covariates}} + \underbrace{\bar g_b(a)^\top A_{\beta\beta}^{-1} \psi_i(\hat\beta)}_{\text{estimation of } \beta}.$$

The variance is estimated by

$$\widehat{\mathrm{Var}}\{\hat\mu_b(a)\} = \frac{1}{N^2} \sum_i \mathrm{IF}_i^2,$$

and covariances between treatment values and subgroups by the corresponding cross products. This is what `adrf()` computes with `vcov = "unconditional"`: each row of its matrix `S` holds $\mathrm{IF}_i$ for every (subgroup, value) pair, the second term is `sandwich::estfun(x) %*% sandwich::bread(x)` times the gradient of the averaged predictions, and the variance is `crossprod(S / N)`.

When $w_i = 1$ and $b_i = 1$, this is the variance of a g-computation potential outcome mean derived by Hansen and Overgaard (2024), who also cover subgroup means such as the effect among the treated. Graubard and Korn (1999) gave the analogous result for predictive margins with survey weights and noted that the "textbook" variance, which keeps only the second term, can be inappropriate when the covariate distribution is itself estimated. The second term alone is `adrf()`'s `vcov = "conditional"`, which treats the covariates as fixed.

### The error fixed in issue 4

The previous code used $(m_i - \hat\mu_b) \, W / W_b$ for the first term, with $W = \sum_i w_i$. It lacked the factor $w_i$ and used $W$ where $N$ belongs. The two expressions agree only when every $w_i = 1$, so the error affected every analysis with non-uniform target weights: user-supplied `wts`, `svyglm` models, and *WeightIt* models with sampling weights. The size of the error depends on how variable the weights are. In the verification (check 1b) it changed the SEs by 0.7% to 2.9%.

The influence function is invariant to rescaling the weights, since $w_i / W_b$ is.

### Clustering

With clusters $g = 1, \dots, G$, the sum of squares is replaced by sums within clusters, with the small-sample factor `sandwich::vcovCL()` uses:

$$\widehat{\mathrm{Var}} = \frac{G}{G-1} \frac{1}{N^2} \sum_g \Big( \sum_{i \in g} \mathrm{IF}_i \Big)^2.$$

## Survey designs (issue 13)

Under design-based inference, $\hat\beta$ solves design-weighted estimating equations and $\hat\mu_b(a)$ is a ratio of design-weighted totals, with $v_i = w_i$ equal to the survey weights. Binder (1983) showed that such estimators are linearized as

$$\hat\theta - \theta \approx \sum_i w_i z_i,$$

with $z_i$ the unit's influence value scaled by the estimated population size, so that the variance of $\hat\theta$ is the **design variance of an estimated total** of $z_i$. That variance is computed from PSU totals within strata, with finite-population corrections, which is what `survey::svyrecvar()` does for every estimator in *survey* (Lumley, 2004). For the ADRF, $N w_i z_i$ is exactly the $\mathrm{IF}_i$ above, with the survey weights in both roles.

Before the fix, `adrf()` summed $\mathrm{IF}_i^2$ over units, which is the design variance only for a with-replacement design with one unit per PSU and no strata (up to the factor $n/(n-1)$). That is why it agreed with the jackknife for a weights-only design and not for a clustered, stratified one.

**The fix for issue 13.** `.svyglm_vcov()` divides each row of `S / N` by the unit's sampling weight and totals the result as design variables with `survey::svytotal()`. Because $\sum_i w_i (\mathrm{IF}_i / N w_i) = \sum_i \mathrm{IF}_i / N$, this is the linearized total, and *survey* chooses the variance estimator for the design: `svyrecvar()` with clusters, strata, finite population corrections, and calibration for linearization designs, the PPS estimators for PPS designs, and the replicate weights for replicate-weight designs. For replicate designs, totaling fixed influence values with each set of replicate weights gives the same answer as linearization (up to how *survey* centers replicates), which the tests check. `cluster` is ignored with a warning for `svyglm` models, since the design carries its own clusters.

The fix also exposed a related bug: `insight::get_weights()` returns `NULL` for `svyglm` models fit to replicate-weight designs, so `adrf()` averaged their predictions without the survey weights. The target weights for `svyglm` models now come from `prior.weights`.

## Estimated balancing weights (issue 12)

Suppose the model weights are estimated, $v_i = v_i(\hat\gamma)$, from a treatment model solving $\sum_i \phi_i(\gamma) = 0$. For a continuous treatment, the stabilized weight is a ratio of normal densities, $v_i = f(A_i; \alpha, \tau) / f(A_i \mid Z_i; \gamma, \sigma)$ (Robins, Hernán, & Brumback, 2000; Naimi et al., 2014). Adding $\phi_i$ to the stack makes the Jacobian block lower-triangular in $(\gamma, \beta, \mu)$, and the coefficient term of the influence function becomes

$$A_{\beta\beta}^{-1} \tilde\psi_i, \qquad \tilde\psi_i = \psi_i + \Big( \frac{1}{N} \sum_j \frac{\partial \psi_j}{\partial \gamma^\top} \Big) A_{\gamma\gamma}^{-1} \phi_i, \qquad A_{\gamma\gamma} = -\frac{1}{N}\sum_j \frac{\partial \phi_j}{\partial \gamma^\top},$$

that is, the outcome model's estimating functions augmented by a projection of the weight model's. The covariate-sampling term is unchanged. Treating the weights as known drops the correction. Lunceford and Davidian (2004) used this joint estimating-equation approach for inverse probability weighted estimators of the average treatment effect. Reifeis and Hudgens (2022) showed that ignoring weight estimation gives conservative SEs for the ATE but can be anticonservative for the ATT, so the correction is not optional in general. Gabriel et al. (2024) studied the estimator `adrf()` computes here (a GLM fit with inverse probability weights, followed by standardization) and showed it is doubly robust when the GLM uses the canonical link.

**Why issue 12 is not a bug.** *WeightIt*'s `estfun.glm_weightit()` returns $\tilde\psi_i$ when the fit has a `weightit` object with M-estimation components (`asympt = TRUE` is the default): it adds `psi_treat %*% solve(-H_treat, t(H_out_treat))` to the outcome model's estimating functions. `adrf()` calls `sandwich::estfun()`, so it uses $\tilde\psi_i$. Check 5b confirms this: `adrf(glm_weightit)` equals the delta method applied to *WeightIt*'s own M-estimation `vcov()`, and both differ from the same model with the weights treated as fixed.

Two caveats remain:

- **The fit's `vcov` setting is ignored.** A `glm_weightit()` fit with `vcov = "HC0"` (weights treated as fixed) still gets the corrected estimating functions in `adrf()`. Passing `asympt = (x$vcov_type == "asympt")` to `estfun()` would respect the user's choice. This is a design decision, not a bug.
- ***WeightIt*'s numerator density is shared across units.** For continuous treatments, *WeightIt*'s numerator is the model-implied marginal density, an average of the fitted conditional densities over all units. A unit's weight therefore depends on every unit, which is not a per-unit estimating function. *WeightIt* handles this by differentiating the whole weight vector numerically. The bootstrap in check 5b, which re-estimates the weights in each replicate, supports the result.

## Bootstrapping with weights (issue 5)

For the bootstrap, each replicate multiplies the unit weights by bootstrap weights $u_i$, so the model must be refit with $v_i u_i$ (not $u_i$ alone, and not $w_i u_i$) and the predictions averaged with $w_i u_i$. When the model weights are estimated, they must be re-estimated in each replicate, with the weighting model fit using $u_i$ times its sampling weights. Before the fix for issue 5, `adrf()` refit with $w_i u_i$, which is correct only when $v_i = w_i$. It now refits with $v_i u_i$, and for models fit with a `weightit` object it re-estimates the weights as *WeightIt*'s own bootstrap does; with the same seed, the two give identical results.

Resampling units is not a valid bootstrap for a complex survey design, which is why bootstrapping `svyglm` models is now an error. Replicate-weight designs are the design-based alternative (Rust & Rao, 1996).

## Verification

All checks use `nhanes3lead` (N = 2,521), the ADRF at logBLL = 0, 1, and 2, and are run by `_dev/variance-verification/verify-variance.R`. *geex* and *stdReg* are not package dependencies; the script requires them to be installed.

### The independent implementations

| Implementation | Source | What it computes | Described in |
|---|---|---|---|
| `stdReg::stdGlm()` | CRAN, *stdReg* 3.4.2 | Regression standardization of a GLM with a stacked sandwich variance; with a weighted GLM, the prior weights are used both to fit and to average, and its residual term is `weights * (pred - est)` | Sjölander (2016) |
| `geex::m_estimate()` | CRAN, *geex* 1.1.1 | Generic stacked M-estimation: I wrote the estimating functions (outcome GLM score with $v_i$, the weighted mean with $w_i$ and $b_i$, and for check 5 the treatment model and weights), and *geex* computes the Jacobian numerically and the sandwich | Saul and Hudgens (2020); Stefanski and Boos (2002) |
| `survey::withReplicates()` | *survey* 4.5 | Design-based jackknife (JKn, JK1): the outcome model is refit and the ADRF recomputed with each set of replicate weights | Rust and Rao (1996); Lumley (2004) |
| `survey::svytotal()` on influence values | *survey* 4.5 | Binder linearization: the influence values computed by hand in the script are totaled with the design's variance estimator | Binder (1983); Lumley (2004) |
| `fwb::fwb()` | *fwb* 1.0.0 | Fractional weighted bootstrap of the whole procedure, 2,000 replicates (Monte Carlo error of an SE about 1.6%) | Greifer (2026) |

*stdReg* computes its middle matrix with `var()`, whose denominator is $n - 1$; `adrf()` and *geex* use $1/n$. Its unclustered results are therefore multiplied by $(n - 1)/n$ in the comparison. With clusters, its `var()` over cluster totals times $G/n^2$ equals the $G/(G-1)$ adjustment above, so no rescaling is needed.

### Results

Standard errors at logBLL = 0, 1, 2.

| Check | Method | 0 | 1 | 2 |
|---|---|---|---|---|
| 1a unweighted logistic | *stdReg* × (n−1)/n | 0.014552 | 0.0080306 | 0.012070 |
| | `adrf()` | 0.014552 | 0.0080306 | 0.012070 |
| 1b weighted linear, $v = w$ | *stdReg* × (n−1)/n | 0.30043 | 0.15537 | 0.26992 |
| | `adrf()` | 0.30043 | 0.15537 | 0.26992 |
| | `adrf()` before the fix | 0.29405 | 0.15095 | 0.27192 |
| 1c weighted logistic, 150 clusters | *stdReg* | 0.037706 | 0.015680 | 0.020537 |
| | `adrf()` | 0.037706 | 0.015680 | 0.020537 |
| 2a logistic, $v \ne w$ | *geex* | 0.029751 | 0.016022 | 0.023374 |
| | `adrf()` | 0.029751 | 0.016022 | 0.023374 |
| | *fwb* | 0.029245 | 0.016237 | 0.023555 |
| 2b as 2a, subgroup Male = 1 | *geex* | 0.040663 | 0.023033 | 0.034580 |
| | `adrf()` | 0.040663 | 0.023033 | 0.034580 |
| 4 stratified cluster design (12 strata, 96 PSUs) | `adrf(svyglm)` | 0.17593 | 0.10276 | 0.18846 |
| | `adrf(svyglm)` before the fix | 0.17059 | 0.10955 | 0.21971 |
| | *survey* JKn | 0.17731 | 0.10306 | 0.18932 |
| | linearization by hand, `svytotal()` | 0.17593 | 0.10276 | 0.18846 |
| 4b weights-only design | `adrf(svyglm)` | 0.17063 | 0.10957 | 0.21975 |
| | `adrf(svyglm)` before the fix | 0.17059 | 0.10955 | 0.21971 |
| | *survey* JK1 | 0.17210 | 0.11010 | 0.22188 |
| 5 stabilized weights, adjusted outcome model | *geex*, weights estimated | 0.13144 | 0.073735 | 0.13916 |
| | *geex*, weights fixed | 0.13202 | 0.073795 | 0.13998 |
| | `adrf(lm(weights = sw))` | 0.13202 | 0.073795 | 0.13998 |
| | *fwb*, weights re-estimated | 0.13318 | 0.074366 | 0.13846 |
| 5c stabilized weights, MSM outcome model | *geex*, weights estimated | 0.15076 | 0.086176 | 0.15465 |
| | *geex*, weights fixed | 0.15409 | 0.086776 | 0.16040 |
| | `adrf(lm(weights = sw))` | 0.15409 | 0.086776 | 0.16040 |
| | *fwb*, weights re-estimated | 0.15474 | 0.086158 | 0.15493 |
| 5b *WeightIt*, MSM outcome model | `adrf(glm_weightit)` | 0.15110 | 0.086206 | 0.15740 |
| | delta method on `vcov(glm_weightit)` | 0.15110 | 0.086206 | 0.15740 |
| | `adrf(lm(weights = w))`, weights fixed | 0.15437 | 0.086889 | 0.16313 |
| | *fwb*, weights re-estimated | 0.15373 | 0.085116 | 0.15647 |

How to read these:

- **Checks 1 and 2 (issue 4).** `adrf()` matches *stdReg* and *geex* exactly in every configuration, and the bootstrap agrees within its Monte Carlo error. The pre-fix formula does not match.
- **Check 4 (issue 13).** The two design-based methods agree with each other (they estimate the same quantity by different routes). Before the fix, `adrf(svyglm)` differed from them by −4% to +16%, in both directions; it now equals the hand-computed linearization and is within 1% of the jackknife. The tests also compare replicate-weight, bootstrap-replicate, calibrated, and domain designs with their jackknife or replicate variances, all within about 1%. For a weights-only design, the fix changes only the factor $n/(n-1)$.
- **Check 5 (estimated weights).** With an outcome model that adjusts for the confounders, accounting for weight estimation changes the SEs by under 1%, less than the bootstrap's Monte Carlo error. With a marginal structural model in the treatment alone, the weights carry all of the adjustment, and treating them as fixed overstates the SE at logBLL = 2 by 3.7% (0.1604 vs. 0.1547). The bootstrap that re-estimates the weights gives 0.1549, matching the full stack. `adrf()` on a weighted `lm()` corresponds to the "weights fixed" row, as it should, since the `lm()` object knows nothing of how its weights were made.
- **Check 5b (issue 12).** `adrf(glm_weightit)` equals *WeightIt*'s M-estimation result and differs from the fixed-weight version, confirming that the correction is already applied.

## References

Binder, D. A. (1983). On the variances of asymptotically normal estimators from complex surveys. *International Statistical Review*, 51(3), 279–292. https://doi.org/10.2307/1402588

Gabriel, E. E., Sachs, M. C., Martinussen, T., Waernbaum, I., Goetghebeur, E., Vansteelandt, S., & Sjölander, A. (2024). Inverse probability of treatment weighting with generalized linear outcome models for doubly robust estimation. *Statistics in Medicine*, 43(3), 534–547. https://doi.org/10.1002/sim.9969

Graubard, B. I., & Korn, E. L. (1999). Predictive margins with survey data. *Biometrics*, 55(2), 652–659. https://doi.org/10.1111/j.0006-341X.1999.00652.x

Greifer, N. (2026). fwb: An R package for the fractional weighted bootstrap. *Journal of Open Source Software*, 11(124), 8880. https://doi.org/10.21105/joss.08880

Hansen, S. N., & Overgaard, M. (2024). Variance estimation for average treatment effects estimated by g-computation. *Metrika*. https://doi.org/10.1007/s00184-024-00962-4

Lumley, T. (2004). Analysis of complex survey samples. *Journal of Statistical Software*, 9(8). https://doi.org/10.18637/jss.v009.i08

Lunceford, J. K., & Davidian, M. (2004). Stratification and weighting via the propensity score in estimation of causal treatment effects: A comparative study. *Statistics in Medicine*, 23(19), 2937–2960. https://doi.org/10.1002/sim.1903

Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014). Constructing inverse probability weights for continuous exposures: A comparison of methods. *Epidemiology*, 25(2), 292–299. https://doi.org/10.1097/EDE.0000000000000053

Reifeis, S. A., & Hudgens, M. G. (2022). On variance of the treatment effect in the treated when estimated by inverse probability weighting. *American Journal of Epidemiology*, 191(6), 1092–1097. https://doi.org/10.1093/aje/kwac014

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal structural models and causal inference in epidemiology. *Epidemiology*, 11(5), 550–560. https://doi.org/10.1097/00001648-200009000-00011

Rust, K. F., & Rao, J. N. K. (1996). Variance estimation for complex surveys using replication techniques. *Statistical Methods in Medical Research*, 5(3), 283–310. https://doi.org/10.1177/096228029600500305

Saul, B. C., & Hudgens, M. G. (2020). The calculus of M-estimation in R with geex. *Journal of Statistical Software*, 92(2). https://doi.org/10.18637/jss.v092.i02

Sjölander, A. (2016). Regression standardization with the R package stdReg. *European Journal of Epidemiology*, 31(6), 563–574. https://doi.org/10.1007/s10654-016-0157-3

Stefanski, L. A., & Boos, D. D. (2002). The calculus of M-estimation. *The American Statistician*, 56(1), 29–38. https://doi.org/10.1198/000313002753631330

### Citation sources

Ten references are in the Zotero library and were taken from it. Binder (1983), Rust and Rao (1996), Lumley (2004), and Saul and Hudgens (2020) are not in the library; their DOIs were confirmed to resolve and their titles, authors, and volumes were taken from Crossref. The claims attributed to Hansen and Overgaard, Graubard and Korn, Gabriel et al., and Reifeis and Hudgens are from their abstracts; the claim attributed to Lunceford and Davidian is from their Section 3.1. The other attributions are at the level of each paper's subject and were not checked against the text.
