# Size of the omnibus test: study for audit issue M19 (2026-11-12)

The calibration study in `_dev/omnibus-pvalue-study/` found that the omnibus test of a flat ADRF rejected 7.5% to 9.75% of the time at the .05 level with the unconditional variance, and 7.7% at n = 2000. This folder holds the studies used to find out why and to compare fixes. The scripts run from the package root; `common.R` holds the shared code.

## Summary

- **Part of M19 was an artifact of the calibration study.** Its outcome model omitted the `X1:X2` term of the data-generating process. Because the treatment is correlated with `X1` and `X2`, the model attributed part of that term to the treatment, so the ADRF the model estimates was not flat (it departs from flat by up to .039), the null was false, and the rejection rate grew with n. With the term included, the test with the unconditional (HC0-type) variance rejects 9.0%, 7.0%, 6.0%, 5.2%, 5.8%, and 5.2% at n = 250, 500, 1000, 2000, 4000, and 8000 (6000 replicates at n = 2000, 2000 at the others), so its error falls like 1/n and is under a percentage point from n = 2000.
- **What remains is the downward bias of HC0, not its variability.** HC0 rescaled to remove its average bias has the size of a test that uses the true covariance (5.6% and 5.8% at n = 250).
- **HC2 removes the bias; HC3 and the jackknife overcorrect.** The unconditional variance with HC2-adjusted influence values rejects 5.6%, 5.4%, and 5.0% at n = 250, 500, and 1000; HC3 and the grouped jackknife reject about 3% at n = 250.
- **The bootstrap and degrees-of-freedom corrections do not fix it.** The fractional weighted bootstrap covariance rejected 7.1% and 8.3% at n = 250 and 500, and the bootstrap distribution of the statistic 6.3% and 7.1%. Weighting the grid by the density of the treatment also changed little.
- **The error is in the covariance of the coefficients, and a Wald test would not avoid it** (2026-11-16, below). With the exact covariance of the coefficients, the omnibus test is at the nominal level from n = 250; the conditional covariance built from the HC0 sandwich has the same size and bias as the unconditional one. For linear outcome models a Wald test tests the same hypothesis, but it is more sensitive to errors in the covariance (12.6% with HC0 and 8.0% with HC2 at n = 250) and is not practical for nonlinear models.
- **HC2 helps with a logistic model and with heteroskedastic errors but leaves some excess** (2026-11-16, below). It is better than HC0 in every scenario where HC0 over-rejects and nearly unbiased throughout, but with a flexible logistic model (7.9% at n = 500) or with errors largest at the extremes of the treatment (7.8% to 8.0% at n = 250 and 500), the variability of the estimate keeps the test above the nominal level.
- **A general substitute for HC2 is the average of HC0 and the grouped jackknife.** To first order in the leverages, HC2 is the average of HC0 and HC3, and the jackknife approximates HC3 while needing only refits of the model. This average rejected 4.7% at n = 250 (against 8.3% for HC0 in the same replicates) and 6.1% at n = 500 (against 7.8%), and 1.3% at the .01 level at both sizes.

## Why the robust covariance is inaccurate

The test statistic is a weighted sum of squared deviations of the estimated ADRF from its weighted mean, so its null distribution is set by the covariance of the estimates along the directions in which the curve can bend. For a polynomial outcome model, those directions are determined by the units at the extremes of the treatment, whose leverages are much larger than average. Two finite-sample errors of the sandwich estimator follow, both concentrated in the same few units.

**Bias.** The residual of a unit with leverage $h_i$ has variance $\sigma^2 (1 - h_i)$ under homoskedasticity, because the fit is pulled toward that unit (Chesher and Jewitt, 1987, on the bias of HC0; Alene, Van Lancker, and Vansteelandt, 2026, call this "own-observation bias"). HC0 uses squared residuals in place of squared errors, so it underestimates the variance most where the leverages are largest. At n = 250, the average leverage is $p / n = 13 / 250 = .052$, but the HC0 variance of the ADRF at the ends of the grid was 17% to 20% too small (4% in the middle), and the variance of the test statistic, which is dominated by the ends, was 16% too small. HC2 divides each squared residual by $1 - h_i$, which is unbiased under homoskedasticity, and its variance of the test statistic was within 2% of the truth at every sample size.

**Variability.** A sandwich estimate is a sum of one squared influence value per unit. Under homoskedasticity, the Satterthwaite degrees of freedom of the estimated variance of a contrast with unit weights $q_i$ are $(\sum_i q_i^2)^2 / \sum_i q_i^4$ (the Bell and McCaffrey degrees of freedom used by Imbens and Kolesár, 2016), which is $n$ divided by the kurtosis of the $q_i$. For polynomial terms in a normal treatment, that kurtosis is large: 3 for the linear term, 15 for the quadratic, and 92 for the cubic. In the study, these degrees of freedom were about $n / 7$ (35, 63, and 110 at n = 250, 500, and 1000). The sandwich estimator is known to be much more variable than the model-based one (Kauermann and Carroll, 2001), and that variability makes a test that treats it as known anti-conservative, in the same way that using $z$ in place of $t$ is. Here, though, the effect was small: removing the bias alone brought the test to the size of the test with the true covariance.

Both errors shrink like $1/n$, but with large constants, which is why they are visible at n = 500. Neither depends on the outcome model being linear; any model whose ADRF at the extremes of the treatment is driven by a few units has them.

## Results by sample size (`run-diagnostic.R`, `analyze-diagnostic.R`)

The data-generating process is the flat scenario of the calibration study (n units, $A = .6 X_1 + .4 X_2 + \epsilon$, outcome not depending on $A$), with the outcome model `Y ~ poly(A, 3) * (X1 + X2) + X1:X2` (13 coefficients) and the ADRF evaluated at 21 points between the population .025 and .975 quantiles of $A$, so the grid is the same in every sample. There were 2000 replicates per sample size (Monte Carlo SE of a .05 rejection rate about .005), plus 4000 more at n = 2000 (`run-diagnostic-2000.R`), whose results are pooled with the first 2000 below; the oracle at n = 2000 is from the first 2000. "Oracle" uses the covariance of the estimates across the replicates. The unconditional variances with HC1, HC2, and HC3 adjust the coefficient part of each unit's influence value; the covariate term is left unadjusted. "df" adds the influence-function degrees of freedom described below.

Rejection rate at the .05 level:

| n | HC0 (current) | HC1 | HC2 | HC3 | Conditional HC0 | Conditional HC3 | HC0, df | HC2, df | Oracle |
|---|---|---|---|---|---|---|---|---|---|
| 250 | .090 | .082 | .056 | .031 | .092 | .032 | .078 | .048 | .059 |
| 500 | .070 | .065 | .054 | .035 | .071 | .036 | .060 | .045 | .052 |
| 1000 | .060 | .060 | .050 | .040 | .061 | .040 | .056 | .043 | .054 |
| 2000 | .052 | .051 | .047 | .040 | .052 | .040 | .050 | .045 | .052 |
| 4000 | .058 | .057 | .054 | .052 | .058 | .052 | .054 | .052 | .048 |
| 8000 | .052 | .052 | .051 | .050 | .052 | .050 | .051 | .051 | .051 |

Average estimated variance of the test statistic relative to the truth:

| n | HC0 | HC1 | HC2 | HC3 |
|---|---|---|---|---|
| 250 | .841 | .887 | 1.010 | 1.328 |
| 500 | .913 | .938 | 1.019 | 1.183 |
| 1000 | .947 | .959 | 1.012 | 1.103 |
| 2000 | .975 | .981 | 1.013 | 1.062 |

HC1 (the factor $n / (n - p)$) corrects the average leverage, not the largest, so it removes little of the bias. The conditional variance behaves like the unconditional one: the covariate term is not where the problem is.

**Influence-function degrees of freedom.** A Satterthwaite approximation for the variance of the test statistic, $\nu = 2 (\sum_i a_i)^2 / \sum_i a_i^2$ with $a_i$ each unit's contribution to the estimated variance of the statistic, needs only the influence values, so it applies to any model. It overstated the degrees of freedom, by a factor of 2.7 at n = 250 and 1.3 at n = 8000, because the $a_i$ are heavy-tailed, and it corrected only part of the variability. Since the variability matters little here, it is not worth adding.

## Candidate fixes (`run-methods.R`, `analyze-methods.R`)

Each replicate (1500 at n = 250 and 500, 1000 at n = 2000) computed these from the same data:

- **Bootstrap covariance:** the fractional weighted bootstrap (`vcov = "fwb"`, R = 499) in place of the analytic covariance, which is what the package does now with a bootstrapped curve.
- **Bootstrap distribution:** the p-value as the share of bootstrap replicates of the centered statistic, $(\hat\theta^* - \hat\theta)' M (\hat\theta^* - \hat\theta)$, at least as large as the observed one, which does not rely on the chi-squared mixture.
- **Density weights:** the grid points weighted by a kernel density estimate of the treatment in addition to the trapezoidal weights, which down-weights the sparse extremes, as the primitive function in Westling's (2022) test does.
- **Grouped jackknife:** units split at random into K groups, with the model refit and the ADRF recomputed without each group. This deletes each group from both the fit and the average over the covariates, so it estimates the unconditional variance and needs nothing but refits. It has K - 1 degrees of freedom.
- **Averages:** (HC0 + HC3)/2, which equals HC2 to first order in the leverages because $1 / (1 - h)$ is about the average of 1 and $1 / (1 - h)^2$, and (HC0 + grouped jackknife with K = 50)/2, its counterpart for any model that can be refit.
- **HC0 rescaled:** HC0 multiplied by the ratio of the true variance of the statistic to its average HC0 estimate, which removes the average bias and leaves only the variability.

Rejection rate at the .05 level:

| Method | n = 250 | n = 500 | n = 2000 |
|---|---|---|---|
| HC0 (current) | .083 | .078 | .073 |
| HC3 | .028 | .042 | .060 |
| HC0, influence-function df | .070 | .069 | .071 |
| Bootstrap covariance | .071 | .083 | .072 |
| Bootstrap distribution | .063 | .071 | .067 |
| HC0, density weights | .078 | .077 | .078 |
| HC3, density weights | .027 | .049 | .064 |
| Bootstrap distribution, density weights | .065 | .072 | .070 |
| Grouped jackknife, K = 50 | .032 | .042 | .057 |
| Grouped jackknife, K = 50, 49 df | .024 | .038 | .050 |
| Grouped jackknife, K = 20, 19 df | .017 | .033 | .047 |
| (HC0 + HC3)/2 | .046 | .057 | .067 |
| (HC0 + grouped jackknife)/2 | .047 | .061 | .066 |
| HC0 rescaled | .056 | .057 | .054 |
| Oracle | .057 | .056 | .056 |

The 1000 replicates at n = 2000 happened to be about 11% more dispersed than usual (the variance of the statistic across them was 3 Monte Carlo SEs above that of the 6000 replicates in the diagnostic runs), which raised every method's rejection rate by about two points in that column; compare methods within it rather than with the diagnostic table.

At the .01 level, HC0 rejected 2.3% and 2.5% at n = 250 and 500, the bootstrap distribution 1.5% and 1.6%, (HC0 + HC3)/2 1.1% and 1.2%, and (HC0 + grouped jackknife)/2 1.3% at both.

**Reading the table.**

- HC0 rescaled to remove its average bias is as accurate as the oracle, so the bias, not the variability, is what needs fixing. Degrees-of-freedom corrections address the wrong error and, applied to an already unbiased estimator like the jackknife, make the test conservative.
- The bootstrap covariance is about as biased as HC0: its average variance of the statistic was 12% and 10% too small at n = 250 and 500, against 14% and 11% for HC0. Resampling units does not remove own-observation bias, because each resample still fits the model to the units it evaluates. The bootstrap distribution of the statistic helps somewhat, for the same reason a percentile interval differs from a Wald one, but leaves most of the error.
- Density weights help little for HC0 because the bias is in the covariance itself, not only at the extremes. They would also change what the test is sensitive to (departures in the sparse regions of the treatment would count for less), so they are a separate design decision from the size problem.
- HC3 and the grouped jackknife overcorrect by about as much as HC0 undercorrects (their variances of the statistic were 35% and 33% too large at n = 250, and 15% and 14% at n = 500), and the average of HC0 and the jackknife was 9% and 2% too large, landing near the nominal level, as the first-order argument predicts.

## Recommendation

**Implemented 2026-11-16** for models with hat values, without clustering and except `svyglm`: the package's unconditional variance now equals the `u2` variance used here (to machine precision on a replicate at n = 250). Models without hat values use HC0 with a message; the jackknife average below was not implemented.

- **Models with hat values (`lm()`, `glm()`, and fits with a `hatvalues()` method):** compute the unconditional variance with the coefficient part of each unit's influence value divided by $\sqrt{1 - h_i}$ (HC2). This is one line in `.adrf_fit_one()`, costs nothing, and was within a percentage point of the nominal level at every sample size studied. The covariate term needs no adjustment.
- **Other models:** the average of the current influence-function variance and a grouped jackknife variance approximates HC2 with only K refits of the model, using the refitting code that bootstrapping already uses. It is the only general method here that came close to the nominal level, but its justification is the first-order argument above and one simulation, so it would need checking in more settings (a GLM, heteroskedastic errors, other hypotheses) before becoming a default.
- **Not worth pursuing:** influence-function degrees of freedom, density weights (as a size fix), and the bootstrap.

HC2 is exactly unbiased only under homoskedasticity, which this study assumes; with heteroskedastic errors it reduces but does not remove the bias (MacKinnon and White, 1985; Imbens and Kolesár, 2016). Only the flat test and one outcome model were studied.

## The source of the error, and the Wald test (`run-wald.R`, `analyze-wald.R`, 2026-11-16)

This study asks three questions: whether the excess rejection comes from the test itself or from the covariance it is given; whether the covariance of the ADRF estimates is a separate problem from the covariance of the outcome model's coefficients; and whether a Wald test of flatness tests the same hypothesis and gives the same inference. It uses the data-generating process and outcome model above, with 2000 null replicates at each of n = 250, 500, 1000, and 2000 and 1000 replicates under each of four alternatives at n = 500. The covariances compared are:

- `u0`: the current unconditional covariance (HC0-type);
- `u2`: the unconditional covariance with HC2-adjusted influence values;
- `c0`, `c2`: the conditional covariance $G \hat{V}_\beta G'$, with $\hat{V}_\beta$ the HC0 or HC2 sandwich covariance of the coefficients;
- `c_cl`: the conditional covariance with the classical covariance of the coefficients, $\hat{\sigma}^2 (X'X)^{-1}$, which is unbiased and has an exact scaled $\chi^2$ distribution here because the errors are homoskedastic and normal; `c_cl_F` uses it with the residual degrees of freedom;
- `oracle`: the covariance of the estimates across the null replicates.

### The test is accurate when the covariance is

Rejection rate at the .05 level (Monte Carlo SE about .005):

| n | Omnibus `u0` | `u2` | `c0` | `c2` | `c_cl` | `c_cl_F` | Oracle |
|---|---|---|---|---|---|---|---|
| 250 | .087 | .060 | .088 | .062 | .048 | .046 | .059 |
| 500 | .078 | .060 | .078 | .062 | .053 | .052 | .052 |
| 1000 | .064 | .057 | .065 | .057 | .052 | .052 | .051 |
| 2000 | .063 | .059 | .063 | .059 | .054 | .054 | .052 |

With the classical covariance, the omnibus test is at the nominal level at every sample size, including n = 250 (and at the .01 level: .007 to .013). The chi-squared mixture and the normality of the estimates are therefore adequate here; the excess rejection comes from estimating the covariance.

The average estimated variance of the statistic, $\mathrm{tr}(M \hat{V})$, relative to the truth was .873, .893, .935, and .964 for `u0` and .868, .890, .933, and .963 for `c0`, against 1.05, 1.00, 1.00, and 1.00 for `u2` and 1.04, 1.00, .99, and 1.00 for `c_cl`.

The HC2 rates are a point or so higher than in the earlier diagnostic study (.056, .054, .050, and .047 there); the two differ by one to two Monte Carlo SEs, and HC2's bias is gone in both, so what remains is the variability of the estimate.

### The ADRF covariance problem is the coefficient covariance problem

The conditional covariance built from the HC0 sandwich covariance of the coefficients (`c0`) has the same size and the same bias as the unconditional covariance (`u0`), and the same is true of `c2` and `u2`. The covariate-sampling term contributes nothing here, because the treatment has no effect for anyone, so the predicted curves of different units differ only by estimation error. The error is entirely in the sandwich estimate of the coefficients' covariance, as it is passed through $G$ to the ADRF: the own-observation bias of HC0 in the coefficient directions that determine the curve's shape. When effects vary across units, the covariate term is not zero, but it is estimated by a sample variance of predicted values, whose bias is the factor $(N - 1)/N$ and does not depend on leverage. So the fix belongs in the coefficient part of the influence values, as recommended above, not in the ADRF-level covariance.

### A Wald test of the same hypothesis

For this outcome model, the ADRF is a cubic polynomial in the treatment, $\theta(a) = G(a) \beta$, so its deviations from its weighted mean over the grid, $C \theta$, lie in a three-dimensional space (the fourth eigenvalue of their covariance is $10^{-16}$ of the first). A flat ADRF over the grid is then exactly the rank-3 linear hypothesis $C G \beta = 0$. The omnibus statistic and a Wald statistic test this same hypothesis but measure the deviations $d = C \hat{\theta}$ differently: the omnibus statistic uses a fixed metric, $T = d' W d$ (the trapezoid-weighted squared deviations over the treatment range), and the estimated covariance enters only through its reference distribution; the Wald statistic weights the deviations by the inverse of their estimated covariance, $d' (C \hat{V} C')^{+} d$, referred to $\chi^2_3$ (or $F_{3, n - p}$).

Rejection rate at the .05 level of the Wald test, and of the Wald test that all nine coefficients involving the treatment are 0:

| n | Wald `u0` | `u2` | `c_cl` | `c_cl_F` | Oracle | All treatment coefficients, classical F | HC0 | HC2 |
|---|---|---|---|---|---|---|---|---|
| 250 | .126 | .080 | .048 | .046 | .054 | .051 | .384 | .229 |
| 500 | .102 | .080 | .050 | .050 | .052 | .051 | .284 | .186 |
| 1000 | .078 | .064 | .048 | .048 | .050 | .051 | .184 | .136 |
| 2000 | .068 | .059 | .052 | .052 | .060 | .040 | .122 | .096 |

- **With the exact covariance, the two tests have the same size.** Both are at the nominal level with the classical covariance.
- **With an estimated robust covariance, the Wald test is worse.** With HC0, it rejected 12.6% at n = 250 against 8.7% for the omnibus test, and HC2 does not fix it (8.0% at n = 250 and 500). The Wald test inverts the estimated covariance, so it is affected both by HC0's bias, which is largest in the least variable direction, and by its variability, since the inverse of an unbiased covariance estimate is biased upward: the average Wald statistic with HC2 was 3.39 at n = 250, against 3 if the covariance were exact. The omnibus statistic does not involve the estimated covariance, which only sets its reference distribution, so it is less sensitive to both errors.
- **The two tests do not give the same inference in a given sample.** Under the null with the classical covariance, their p-values correlated .91 and they reached different decisions at the .05 level in 3.8% of datasets; of the datasets in which either rejected, both rejected in only 45%. Their power depends on the shape of the departure from flatness (n = 500, classical covariance with F reference):

  | Effect on the outcome | Omnibus | Wald | All treatment coefficients |
  |---|---|---|---|
  | $.1 A$ | .364 | .311 | .267 |
  | $.1 A^2$ | .607 | .644 | .883 |
  | $.3\,\mathrm{logit}^{-1}(4(A - 1))$ | .311 | .300 | .214 |
  | $.03 A^3$ | .471 | .471 | .784 |

  The omnibus test weights a departure by its size over the treatment range, and the Wald test by its size relative to its uncertainty. Here the three eigenvalues of the covariance of the deviations are within a factor of about 5 of each other (relative eigenvalues 1, .68, .22), so the two metrics are similar and so is their power; they would differ more for a more flexible outcome model or a less evenly spread treatment.
- **The test of all treatment coefficients tests a different hypothesis.** Its null is that the treatment has no effect for any covariate pattern, which is stronger than a flat ADRF: it rejects when effects vary across units but average to a flat curve. When the effect is the same for everyone, as in these alternatives, it can be more powerful than either ADRF test, because it uses the restriction that the interaction coefficients are 0. With robust covariances it is far worse in size (38% with HC0 at n = 250), because it inverts a 9-dimensional covariance that includes the poorly estimated interaction directions.
- **For a nonlinear outcome model, a Wald test is not a practical alternative.** With a logistic model, the ADRF is not a polynomial and flatness is not a linear restriction on the coefficients. The covariance of the deviations is numerically full rank, with eigenvalues trailing off (relative eigenvalues 1, .61, .20, .0025, .0004, .0001, ... in one sample of 1000), so a Wald test would have to choose a rank and invert eigenvalues estimated with little precision. The omnibus test needs no inversion.

### Conclusions

- The size problem is the bias of the sandwich estimate of the coefficients' covariance, propagated to the ADRF; the omnibus test is accurate with an accurate covariance.
- A Wald test would not fix it: it tests the same hypothesis for linear outcome models but is more sensitive to errors in the covariance, so it is worse with robust covariances, and it is not well defined for nonlinear models.
- HC2 removes the bias but leaves the test a point or so above the nominal level at these sample sizes in this run, from the variability of the estimate.

## HC2 with a logistic model and with heteroskedastic errors (`run-glm-het.R`, `analyze-glm-het.R`, 2026-11-16)

After the HC2 adjustment became the package's default, this study checked it outside the setting it was chosen in. Each cell has 2000 replicates (Monte Carlo SE of a .05 rejection rate about .005); the ADRF is flat in all of them, and the outcome model is the 13-coefficient model above.

- **Logistic model:** a binary outcome with $P(Y = 1) = \mathrm{logit}^{-1}(-.5 + .8 X_1 - .5 X_2 + .3 X_1 X_2)$, not depending on the treatment, fit with `glm(family = binomial)` at n = 500, 1000, and 2000. HC0 and HC2 come from `adrf()` itself, with the adjustment switched off for HC0. The flat test gave the same p-value with and without the default logit transformation (to $10^{-13}$), because the transformed flat projection is the constant on the response scale that minimizes the weighted squared deviations, which is the weighted mean.
- **Heteroskedastic errors**, with the linear model at n = 250, 500, and 1000, the errors normal with an average variance near 1 and SD:
  - *increasing*: $.43 + .52 |A - .2|$, three times as large at the extremes of the treatment as in its middle, where the leverages are largest;
  - *decreasing*: $1.6 / (1 + |A - .2|)$, a third as large at the extremes;
  - *covariate*: $\exp(.5 X_1 - .25)$, depending on a covariate but not on the treatment.

Rejection rate at the .05 level, and the average estimated variance of the statistic relative to the truth:

| Scenario | n | HC0 | HC2 | HC3 | Classical | Oracle | Variance ratio, HC0 | HC2 |
|---|---|---|---|---|---|---|---|---|
| Logistic | 500 | .098 | .079 | | | .051 | .86 | .97 |
| | 1000 | .076 | .066 | | | .045 | .89 | .96 |
| | 2000 | .065 | .060 | | | .055 | .95 | .99 |
| Increasing | 250 | .128 | .078 | .031 | .260 | .056 | .75 | .97 |
| | 500 | .104 | .080 | .043 | .284 | .053 | .80 | .94 |
| | 1000 | .069 | .050 | .034 | .292 | .050 | .91 | 1.01 |
| Decreasing | 250 | .052 | .036 | .023 | .002 | .057 | .97 | 1.10 |
| | 500 | .052 | .045 | .034 | .002 | .056 | .98 | 1.04 |
| | 1000 | .050 | .044 | .038 | .000 | .052 | 1.02 | 1.05 |
| Covariate | 250 | .093 | .057 | .028 | .068 | .061 | .80 | .99 |
| | 500 | .072 | .045 | .029 | .060 | .050 | .84 | .98 |
| | 1000 | .070 | .052 | .040 | .066 | .056 | .90 | 1.00 |

"Classical" is the conditional covariance with $\hat\sigma^2 (X'X)^{-1}$ and "Oracle" uses the covariance of the estimates across the replicates; HC3 and the classical covariance were computed only for the linear models.

- **The oracle test is at the nominal level in every cell**, including the logistic model, so the chi-squared mixture remains adequate and the errors are again in the covariance.
- **HC2 is better than HC0 in every scenario except one, where HC0 was already accurate.** It nearly removes the bias of the variance of the statistic everywhere (ratios .94 to 1.05 against .75 to 1.02 for HC0). When the errors are smaller at the extremes, HC0's bias is offset (the high-leverage units contribute little variance), HC0 is at the nominal level, and HC2 is slightly conservative (3.6% at n = 250, 4.4% to 4.5% from n = 500).
- **HC2 does not fully fix the size when the high-leverage units are also noisy.** With errors larger at the extremes, HC2 rejects 7.8% and 8.0% at n = 250 and 500 (12.8% and 10.4% for HC0), and with the logistic model, whose largest leverage was .45 in a typical sample of 500, 7.9% and 6.6% at n = 500 and 1000. HC2's variance of the statistic is close to unbiased in both, so what remains is its variability: the estimate rests on a few units with high leverage and large residual variance, so it has few effective degrees of freedom, the second error described above, which mattered little in the homoskedastic study.
- **HC3 overcorrects in every case** (2.3% to 4.3%), as before.
- **The classical covariance is badly wrong under heteroskedasticity** (26% to 29% when errors grow at the extremes, 0.2% when they shrink), as expected; the conditional variance in `adrf()` uses the model's own covariance by default, which is the classical one for `lm()`.

**Conclusion.** The HC2 default is a clear improvement over HC0 across these scenarios and never much worse, but with a flexible logistic model at n of 1000 or less, or with error variance concentrated at the extremes of the treatment at n of 500 or less, the omnibus test still rejects 6.6% to 8.0% at the .05 level. The remaining error is variability rather than bias, so a fix would be a degrees-of-freedom correction suited to heavy-tailed contributions or a different reference distribution, not a different leverage adjustment; the influence-function degrees of freedom tried earlier overstated the degrees of freedom by a factor of 2.7 at n = 250 and would not suffice.

## The calibration study's misspecified model (`misspecified/`)

The partial results of the first run of `run-diagnostic.R`, with the calibration study's model (no `X1:X2`), are kept in `misspecified/`. With 2000 replicates at n = 250 to 2000, its oracle test rejected 6.0%, 5.4%, 5.9%, and 6.3%; the noncentrality of the statistic relative to its null mean was .024 at n = 500 and .11 at n = 2000 and grows in proportion to n. The README of `_dev/omnibus-pvalue-study/` has been corrected.

## References

Alene, M., Van Lancker, K., & Vansteelandt, S. (2026). On the use of G-computation in small randomized controlled trials with many covariates. arXiv:2609.37630. https://doi.org/10.48550/arXiv.2609.37630

Chesher, A., & Jewitt, I. (1987). The bias of a heteroskedasticity consistent covariance matrix estimator. *Econometrica*, 55(5), 1217–1222. https://doi.org/10.2307/1911269

Imbens, G. W., & Kolesár, M. (2016). Robust standard errors in small samples: Some practical advice. *Review of Economics and Statistics*, 98(4), 701–712. https://doi.org/10.1162/REST_a_00552

Kauermann, G., & Carroll, R. J. (2001). A note on the efficiency of sandwich covariance matrix estimation. *Journal of the American Statistical Association*, 96(456), 1387–1396. https://doi.org/10.1198/016214501753382309

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305–325. https://doi.org/10.1016/0304-4076(85)90158-7

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *Journal of Business & Economic Statistics*, 36(4), 672–683. https://doi.org/10.1080/07350015.2016.1247004

Westling, T. (2022). Nonparametric tests of the causal null with nondiscrete exposures. *Journal of the American Statistical Association*, 117(539), 1551–1562. https://doi.org/10.1080/01621459.2020.1865168
