# Size of the omnibus test: study for audit issue M19 (2026-11-12)

The calibration study in `_dev/omnibus-pvalue-study/` found that the omnibus test of a flat ADRF rejected 7.5% to 9.75% of the time at the .05 level with the unconditional variance, and 7.7% at n = 2000. This folder holds the studies used to find out why and to compare fixes. The scripts run from the package root; `common.R` holds the shared code.

## Summary

- **Part of M19 was an artifact of the calibration study.** Its outcome model omitted the `X1:X2` term of the data-generating process. Because the treatment is correlated with `X1` and `X2`, the model attributed part of that term to the treatment, so the ADRF the model estimates was not flat (it departs from flat by up to .039), the null was false, and the rejection rate grew with n. With the term included, the test with the unconditional (HC0-type) variance rejects 9.0%, 7.0%, 6.0%, 5.2%, 5.8%, and 5.2% at n = 250, 500, 1000, 2000, 4000, and 8000 (6000 replicates at n = 2000, 2000 at the others), so its error falls like 1/n and is under a percentage point from n = 2000.
- **What remains is the downward bias of HC0, not its variability.** HC0 rescaled to remove its average bias has the size of a test that uses the true covariance (5.6% and 5.8% at n = 250).
- **HC2 removes the bias; HC3 and the jackknife overcorrect.** The unconditional variance with HC2-adjusted influence values rejects 5.6%, 5.4%, and 5.0% at n = 250, 500, and 1000; HC3 and the grouped jackknife reject about 3% at n = 250.
- **The bootstrap and degrees-of-freedom corrections do not fix it.** The fractional weighted bootstrap covariance rejected 7.1% and 8.3% at n = 250 and 500, and the bootstrap distribution of the statistic 6.3% and 7.1%. Weighting the grid by the density of the treatment also changed little.
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

- **Models with hat values (`lm()`, `glm()`, and fits with a `hatvalues()` method):** compute the unconditional variance with the coefficient part of each unit's influence value divided by $\sqrt{1 - h_i}$ (HC2). This is one line in `.adrf_fit_one()`, costs nothing, and was within a percentage point of the nominal level at every sample size studied. The covariate term needs no adjustment.
- **Other models:** the average of the current influence-function variance and a grouped jackknife variance approximates HC2 with only K refits of the model, using the refitting code that bootstrapping already uses. It is the only general method here that came close to the nominal level, but its justification is the first-order argument above and one simulation, so it would need checking in more settings (a GLM, heteroskedastic errors, other hypotheses) before becoming a default.
- **Not worth pursuing:** influence-function degrees of freedom, density weights (as a size fix), and the bootstrap.

HC2 is exactly unbiased only under homoskedasticity, which this study assumes; with heteroskedastic errors it reduces but does not remove the bias (MacKinnon and White, 1985; Imbens and Kolesár, 2016). Only the flat test and one outcome model were studied.

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
