# Choice of omnibus test statistic (2026-11-16)

Tests of equivalent null hypotheses (a flat ADRF, a reference curve equal to 0, an AMEF equal to 0; a linear ADRF, a flat AMEF) use different statistics and can reach different conclusions on the same data (`_dev/omnibus-size-study/README.md`, equivalent hypotheses). This study compares the statistics directly, to find one to recommend for each hypothesis. Only the statistic varies: every test in a replicate uses the same covariance matrix of the ADRF estimates.

## Design

**Data.** Following the example in the maintainer's `Untitled1.R`, each unit's dose-response curve includes the heterogeneity term

$$(2X_1 - 1)\{(A^2 - 1) + .3A^3\} + .5 X_2 A,$$

which averages to 0 at every treatment value ($X_1 \sim \mathrm{Bernoulli}(.5)$ independent of $A$, $E[X_2] = 0$). Each unit's curve is therefore cubic and differs by covariate, while the ADRF has only the shape added to the outcome: flat, linear ($\delta a$), quadratic ($\delta a^2$), or cubic ($\delta a^3$). The outcome is $Y = g(A) + \text{heterogeneity} + X_1 + X_2 + 3\varepsilon$, and the outcome model, `Y ~ poly(A, 3) * (X1 + X2)`, is correctly specified. The effect sizes $\delta$ were calibrated in a pilot for roughly 50% power of the flat ADRF test at n = 1000 and scaled by $\sqrt{1000/n}$ so power is comparable at n = 300.

**Treatment distributions.** *Uniform*: $A \sim U(-2, 2)$, independent of the covariates, with the grid on $(-1.9, 1.9)$ as in `Untitled1.R`. *Skewed*: $A = 4\,\mathrm{Beta}(2, 4) - 2$, confounded by $X_2$ (through a Gaussian copula with correlation .6), with the grid between the population .025 and .975 quantiles, $(-1.79, .87)$. The grid has 21 points.

**Covariance.** The oracle: the covariance of the ADRF estimates across 4000 independent replicates of each cell (`run-oracle.R`). Each test was also computed with the package's default HC2 unconditional covariance as a secondary lens; those results are similar and are summarized at the end.

**Tests.** All are quadratic forms in the ADRF estimates on the grid, with p-values from Imhof's method, except where noted. *Uniform* weighting uses the trapezoidal weights of the grid (what the package does now); *empirical* weighting uses the share of the sample's treatment values in each grid point's bin, so the integral is over the sample's treatment distribution.

- Flat ADRF (`adrf_flat`): weighted squared deviations of the ADRF from its weighted mean.
- Reference curve equal to 0 (`ref_*`), with the reference at the 0, .25, .5, .75, and 1 quantiles of the treatment within the grid.
- AMEF equal to 0 (`amef_zero`).
- Westling-type primitive tests: the primitive of the centered ADRF under the empirical treatment distribution, $\Omega(a_k) = \sum_{j \le k} \hat F_j (\hat\theta_j - \hat\gamma)$, tested by its empirically weighted integrated square (`prim_l2`) and by its maximum absolute value (`prim_max`, p-value by simulation). These use Westling's (2022) primitive applied to the g-computation ADRF, not his one-step estimator.
- Wald test of flatness (`wald_flat`), on 3 df.
- For linearity: linear ADRF (`adrf_linear`), flat AMEF (`amef_flat`), and the Wald test on 2 df (`wald_linear`).

1000 replicates per cell (Monte Carlo SE about .007 at .05), for 2 treatment distributions × 4 shapes × n = 300 and 1000.

## Results (oracle covariance)

**Size.** Every test held its size: the largest rejection rate of a true null across the cells was .067 for the flat-family tests and .069 for the linear-family tests, and most were within two Monte Carlo SEs of .05. The differences below are therefore differences in power alone.

**Tests of no effect (flat ADRF), power at n = 1000:**

| Test | Uniform: linear | quadratic | cubic | Skewed: linear | quadratic | cubic |
|---|---|---|---|---|---|---|
| Flat ADRF, uniform weights | .525 | .520 | .631 | .502 | .495 | .362 |
| Flat ADRF, empirical weights | .525 | .522 | .629 | .547 | .635 | .469 |
| Reference = 0, uniform weights (range over 5 references) | .25 to .45 | .25 to .65 | .31 to .59 | .20 to .60 | .06 to .59 | .09 to .53 |
| AMEF = 0, empirical weights | .078 | .375 | .359 | .132 | .351 | .346 |
| Primitive, integrated square | .654 | .144 | .594 | .631 | .676 | .463 |
| Primitive, maximum | .596 | .217 | .504 | .603 | .685 | .471 |
| Wald | .532 | .516 | .661 | .471 | .577 | .413 |

Summarized over the three alternatives at n = 1000 (mean and minimum power):

| Test | Uniform: mean | min | Skewed: mean | min |
|---|---|---|---|---|
| Flat ADRF, uniform weights | .559 | .520 | .453 | .362 |
| **Flat ADRF, empirical weights** | **.559** | **.522** | **.550** | **.469** |
| Reference = 0, best and worst reference | .483, .329 | .370, .249 | .468, .112 | .352, .055 |
| AMEF = 0, empirical weights | .271 | .078 | .276 | .132 |
| Primitive, integrated square | .464 | .144 | .590 | .463 |
| Primitive, maximum | .439 | .217 | .586 | .471 |
| Wald | .570 | .516 | .487 | .413 |

**Tests of linearity, power at n = 1000:**

| Test | Uniform: quadratic | cubic | Skewed: quadratic | cubic |
|---|---|---|---|---|
| Linear ADRF, uniform weights | .615 | .157 | .314 | .174 |
| **Linear ADRF, empirical weights** | **.614** | **.156** | **.289** | **.222** |
| Flat AMEF, uniform weights | .420 | .220 | .202 | .150 |
| Flat AMEF, empirical weights | .426 | .219 | .194 | .193 |
| Wald | .575 | .179 | .248 | .180 |

**Disagreement.** Pooled over all cells, the uniform-weighted flat ADRF test reached a different decision at .05 than the empirically weighted one in 6.7% of replicates, than Wald in 7.1%, than the primitive tests in 17% to 18%, than the reference-curve tests in 13% to 25%, and than the AMEF test in 24% to 25%. For linearity, the linear ADRF test disagreed with its empirically weighted version in 2.8%, with Wald in 3.4%, and with the flat AMEF test in 8% to 9%.

## Reading the results

- **Testing the ADRF itself beats testing a derived curve.** The flat ADRF test had higher power than the AMEF test against every alternative, by a factor of 1.4 to 7 (the largest against a linear effect, .525 against .078), and than the reference-curve test for most references; the linear ADRF test had higher power than the flat AMEF test in 6 of the 8 comparisons (all but the cubic alternative with the uniform treatment, .156 against .219). Differentiating the curve amplifies the noise at the ends of the grid, and a reference adds the noise of a single point.
- **The reference-curve test depends on the reference in ways a user cannot anticipate.** With the skewed treatment, a reference at the sparse upper end gave 6% power against a quadratic effect (essentially none), and a reference at the .75 quantile 69%. The same reference can be among the best against one alternative and the worst against another: with the uniform treatment, a reference at the .25 quantile had the second-highest power against a linear effect (.43) and the lowest against a quadratic one (.25).
- **Empirical weighting never hurt and often helped.** With the uniform treatment it changed nothing (the weights are nearly equal), and with the skewed treatment it raised the flat test's power from .50 to .55 (linear), .50 to .64 (quadratic), and .36 to .47 (cubic), because uniform weights give the sparse, noisy upper end of the grid as much weight as the dense part. It is also the weighting under which the test statistic is the average squared departure over the population's treatment values.
- **The primitive tests have a blind spot.** They were the most powerful against monotone departures and against the skewed alternatives, but against a quadratic effect symmetric about the middle of the grid their power fell to .14 to .22, because the positive and negative parts of the centered ADRF cancel in the integral. Westling notes that tests with this construction must have low power against some smooth alternatives.
- **The Wald test is competitive under the oracle** (mean power .570 and .487, against .559 and .550 for the empirically weighted flat test), but it is invariant only because it weights departures by their precision, and with an estimated robust covariance it over-rejects more than the omnibus test (`_dev/omnibus-size-study/README.md`).

## Recommendation

**Test hypotheses on the ADRF itself, with the deviations weighted by the empirical distribution of the treatment.** The empirically weighted ADRF test had the best or near-best mean power and the best or near-best worst-case power for both hypotheses in both treatment distributions, held its size, and has no blind spot or arbitrary choice. Concretely:

- For "does the treatment have an effect", test whether the ADRF is flat (`summary(adrf, "flat")`), rather than whether a reference curve or the AMEF is 0.
- For "is the effect linear", test whether the ADRF is linear (`summary(adrf, "linear")`), rather than whether the AMEF is flat.
- Reference curves and AMEFs remain the right objects for estimating and plotting effects relative to a reference or slopes; their omnibus tests are valid but less powerful, and the reference-curve test's power depends on the reference.

Adopting empirical weighting would mean changing the weights in `summary.effect_curve()` (and `curve_projection()`, which uses the same trapezoidal weights) from the trapezoidal weights of the grid to the sample's treatment distribution on the grid, which needs the treatment values (or their distribution on the grid) to be stored in the `effect_curve` object. That is a design decision and has not been implemented.

**Implemented 2026-11-16** as the `integrate` argument of `summary.effect_curve()` and `curve_projection()` (`"uniform"`, the default, or `"empirical"`), with the density estimated by `stats::density()` from all units used, as planned in `_dev/plan-density-weights.md`.

## Questions about empirical weighting

**Does the randomness of the empirical weights need to be accounted for?** Not to first order. Under the null hypothesis the deviations of the ADRF from flat are $O_p(n^{-1/2})$, so the statistic is $O_p(n^{-1})$, and the weights enter it only multiplied by squared deviations: perturbing the weights by $O_p(n^{-1/2})$ changes the statistic by a relative $O_p(n^{-1/2})$, which vanishes. The simulation recomputed the weights in every replicate and treated them as fixed, and the empirically weighted tests held their size (largest .061). The package already treats a data-dependent quantity this way: the default grid runs between sample quantiles of the treatment (`range = .95`), which are fixed once computed.

**Does it change the hypothesis?** Not the null hypothesis: with positive weight at every grid point, the statistic is 0 if and only if the ADRF is flat (or linear) at every grid point, whatever the weights. It changes what the statistic measures under the alternative: with empirical weights it estimates $\mathrm{Var}_F\{\theta(A)\}$, the variance of the ADRF over the population's treatment values within the range, rather than the variance of the ADRF over the range as if the treatment were uniform. Two caveats follow. A grid point whose bin holds no units gets weight 0, so the test says nothing about the ADRF there; a smoothed density (e.g., a kernel density estimate) in place of bin counts would avoid that. And for projections (`curve_projection()`), the weights define the projection itself (the best approximation in $L^2(F)$ rather than in $L^2$ over the range), so they would change its coefficients, not only its test.

**Doesn't the null distribution already account for the less certain regions?** It accounts for them in the size, not in the power. A noisy region contributes large eigenvalues to the null distribution of the statistic, which raises the critical value; a departure in a precise region then has to clear that higher bar, so the noise of one region dilutes the evidence from another. Down-weighting the noisy region in the statistic itself avoids that dilution. The Wald statistic does this fully, by weighting with the inverse covariance; empirical weighting does it approximately and without inverting anything, because the variance of the ADRF at a treatment value is roughly inversely proportional to the density of the treatment there. With the skewed treatment, the variance of the ADRF varied 19-fold across the grid, the empirical weights tracked its inverse (correlation .93 on the log scale), and the coefficient of variation of the weighted variances fell from .88 with uniform weights to .46. With the uniform treatment, the variance still varied 6-fold (it rises at the ends of the range, where the polynomial is least constrained), which empirical weights, nearly constant there, do not correct.

## HC2 covariance

With the package's HC2 covariance in place of the oracle, the pattern is the same (mean power for the empirically weighted flat test .559 and .555, uniform .557 and .465, Wald .577 and .505). Sizes were larger in a few cells, mostly with the skewed treatment: up to .080 for the reference at the upper end and .077 for Wald, against up to .066 for the flat ADRF tests (.063 with empirical weights) and up to .078 for the linear ADRF tests (.074 with empirical weights), which reflects the remaining small-sample error of HC2 described in `_dev/omnibus-size-study/README.md`.

## Limitations

One outcome model family (a cubic polynomial with interactions, correctly specified), two treatment distributions, polynomial alternatives, and one effect size per alternative calibrated for about 50% power. A misspecified model, nonpolynomial alternatives (e.g., a threshold), or a binary outcome were not studied.

## Files

`common.R` (data, tests, run helpers), `run-oracle.R` (stage 1, oracle covariances, `oracle-cov.rds`), `run-tests.R` (stage 2, `tests.rds`), `analyze.R` (`summary.rds`).

## Reference

Westling, T. (2022). Nonparametric tests of the causal null with nondiscrete exposures. *Journal of the American Statistical Association*, 117(539), 1551–1562. https://doi.org/10.1080/01621459.2020.1865168
