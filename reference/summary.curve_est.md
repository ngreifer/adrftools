# Compute points on an effect curve

[`summary()`](https://rdrr.io/r/base/summary.html) computes estimates
and confidence intervals for specified points on the supplied effect
curve.

## Usage

``` r
# S3 method for class 'curve_est'
summary(
  object,
  conf_level = 0.95,
  simultaneous = TRUE,
  null = NULL,
  transform = TRUE,
  ci.type = "perc",
  df = NULL,
  ...
)

# S3 method for class 'curve_est'
coef(object, ...)

# S3 method for class 'curve_est'
vcov(object, ...)
```

## Arguments

- object:

  a
  [`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  object; the output of a call to an
  [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  object.

- conf_level:

  the desired confidence level. Set to 0 to omit confidence intervals.
  Default is .95.

- simultaneous:

  `logical`; whether the computed p-values and confidence intervals
  should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous
  (also known as uniform) intervals jointly cover all specified
  estimates at the desired confidence level, whereas pointwise
  confidence intervals only cover each estimate at the desired level.
  Simultaneous p-values are inversions of the simultaneous confidence
  intervals. Default is `TRUE`. See Details.

- null:

  the null value for the hypothesis tests. Default is to use a null
  value of 0 when the effect curve is an AMEF, a curve contrast, or a
  reference effect curve, and to omit hypothesis tests otherwise. Set to
  `NA` to manually omit hypothesis tests.

- transform:

  whether to compute intervals and perform tests on the transformed
  estimates. Allowable options include `TRUE`, `FALSE`, a family object
  (e.g., [`binomial()`](https://rdrr.io/r/stats/family.html)), a link
  object (e.g., the output of
  [`make.link()`](https://rdrr.io/r/stats/make.link.html)), or one of
  the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`. Ignored
  unless `object` is an ADRF and Wald inference is used (i.e., ignored
  with percentile intervals, the default with bootstrapping and Bayesian
  models). See Details.

- ci.type:

  string; when bootstrapping or Bayesian inference is used in the
  original effect curve, which type of confidence interval is to be
  computed. For bootstrapping, allowable options include `"perc"` for
  percentile intervals, `"wald"` for Wald intervals, and other options
  allowed by
  [`fwb::summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.html)
  . When `simultaneous = TRUE`, only `"perc"` and `"wald"` are allowed.
  For Bayesian models, allowable options include `"perc"` for
  equi-tailed intervals and `"wald"` for Wald intervals. Default is
  `"perc"`. Ignored when bootstrapping is not used and the model is not
  Bayesian.

- df:

  the "denominator" degrees of freedom to use for the tests and critical
  test statistics for confidence intervals. Default is to use the
  residual degrees of freedom from the original model if it is a linear
  model and `Inf` otherwise, or, for models fit to multiply imputed
  data, degrees of freedom that account for the number of imputations
  (see
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)).

- ...:

  ignored.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.curve_est`, which inherits from
[`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md).
This is a `data.frame` with columns for the treatment, estimates, and
uncertainty measures (p-values, confidence intervals, etc.).

## Details

### Transform

The usual confidence intervals and tests assume the estimates along the
effect curve are normally distributed (or t-distributed when `df` is not
`Inf`). However, when the outcome is bounded (e.g., a probability
bounded between 0 and 1), this assumption may not be valid for the ADRF
in finite samples. `transform` transforms the estimates to ones that are
unbounded and computes the corresponding distribution of transformed
estimates using the delta method. By default, if a generalized linear
model is used for the outcome with a non-identity link function, the
estimates are transformed by the link function to be on an unbounded
scale. Note this is not the same as using the linear predictor for the
effect curve; this is simply a transformation of the estimated points
along the curve already computed. Confidence intervals are computed
using the transformed estimates before being back-transformed to ensure
they are within the bounds of the outcome. When `null` is a number, that
number is also transformed. When a transformation other than the
identity is used, including the default link function of a generalized
linear model, standard errors are not reported, because the standard
errors used for tests and confidence intervals are those of the
transformed estimates. `transform` is only used with Wald inference;
percentile intervals, the default with bootstrapping and Bayesian
models, already respect the bounds of the outcome because they are
invariant to transformation.

When `transform` is `TRUE`, the link function of the outcome model is
used; when a family or link object is supplied, its link function is
used. Only transformations with a known inverse and derivative are
allowed. A value of `null` outside the domain of the transformation
(e.g., a null value of 0 or 1 with a logit transformation) gives an
error; set `transform = FALSE` to test it on the original scale. A
transformation supplied to `transform` that cannot be applied to every
estimate, or, with Wald intervals, to every bootstrap or posterior draw,
also gives an error. When `transform = TRUE` and the link function of
the outcome model cannot be applied to every draw (e.g., a draw between
grid points just outside 0 and 1), the estimates are left untransformed
with a warning.

### Simultaneous confidence intervals and tests

Simultaneous confidence intervals ensure all estimates, not just a given
individual point, are contained within their confidence intervals at the
given confidence level. These are wider than pointwise intervals to
reflect that they are covering multiple estimates, which otherwise would
decrease the true coverage rate from that specified.
[`summary()`](https://rdrr.io/r/base/summary.html) uses the "sup-t"
simultaneous confidence interval, which is the smallest one-parameter
interval that covers all estimates at the desired rate. Simultaneous
hypothesis tests are performed by inverting the simultaneous confidence
intervals; the p-value for each test is the complement of the smallest
confidence level for which a simultaneous confidence interval accounting
for the other tests contains the null value. The widths of the
confidence intervals and the p-values depend on how many and which
estimates are computed. With Wald inference, the simultaneous critical
value and p-values are computed from the multivariate t (or normal)
distribution of the estimates using
[`mvtnorm::qmvt()`](https://rdrr.io/pkg/mvtnorm/man/qmvt.html) and
[`mvtnorm::pmvt()`](https://rdrr.io/pkg/mvtnorm/man/pmvt.html) , which
use a randomized algorithm, so they can differ slightly between calls;
use [`set.seed()`](https://rdrr.io/r/base/Random.html) for reproducible
results. With bootstrapping or Bayesian models and `ci.type = "perc"`,
the simultaneous intervals are the equal-tailed percentile intervals
with the tail probability calibrated so that they jointly contain the
desired proportion of the bootstrap or posterior draws, as described by
Montiel Olea and Plagborg-Møller (2019).

## References

Montiel Olea, J. L., & Plagborg-Møller, M. (2019). Simultaneous
confidence bands: Theory, implementation, and an application to SVARs.
*Journal of Applied Econometrics*, 34(1), 1–17.
[doi:10.1002/jae.2656](https://doi.org/10.1002/jae.2656)

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  computing the ADRF

- [`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  for information on the output of an effect curve

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting an effect curve

- [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  for testing omnibus hypotheses about an effect curve

## Examples

``` r
data("nhanes3lead")

fit <- glm(Block >= 12 ~ poly(logBLL, 3) *
             Male * (Age + Race + PIR + NICU +
                       Smoke_Pregnant),
           data = nhanes3lead,
           family = binomial)

# ADRF of logBLL on P(Block >= 12)
adrf1 <- adrf(fit, treat = "logBLL")

# Estimates along ADRF with simultaneous CIs computed
# from transformed estimates
adrf1(logBLL = c(0, 1, 2)) |>
  summary()
#>         ADRF Estimates
#> ───────────────────────────────
#>  logBLL Estimate CI Low CI High
#>       0   0.2274 0.1858  0.2752
#>       1   0.1925 0.1669  0.2210
#>       2   0.0739 0.0377  0.1399
#> ───────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (z* = 2.386)

# Estimates along ADRF with pointwise CIs computed
# from transformed estimates
adrf1(logBLL = c(0, 1, 2)) |>
  summary(simultaneous = FALSE)
#>         ADRF Estimates
#> ───────────────────────────────
#>  logBLL Estimate CI Low CI High
#>       0   0.2274 0.1928  0.2662
#>       1   0.1925 0.1712  0.2157
#>       2   0.0739 0.0425  0.1253
#> ───────────────────────────────
#> Inference: unconditional, pointwise
#> Confidence level: 95% (z* = 1.96)

# Estimates along ADRF with simultaneous CIs computed
# from original estimates
adrf1(logBLL = c(0, 1, 2)) |>
  summary(transform = FALSE)
#>               ADRF Estimates
#> ──────────────────────────────────────────
#>  logBLL Estimate Std. Error CI Low CI High
#>       0   0.2274     0.0187 0.1827  0.2721
#>       1   0.1925     0.0113 0.1654  0.2196
#>       2   0.0739     0.0204 0.0251  0.1226
#> ──────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (z* = 2.386)

# Estimates along ADRF with simultaneous CIs computed
# from transformed estimates, hypothesis tests against
# null of .1
adrf1(logBLL = c(0, 1, 2)) |>
  summary(null = .1)
#>                 ADRF Estimates
#> ───────────────────────────────────────────────
#>  logBLL Estimate      z  P-value CI Low CI High
#>       0   0.2274  9.138 < 0.0001 0.1858  0.2752
#>       1   0.1925 10.457 < 0.0001 0.1669  0.2210
#>       2   0.0739 -1.109   0.6029 0.0377  0.1399
#> ───────────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (z* = 2.385)
#> Null value: 0.1
```
