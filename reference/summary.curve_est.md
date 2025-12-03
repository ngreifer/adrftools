# Compute Points on an Effect Curve

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

  an
  [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  object; the output of a call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) or
  a function that modifies it.

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
  estimates. Allowable options include `TRUE`, `FALSE`, or a function
  specifying a transformation. Ignored unless `object` is an ADRF. See
  Details.

- ci.type:

  string; when bootstrapping or Bayesian inference is used in the
  original effect curve, which type of confidence interval is to be
  computed. For bootstrapping, allowable options include `"perc"` for
  percentile intervals, `"wald"` for Wald intervals, and other options
  allowed by
  [`fwb::summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.html).
  When `simultaneous = TRUE`, only `"perc"` and `"wald"` are allowed.
  For Bayesian models, allowable options include `"perc"` for
  equi-tailed intervals and `"wald"` for Wald intervals. Default is
  `"perc"`. Ignored when bootstrapping is not used and the model is not
  Bayesian.

- df:

  the "denominator" degrees of freedom to use for the tests and critical
  test statistics for confidence intervals. Default is to use the
  residual degrees of freedom from the original model if it is a linear
  model and `Inf` otherwise.

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
effect curve; this is simple a transformation of the estimated points
along the curve already computed. Confidence intervals are computed
using the transformed estimates before being back-transformed to ensure
they are within the bounds of the outcome. When `null` is a number, that
number is also transformed. When `transform` is specified, standard
errors are not reported (i.e., because the standard errors used for
tests and confidence intervals are those of the transformed estimates).

### Simultaneous confidence intervals and tests

Simultaneous confidence intervals ensure all estimates, not just a given
individual point, are contained within the union of confidence intervals
at the given confidence level. These are wider than pointwise intervals
to reflect that they are covering multiple estimates, which otherwise
would decrease the true coverage rate from that specified.
[`summary()`](https://rdrr.io/r/base/summary.html) uses the "sup-t"
simultaneous confidence interval, which is the smallest one-parameter
interval that covers all estimates at the desired rate. Simultaneous
hypothesis tests are performed by inverting the simultaneous confidence
intervals; the p-value for each test is the complement of the smallest
confidence level for which a simultaneous confidence intervals
accounting for other tests contains the null value. The widths of the
confidence intervals and the p-values depend on how many and which
estimates are computed.

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  computing the ADRF

- [`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  for information on the output of an effect curve

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting an effect curve

- [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  for testing omnibus hypotheses about a effect curve

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
#>       0   0.2274 0.1867  0.2740
#>       1   0.1925 0.1672  0.2206
#>       2   0.0739 0.0386  0.1368
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
#>       0   0.2274 0.1936  0.2653
#>       1   0.1925 0.1715  0.2154
#>       2   0.0739 0.0434  0.1230
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
#>       0   0.2274     0.0183 0.1838  0.2711
#>       1   0.1925     0.0112 0.1658  0.2192
#>       2   0.0739     0.0197 0.0269  0.1209
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
#>       0   0.2274  9.349 < 0.0001 0.1867  0.2740
#>       1   0.1925 10.599 < 0.0001 0.1672  0.2206
#>       2   0.0739 -1.151   0.5741 0.0386  0.1368
#> ───────────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (z* = 2.386)
#> Null value: -2.19722457733622
```
