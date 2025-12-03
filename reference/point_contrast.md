# Contrast Point Estimates Along an Effect Curve

`point_contrast()` computes pairwise contrasts of estimates from an
effect curve.

## Usage

``` r
point_contrast(object)

# S3 method for class 'curve_est_contrast'
summary(
  object,
  conf_level = 0.95,
  simultaneous = TRUE,
  null = 0,
  ci.type = "perc",
  df = NULL,
  ...
)
```

## Arguments

- object:

  for `point_contrast()`, a
  [`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  object; the output of a an `effect_curve` object. For
  [`summary()`](https://rdrr.io/r/base/summary.html), a
  `curve_est_contrast` object; the output of a call to
  `point_contrast()`.

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

  the null value for hypothesis tests. Default is 0. Set to `NA` to omit
  tests.

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

`point_contrast()` returns an object of class `curve_est_contrast`,
which is like a `curve_est` object but with its own
[`summary()`](https://rdrr.io/r/base/summary.html) method.

## Details

`point_contrast()` computes all pairwise contrasts between effect curve
estimates. Because pairwise contrasts are a linear operation over the
original estimates, the delta method can be used to perform Wald
inference for the contrasts. When `by` was specified in the original
call to
[`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) or
the effect curve is a `contrast_curve` object resulting from
[`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md),
pairwise contrasts occur only within subgroups or within subgroup
contrasts, respectively. To compare points on an effect curve to a
single point, use
[`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md).

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  computing the ADRF

- [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  for comparing points on an effect curve to a single point

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for inference on individual points on an effect curve

- [`marginaleffects::hypotheses()`](https://marginaleffects.com/man/r/hypotheses.html)
  for general hypotheses on `curve_est` (and other) objects

## Examples

``` r
data("nhanes3lead")

fit <- lm(Math ~ poly(logBLL, 5) *
            (Male + Age + Race + PIR +
               Enough_Food),
          data = nhanes3lead)

# ADRF of logBLL on Math, unconditional
# inference
adrf1 <- adrf(fit, treat = "logBLL")

# Differences among ADRF estimates at given points
adrf1(logBLL = c(0, 1, 2)) |>
  point_contrast() |>
  summary()
#>                           ADRF Point Contrasts
#> ────────────────────────────────────────────────────────────────────────
#>                         Term Estimate Std. Error      t  P-value  CI Low
#>  [logBLL = 1] - [logBLL = 0]  -0.4579     0.2662 -1.720   0.1962 -1.0810
#>  [logBLL = 2] - [logBLL = 0]  -1.4403     0.3074 -4.685 < 0.0001 -2.1601
#>  [logBLL = 2] - [logBLL = 1]  -0.9824     0.2610 -3.764   0.0005 -1.5935
#> ────────────────────────────────────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (t* = 2.341, df = 2473)
#> Null value: 0
```
