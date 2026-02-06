# Effect curve objects

An `effect_curve` object is a function that takes in values of the
treatment and produces estimates of the effect curve at those values.
`effect_curve` objects are produces by
[`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) and
functions that modify effect curves, such as
[`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md),
[`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md),
[`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md),
and
[`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md).
The output of an `effect_curve` object is a `curve_est` object
containing the effect curve estimates. This page describes
`effect_curve` and `curve_est` objects.

### Usage

    f <- adrf(x, ...)

    f({treat}, subset = NULL)

## Arguments

- {treat}:

  the values of the treatment at which to evaluate the effect curve.

- subset:

  an optional logical expression indicating the subset of the subgroups
  for which to compute estimates. Can only be used when `by` was
  supplied to the original call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md),
  and only to refer to variables defining subgroups.

- x:

  a `curve_est` object; the output of an `effect_curve` object call.

- digits:

  the number of digits to display.

- ...:

  arguments passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html).

## Value

A call to an `effect_curve` object returns a `curve_est` object, which
is a `data.frame` containing a column for the treatment and a column for
the effect curve estimates. `curve_est` objects have
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md),
[`ceof()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md),
and
[`vcov()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
methods.

## Details

An `effect_curve` object contains a set of grid points on which the
effect curve is initially evaluated. The effect curve estimates produced
by a call to the `effect_curve` object are interpolated using 3rd-degree
local polynomial regression with a Gaussian kernel and bandwidth equal
to half the distance between grid points, unless they coincide with the
grid points; this means the produced estimates are linear combinations
of the grid point estimates.

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  generating an effect curve

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for performing inference on effect curve estimates

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting the effect curve

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

adrf1
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>  - inference: unconditional
#> 

# Compute estimates along the ADRF
adrf1(logBLL = c(0, 1, 2))
#>  ADRF Estimates
#> ────────────────
#>  logBLL Estimate
#>       0    8.470
#>       1    8.012
#>       2    7.029
#> ────────────────

# Perform inference on the estimates
adrf1(logBLL = c(0, 1, 2)) |>
  summary()
#>               ADRF Estimates
#> ──────────────────────────────────────────
#>  logBLL Estimate Std. Error CI Low CI High
#>       0    8.470     0.2196  7.946   8.993
#>       1    8.012     0.0985  7.777   8.247
#>       2    7.029     0.2298  6.482   7.577
#> ──────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (t* = 2.384, df = 2473)

# ADRF within groups defined by `Male`
adrf2 <- adrf(fit, treat = "logBLL",
              by = ~Male)

adrf2
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>  - by: Male
#>  - inference: unconditional
#> 

# Estimates in both groups
adrf2(logBLL = c(0, 1, 2))
#>    ADRF Estimates
#> ─────────────────────
#>  Male logBLL Estimate
#>     0      0    8.367
#>     0      1    8.406
#>     0      2    7.281
#>     1      0    8.569
#>     1      1    7.631
#>     1      2    6.786
#> ─────────────────────

# Estimates in one group
adrf2(logBLL = c(0, 1, 2), subset = Male == 1)
#>    ADRF Estimates
#> ─────────────────────
#>  Male logBLL Estimate
#>     1      0    8.569
#>     1      1    7.631
#>     1      2    6.786
#> ─────────────────────
```
