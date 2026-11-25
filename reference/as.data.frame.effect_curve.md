# Extract estimates and confidence bands of an effect curve

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the estimates of an effect curve at its grid points and the limits of
their confidence bands as a data frame. These are the values
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) displays, so
the data frame can be used to build a custom plot of the effect curve.

## Usage

``` r
# S3 method for class 'effect_curve'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  conf_level = 0.95,
  simultaneous = TRUE,
  subset = NULL,
  outcome = NULL,
  proj = NULL,
  transform = TRUE,
  ci.type = "perc",
  df = NULL,
  ...
)
```

## Arguments

- x:

  an `effect_curve` object; the output of a call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) or
  functions that modify it.

- row.names:

  `NULL` or a character vector of row names for the returned data frame.
  Default is `NULL` to use the row numbers.

- optional:

  ignored.

- conf_level:

  the desired confidence level. Default is .95 for 95% confidence bands.

- simultaneous:

  `logical`; whether the computed confidence bands should be
  simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known
  as uniform) bands cover the full line at the desired confidence level,
  whereas pointwise confidence bands only cover each point at the
  desired level. Default is `TRUE`. See Details.

- subset:

  an optional logical expression indicating the subset of the subgroups
  to plot. Can only be used when `by` was supplied to the original call
  to [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md),
  and only to refer to variables defining subgroups.

- outcome:

  an optional character vector naming the categories of the outcome to
  plot. Can only be used when the model gives a separate prediction for
  each category of the outcome (e.g., multinomial and ordinal models);
  see the section *Categorical outcomes* in Details of
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- proj:

  an optional `curve_projection` object, the output of a call to
  [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md).
  Supplying this adds the projection curve to the effect curve plot.

- transform:

  whether to compute intervals on the transformed estimates. Allowable
  options include `TRUE`, `FALSE`, a family object (e.g.,
  [`binomial()`](https://rdrr.io/r/stats/family.html)), a link object
  (e.g., the output of
  [`make.link()`](https://rdrr.io/r/stats/make.link.html)), or one of
  the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`. Ignored
  unless `x` is an ADRF and Wald inference is used (i.e., ignored with
  percentile intervals, the default with bootstrapping and Bayesian
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

  the "denominator" degrees of freedom to use for the critical test
  statistics for confidence bands. Default is to use the residual
  degrees of freedom from the original model if it is a linear model and
  `Inf` otherwise, or, for models fit to multiply imputed data, degrees
  of freedom that account for the number of imputations (see
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)).

- ...:

  ignored.

## Value

A data frame with one row for each grid point of the effect curve within
each subgroup or contrast, ordered by subgroup or contrast and then by
treatment value. It has the following columns:

- the variables supplied to `by` in the call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  (when supplied), with one value for each subgroup, or `contrast`, the
  label of each contrast, when `x` is the output of a call to
  [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md).
  For models that give a separate prediction for each category of the
  outcome, a column (usually named `outcome`) contains the category.

- a column named for the treatment (or `shift` for a shift curve)
  containing the grid points.

- `estimate`: the estimate of the effect curve at each grid point.

- `conf.low` and `conf.high`: the lower and upper limits of the
  confidence band at each grid point. These are omitted when
  `vcov = "none"` in the call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- `proj_estimate`: the estimate of the projection at each grid point,
  when `proj` is supplied.

## Details

The confidence bands are computed exactly as they are by
[`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md),
with the same arguments, so a call to
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) and a
call to [`plot()`](https://rdrr.io/r/graphics/plot.default.html) with
the same arguments give the same band (with Wald inference and
`simultaneous = TRUE`, the critical value is computed using a randomized
algorithm, so [`set.seed()`](https://rdrr.io/r/base/Random.html) should
be called before each for the bands to be identical). To compute
estimates and confidence intervals at specific treatment values rather
than at the grid points, call the effect curve with those values and use
[`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
on the result.

## See also

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting an effect curve

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for confidence intervals at specific treatment values

## Examples

``` r
data("nhanes3lead")

fit <- lm(Math ~ poly(logBLL, 3) * Male + Age + Race +
            PIR + NICU + Smoke_Pregnant,
          data = nhanes3lead)

# ADRFs of logBLL on Math within subgroups
adrf1 <- adrf(fit, treat = "logBLL", by = ~Male,
              n = 21)

# Estimates and simultaneous 95% confidence bands
# at the grid points
set.seed(123)
bands <- as.data.frame(adrf1)

head(bands)
#>   Male      logBLL estimate conf.low conf.high
#> 1    0 -0.35667494 8.825445 8.018613  9.632278
#> 2    0 -0.21760106 8.793493 8.167224  9.419762
#> 3    0 -0.07852718 8.754384 8.246975  9.261794
#> 4    0  0.06054671 8.708593 8.269553  9.147633
#> 5    0  0.19962059 8.656594 8.253283  9.059906
#> 6    0  0.33869447 8.598863 8.216720  8.981006

# Pointwise bands for one subgroup, with the linear
# projection of the ADRF
as.data.frame(adrf1, simultaneous = FALSE,
              subset = Male == 1,
              proj = curve_projection(adrf1, "linear")) |>
  head()
#>   Male      logBLL estimate conf.low conf.high proj_estimate
#> 1    1 -0.35667494 8.321380 7.767064  8.875696      8.682132
#> 2    1 -0.21760106 8.344383 7.909623  8.779144      8.571710
#> 3    1 -0.07852718 8.343474 7.979274  8.707675      8.461289
#> 4    1  0.06054671 8.320511 7.990647  8.650374      8.350867
#> 5    1  0.19962059 8.277350 7.962885  8.591815      8.240446
#> 6    1  0.33869447 8.215852 7.911953  8.519751      8.130024
```
