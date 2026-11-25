# Contrast multiple subgroup effect curves

`curve_contrast()` computes the difference between effect curves across
levels of a subgrouping variable supplied to `by` in the original call
to [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

## Usage

``` r
curve_contrast(x)
```

## Arguments

- x:

  an
  [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  object; the output of a call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)
  with `by` supplied.

## Value

An object of class `contrast_curve`, which inherits from
[`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md),
with additional information about the groups being contrasted.

## Details

`curve_contrast()` creates a new effect curve corresponding to the
difference between effect curves in two groups. When multiple subgroups
are specified by `by` in the original call to
[`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md), all
pairwise comparisons are included. Use the `subset` argument in the
original function call to restrict comparisons to fewer groups.

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  computing the ADRF

- [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  for comparing effect curves to a point along the curve

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting the effect curve contrasts

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for performing tests of effect curve contrasts at specific points

- [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  for performing omnibus tests of effect curve contrasts (e.g., whether
  the contrast curve differs from 0)

## Examples

``` r
data("nhanes3lead")

fit <- lm(Math ~ poly(logBLL, 5) *
            Male * Smoke_in_Home *
            (Age + Race + PIR),
          data = nhanes3lead)

# ADRFs in Race subgroups, excluding Other
adrf_by <- adrf(fit, treat = "logBLL",
                by = ~Race,
                subset = Race != "Other")

adrf_by
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>    + grid points: 51
#>  - by: Race (3 subgroups)
#>  - sample size: 2,410
#>  - inference: unconditional (HC2)
#> 

# Contrast subgroup ADRFs
adrf_contrast <- curve_contrast(adrf_by)

adrf_contrast
#> An <effect_curve> object
#> 
#>  - curve type: ADRF contrast
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>    + grid points: 51
#>  - contrasts: [Race = "Hispanic"] - [Race = "Black"], [Race = "White"] - [Race = "Black"], [Race = "White"] - [Race = "Hispanic"]
#>  - sample size: 2,410
#>  - inference: unconditional (HC2)
#> 

# Plot contrast ADRFs
plot(adrf_contrast, simultaneous = FALSE)


# Compute and test difference between subgroup
# ADRFs at specific points
adrf_contrast(logBLL = c(0, 2)) |>
  summary()
#>                          ADRF Contrast Estimates
#> ─────────────────────────────────────────────────────────────────────────
#>                                Contrast logBLL Estimate Std. Error      t
#>  [Race = "Hispanic"] - [Race = "Black"]      0  -0.3926     0.5210 -0.754
#>  [Race = "Hispanic"] - [Race = "Black"]      2   0.2837     0.3862  0.735
#>     [Race = "White"] - [Race = "Black"]      0   1.7703     0.5556  3.186
#>     [Race = "White"] - [Race = "Black"]      2  -0.3083     0.6811 -0.453
#>  [Race = "White"] - [Race = "Hispanic"]      0   2.1629     0.5409  3.999
#>  [Race = "White"] - [Race = "Hispanic"]      2  -0.5921     0.7006 -0.845
#> ─────────────────────────────────────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (t* = 2.59, df = 2377)
#> Null value: 0

# Test if ADRF differences are present
summary(adrf_contrast)
#>                   Omnibus Curve Test
#> ───────────────────────────────────────────────────────
#> H₀: ADRF contrast is 0 for values of logBLL between
#>     -0.3567 and 2.4248
#> 
#>                                Contrast P-value
#>  [Race = "Hispanic"] - [Race = "Black"]  0.8987
#>     [Race = "White"] - [Race = "Black"]  0.0058
#>  [Race = "White"] - [Race = "Hispanic"]  0.0025
#> ───────────────────────────────────────────────────────
#> Computed using the Imhof approximation
```
