# Compute the average derivative effect

`ade()` computes the average derivative effect (ADE) from a shift curve,
which is the derivative of the shift curve at a shift of 0.

## Usage

``` r
ade(x, subset = NULL, eps = 1e-05)
```

## Arguments

- x:

  a `shift_curve` object; the output of a call to
  [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md),
  or of
  [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  or
  [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md)
  applied to one.

- subset:

  an optional logical expression indicating the subset of the subgroups
  for which to compute the ADE. Can only be used when `by` was supplied
  to the original call to
  [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md),
  and only to refer to variables defining subgroups.

- eps:

  numeric; the step size to use when calculating the numerical
  derivative, relative to the range of shifts at which the shift curve
  was computed: the step is `eps` times a quarter of that range. Default
  is `1e-5` (.00001).

## Value

An object of class `ade_est`, which is like a
[`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
object but with one estimate per subgroup (or per contrast) and no
column for the shift.
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html), and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) work as they do for
`curve_est` objects; see
[`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md).

## Details

With an additive shift, the ADE is \\E\[\partial \mu(A, X) / \partial
a\]\\, where \\\mu(a, x)\\ is the expected outcome given treatment \\a\\
and covariates \\x\\: the slope of the outcome with respect to the
treatment for each unit at the treatment it actually received, averaged
over units. It is the change in the mean outcome per unit of an
infinitesimal shift of every unit's treatment (Hines, Diaz-Ordaz, &
Vansteelandt, 2026). Because it involves only treatment values near
those actually received, it requires positivity only locally, whereas
the ADRF and its derivative, the AMEF, require that every unit could
have received each treatment value in the range.

The ADE differs from the AMEF averaged over the treatment. The AMEF at a
treatment value \\a\\ averages the slope at \\a\\ over the distribution
of the covariates in the whole sample, whereas the ADE averages each
unit's slope at its own treatment value. The two coincide when the
treatment is unrelated to the covariates (e.g., in a randomized trial)
but differ otherwise. A single-number summary of the ADRF itself is the
slope of its linear projection; see
[`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md).

With `shift_fun = "multiplicative"` in
[`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md),
the ADE is \\E\[A \\ \partial \mu(A, X) / \partial a\]\\, the change in
the mean outcome per unit proportional increase in every unit's
treatment; multiplying it by .01 gives the change for a 1% increase.
More generally, for a function \\g(a, \delta)\\ supplied to `shift_fun`,
the ADE averages each unit's slope weighted by the derivative of \\g\\
with respect to \\\delta\\ at 0. `ade()` can only be used when \\g(a, 0)
= a\\ for every unit, i.e., when a shift of 0 leaves every unit's
treatment unchanged, and when the range of shifts in `x` includes 0.

The ADE is computed as the derivative of the interpolated shift curve at
0, using the same method as
[`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md), and
its standard error is computed from that of the shift curve using the
delta method, or from the bootstrap replicates or posterior draws of the
shift curve when those are used. For smooth outcome models, it agrees
closely with the ADE computed directly from the predictions (e.g., by
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html))
and does not depend on the range of shifts. For outcome models whose
predictions are step functions of the treatment (e.g., Bayesian additive
regression trees with hard decision rules), the derivative of each
unit's prediction is 0 almost everywhere, and the ADE computed from the
shift curve reflects the change in the mean outcome over the spacing of
the grid of shifts, which depends on `shift` and `n` in the call to
[`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md);
in that case, the effects of specific shifts computed using
[`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
are easier to interpret than the ADE. Bayesian additive regression trees
with soft (smooth) decision rules have predictions that are smooth in
the treatment, so the ADE does not depend on the grid.

The default unconditional standard errors account for sampling from the
population, which matters when the slope varies across units. The ADE of
a contrast between subgroups (i.e., when `x` is the output of
[`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md))
is the difference between the ADEs in those subgroups. When `x` contains
shift curves for several subgroups,
[`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
computes simultaneous confidence intervals across the subgroups by
default; set `simultaneous = FALSE` for pointwise intervals.

## References

Hines, O. J., Diaz-Ordaz, K., & Vansteelandt, S. (2026). Parameterizing
the effect of a continuous treatment using average derivative effects.
*Biometrika*, 113(2), asag012.
[doi:10.1093/biomet/asag012](https://doi.org/10.1093/biomet/asag012)

## See also

- [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md)
  for estimating a shift curve

- [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md) for
  computing the derivative of an effect curve

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for performing inference on the ADE

- [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  for summarizing an ADRF with the slope of its linear projection

- [`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
  for computing the ADE with conditional inference

## Examples

``` r
data("nhanes3lead")

fit <- lm(Math ~ poly(logBLL, 5) *
            Male * (Age + Race + PIR +
                      Enough_Food),
          data = nhanes3lead)

# Shift curve of logBLL on Math, unconditional
# inference
sc <- shift_curve(fit, treat = "logBLL")

# ADE of logBLL on Math
ade(sc) |>
  summary()
#>         Average Derivative Effect Estimates
#> ────────────────────────────────────────────────────
#>  Estimate Std. Error       t P-value  CI Low CI High
#>   -0.5335     0.2142 -2.4907  0.0128 -0.9536 -0.1135
#> ────────────────────────────────────────────────────
#> Inference: unconditional, pointwise
#> Confidence level: 95% (t* = 1.961, df = 2437)
#> Null value: 0

# ADEs within subgroups and their difference
sc_by <- shift_curve(fit, treat = "logBLL",
                     by = ~Male)

ade(sc_by) |>
  summary()
#>            Average Derivative Effect Estimates
#> ─────────────────────────────────────────────────────────
#>  Male Estimate Std. Error       t P-value  CI Low CI High
#>     0  -0.5520     0.3553 -1.5539  0.2261 -1.3470  0.2429
#>     1  -0.5157     0.2440 -2.1130  0.0682 -1.0618  0.0305
#> ─────────────────────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (t* = 2.238, df = 2437)
#> Null value: 0

curve_contrast(sc_by) |>
  ade() |>
  summary()
#>                Average Derivative Effect Contrast Estimates
#> ───────────────────────────────────────────────────────────────────────────
#>                 Contrast Estimate Std. Error      t P-value  CI Low CI High
#>  [Male = 1] - [Male = 0]   0.0364      0.431 0.0844  0.9327 -0.8088  0.8816
#> ───────────────────────────────────────────────────────────────────────────
#> Inference: unconditional, pointwise
#> Confidence level: 95% (t* = 1.961, df = 2437)
#> Null value: 0
```
