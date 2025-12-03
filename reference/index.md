# Package index

## Compute the ADRF

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) :
  Estimate an average dose-response function (ADRF)

## Characterize and test an effect curve

- [`plot(`*`<effect_curve>`*`)`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  : Plot an effect curve
- [`summary(`*`<effect_curve>`*`)`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  [`print(`*`<summary.effect_curve>`*`)`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  : Test omnibus hypotheses about an effect curves
- [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  [`summary(`*`<curve_projection>`*`)`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  [`coef(`*`<curve_projection>`*`)`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  [`vcov(`*`<curve_projection>`*`)`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  [`anova(`*`<curve_projection>`*`)`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  : Project an effect curve onto a simpler model

## Compute and contrast points on an effect curve

- [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  [`effect_curve-class`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  [`curve_est`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  [`curve_est-class`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  : Effect curve objects
- [`summary(`*`<curve_est>`*`)`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  [`coef(`*`<curve_est>`*`)`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  [`vcov(`*`<curve_est>`*`)`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  : Compute points on an effect curve
- [`point_contrast()`](https://ngreifer.github.io/adrftools/reference/point_contrast.md)
  [`summary(`*`<curve_est_contrast>`*`)`](https://ngreifer.github.io/adrftools/reference/point_contrast.md)
  : Contrast point estimates along an effect curve

## Compute other effect curves

- [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md) :
  Estimate the average marginal effect function (AMEF)
- [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  : Contrast an effect curve with a reference point
- [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md)
  : Contrast multiple subgroup effect curves

## Example data

- [`nhanes3lead`](https://ngreifer.github.io/adrftools/reference/nhanes3lead.md)
  : Data from NHANES III on blood lead levels and cognitive outcomes
