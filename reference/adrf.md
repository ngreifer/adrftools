# Estimate an average dose-response function (ADRF)

Estimates the average dose-response function (ADRF) for a fitted model
object.

## Usage

``` r
adrf(x, ...)

# Default S3 method
adrf(
  x,
  treat,
  vcov = "unconditional",
  cluster = NULL,
  type = "response",
  data = NULL,
  subset = NULL,
  by = NULL,
  wts = NULL,
  range = 0.95,
  n = 51,
  fwb.args = list(),
  ...
)
```

## Arguments

- x:

  a fitted model object (e.g., from
  [`lm()`](https://rdrr.io/r/stats/lm.html) or
  [`glm()`](https://rdrr.io/r/stats/glm.html)).

- ...:

  further arguments passed to
  [`marginaleffects::get_predict()`](https://rdrr.io/pkg/marginaleffects/man/get_predict.html).

- treat:

  a string specifying the name of the treatment variable.

- vcov:

  how the covariance matrix of the estimates should be computed. If
  `"unconditional"` (the default for frequentist models), use the
  sandwich estimator including sampling uncertainty, with the HC2
  small-sample adjustment when it is available (see Details). If
  `"boot"` or `"fwb"`, use the traditional or fractional weighted
  bootstrap, respectively (both of which require the
  [fwb](https://CRAN.R-project.org/package=fwb) package to be
  installed). Otherwise, may be a covariance matrix or other allowed
  input to the `vcov` argument of
  [`marginaleffects::get_vcov()`](https://rdrr.io/pkg/marginaleffects/man/get_vcov.html).
  Can also be `"none"` to avoid computing the uncertainty. For Bayesian
  models, only `"posterior"`, which uses the posterior of the estimates,
  and `"none"` are allowed. For models fit to multiply imputed data,
  `"boot"` and `"fwb"` are not allowed.

- cluster:

  an optional data frame or one-sided formula with the clustering terms
  for cluster-robust inference. When `vcov` is `"unconditional"`, the
  CR2 small-sample adjustment is used when it is available (see
  Details). When `vcov` is `"conditional"`, the covariance matrix of the
  model coefficients is computed using
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  with its default type; `vcov` can also be `"HC0"`, `"HC1"`, `"HC2"`,
  or `"HC3"` to request that type from
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  (`"HC2"` and `"HC3"` are available only for
  [`lm()`](https://rdrr.io/r/stats/lm.html) and
  [`glm()`](https://rdrr.io/r/stats/glm.html) models). Other strings and
  functions supplied to `vcov` cannot be used with `cluster`.

- type:

  character string indicating the type of prediction. Passed to
  [`marginaleffects::get_predict()`](https://rdrr.io/pkg/marginaleffects/man/get_predict.html).
  Default is `"response"` for predictions on the scale of the outcome
  variable. Other options might include `"link"` for the linear
  predictor. This argument is ignored for `lm` objects.

- data:

  an optional data frame containing the observations originally used to
  fit the outcome model supplied to `x`. This is needed when the dataset
  cannot be recovered from the model (e.g., when the model was fit
  inside a function or is not supported by insight), and it is used to
  find variables named in `subset`, `by`, `wts`, and `cluster` that are
  not in the dataset used to fit the model. In most cases, this does not
  need to be supplied.

- subset:

  an optional logical expression indicating the subset of data to use
  for estimation. Will be evaluated in the environment of the original
  dataset supplied to the model fitting function. Can also be a vector
  of row numbers or row names; logical vectors and row numbers refer to
  the units used to fit the model (i.e., excluding any dropped because
  of missing values), and row names refer to the original dataset.

- by:

  optional variable(s) over which to group the estimation. Can be a
  character vector or one-sided formula.

- wts:

  optional weights used to generalize the effect curve to a weighted
  population. Can be a numeric vector, the name of a variable in the
  original dataset, `TRUE` to use the weights used to fit the model, or
  `FALSE` for no weights. The default depends on the model; see Details.

- range:

  numeric; a numeric vector corresponding either to the lower and upper
  bounds of the treatment values for which to compute the effect curve
  or a single number corresponding to the middle quantile of the
  treatment. Default is .95 to use the .025 and .975 quantiles of the
  treatment. See Details.

- n:

  integer specifying the number of equally spaced grid points on which
  to compute the effect curve anchor points. Must be between 2 and 1000.
  Default is 51; higher numbers increase computation time and size of
  the resulting object but improve accuracy.

- fwb.args:

  an optional list of arguments to be passed to
  [`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.html) when
  `vcov` is `"boot"` or `"fwb"`.

## Value

An object of class `<effect_curve>`. This object is a function with
attributes. See
[`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
for details on this function and its outputs.

## Details

`adrf()` estimates the ADRF by computing average predicted outcomes in
the sample for counterfactual treatment values, optionally stratified by
grouping variables and accounting for estimation uncertainty via
unconditional or conditional variance estimation or bootstrapping.
Unconditional variance estimation and bootstrapping treat the sample as
random. When `vcov = "unconditional"`, the variance is computed using
the formulas in Hansen et al. (2024), which involves augmenting the
influence function with a term to account for sampling from the
superpopulation. Unconditional variance estimation requires
[`sandwich::estfun()`](https://zeileis.codeberg.page/sandwich/reference/estfun.html)
and
[`sandwich::bread()`](https://zeileis.codeberg.page/sandwich/reference/bread.html)
methods for the supplied object to be available. For models fit with
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) or
[`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) , these are
computed from the penalized fit instead, treating the smoothing
parameters as fixed, so the variance does not account for their
selection, which can make it noticeably smaller than a bootstrap
variance that does; models with extended families (e.g.,
[`mgcv::nb()`](https://rdrr.io/pkg/mgcv/man/negbin.html)) are not
supported.

In the unconditional variance, each unit's contribution to the influence
function of the model coefficients is divided by \\\sqrt{1 - h_i}\\,
where \\h_i\\ is the unit's leverage (i.e., hat value), as in the HC2
estimator of MacKinnon and White (1985). Without this adjustment, the
variance would share the downward bias of the HC0 sandwich estimator,
which is largest for estimates that depend on units with high leverage,
such as the ADRF at the extremes of the treatment, and which can make
the omnibus tests in
[`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
reject a true null hypothesis too often in small samples. The leverages
are computed from the fitted model for
[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) models, taken from the fit
for [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) models (in
which they are the leverages of the penalized fit, which sum to the
effective degrees of freedom), and extracted using
[`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html) for
other models; when they cannot be found, the HC1 adjustment is used
instead with a message, which multiplies the variance contributed by the
model coefficients by \\n/(n-p)\\, where \\n\\ is the number of units
and \\p\\ the number of estimated parameters (the effective degrees of
freedom for [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) and
[`mgcv::bam()`](https://rdrr.io/pkg/mgcv/man/bam.html) models), and
corrects the bias only for the average leverage. The adjustment is not
used with `svyglm` models, whose variance is computed from the survey
design.

When `cluster` is supplied, the unconditional variance uses the CR2
adjustment described by Pustejovsky and Tipton (2018), the cluster-level
analog of HC2: within each cluster, the vector of residuals is
multiplied by \\(I - H\_{gg})^{-1/2}\\, where \\H\_{gg}\\ is the block
of the hat matrix for the units in the cluster, in place of the usual
\\G/(G-1)\\ adjustment, where \\G\\ is the number of clusters. With
clusters of one unit, it is the HC2 adjustment. CR2 is available only
for [`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) models (including those
whose influence function is that of an ordinary
[`glm()`](https://rdrr.io/r/stats/glm.html) fit, such as
`glm_weightit()` models with fixed weights); for other models, the CR1
adjustment, \\\frac{G}{G-1}\frac{n-1}{n-p}\\ (the default of
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)),
is used with a message. For generalized linear models, the residuals and
hat matrix are those of the final iteratively reweighted least squares
fit, so the adjusted variance can differ from that of
`sandwich::vcovCL(type = "HC2")` when the working weights vary within
clusters.

When a `<mira>` object from mice or a `<mimira>` object from MatchThem
is supplied, analyses are applied to each imputed dataset and pooled
using Rubin's rules. The degrees of freedom used for tests and
confidence intervals are those of Reiter (2007), or of Li, Raghunathan,
and Rubin (1991) when the complete-data degrees of freedom are infinite
(e.g., for generalized linear models). These account for the number of
imputations and the fraction of missing information, which is assumed to
be the same for all points on the curve, so a single value is used for
every point and for joint inference. Bootstrapping is not allowed with
such objects.

### `range`

The `range` argument controls for which range of the treatment the
effect curve is to be evaluated. It can be supplied either as two
numbers corresponding to the lower and upper bounds for the treatment
(e.g., `range = c(0, 10)`) or as a single number corresponding to the
middle quantile of the treatment (e.g., `range = .9`, which uses the .05
and .95 quantiles of the treatment as the bounds). The default is .95 to
use the .025 and .975 quantiles of the treatment. When supplied as a
quantile, the quantiles are evaluated incorporating the weights supplied
to `wts`.

A reason not to use the full treatment range (e.g., by setting
`range = 1`) is that there is likely very little certainty about the
effect curve at the treatment extremes. This uncertainty can muddy tests
of the effect curve. However, limiting the treatment range means
inferences about the effect curve are less generalizable to more extreme
values of the treatment. Note that this does not change the data used to
fit the effect curve, just the points along the effect curve for which
inference and estimation are to take place.

### `wts`

The `wts` argument supplies weights that define the population over
which the predicted outcomes are averaged, such as sampling weights.
When supplied as a numeric vector, it must have one value for each unit
used to fit the model (or, with multiply imputed data, either one value
for each unit in each imputed dataset or one value for each row of the
imputed datasets stacked together). When supplied as a string, it is the
name of a variable in the original dataset. `wts = TRUE` uses the
weights used to fit the model, as extracted by
[`insight::get_weights()`](https://easystats.github.io/insight/reference/get_weights.html),
and `wts = FALSE` weights all units equally.

When `wts` is not supplied, the default depends on the model. When a
`<svyglm>` object from survey is supplied, the survey weights are used.
When a `<glm_weightit>` object, etc., from WeightIt is supplied, the
sampling weights are used if `s.weights` were supplied in the original
call to `weightit()`; the balancing weights are not, because they are
used to fit the outcome model rather than to define the target
population. For other models, all units are weighted equally, even when
the model was fit with weights. See
[`vignette("adrftools")`](https://ngreifer.github.io/adrftools/articles/adrftools.md)
for more details on using the `wts` argument.

## References

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample
significance levels from multiply imputed data using moment-based
statistics and an *F* reference distribution. *Journal of the American
Statistical Association*, 86(416), 1065–1073.
[doi:10.1080/01621459.1991.10475152](https://doi.org/10.1080/01621459.1991.10475152)

MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
*Journal of Econometrics*, 29(3), 305–325.
[doi:10.1016/0304-4076(85)90158-7](https://doi.org/10.1016/0304-4076%2885%2990158-7)

Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics*, 36(4),
672–683.
[doi:10.1080/07350015.2016.1247004](https://doi.org/10.1080/07350015.2016.1247004)

Reiter, J. P. (2007). Small-sample degrees of freedom for
multi-component significance tests with multiple imputation for missing
data. *Biometrika*, 94(2), 502–508.
[doi:10.1093/biomet/asm028](https://doi.org/10.1093/biomet/asm028)

## See also

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting the ADRF

- [`summary.effect_curve()`](https://ngreifer.github.io/adrftools/reference/summary.effect_curve.md)
  for testing hypotheses about the ADRF

- [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md)
  for computing point estimates along the ADRF

- [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  for projecting a simpler model onto the ADRF

- [`reference_curve()`](https://ngreifer.github.io/adrftools/reference/reference_curve.md)
  for computing the difference between each point on the ADRF and a
  specific reference point

- [`curve_contrast()`](https://ngreifer.github.io/adrftools/reference/curve_contrast.md)
  for contrasting ADRFs computed within subgroups

- [`amef()`](https://ngreifer.github.io/adrftools/reference/amef.md) for
  computing the average marginal effect function (AMEF), the derivative
  of the ADRF

- [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md)
  for estimating the mean outcome when every unit's treatment is shifted
  from its observed value

- [`marginaleffects::avg_predictions()`](https://rdrr.io/pkg/marginaleffects/man/predictions.html)
  for computing average adjusted predictions for fitted models (similar
  to the ADRF)

## Examples

``` r
data("nhanes3lead")

fit <- lm(Math ~ poly(logBLL, 5) *
            Male * (Age + Race + PIR +
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

## Plot the ADRF
plot(adrf1)


## ADRF estimates at given points
adrf1(logBLL = c(0, 1, 2)) |>
  summary()
#>               ADRF Estimates
#> ──────────────────────────────────────────
#>  logBLL Estimate Std. Error CI Low CI High
#>       0    8.413     0.2312  7.861   8.964
#>       1    8.005     0.1005  7.765   8.245
#>       2    7.103     0.2495  6.508   7.698
#> ──────────────────────────────────────────
#> Inference: unconditional, simultaneous
#> Confidence level: 95% (t* = 2.385, df = 2437)

# ADRF of logBLL on Math, unconditional
# inference; manual range
adrf2 <- adrf(fit, treat = "logBLL",
              range = c(0, 2))

adrf2
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: 0 to 2
#>  - inference: unconditional
#> 

plot(adrf2)

# ADRF of logBLL on Math, bootstrap
# inference
adrf_b <- adrf(fit, treat = "logBLL",
               vcov = "fwb")
#> 
  |                                                  | 0 % ~calculating  
  |+                                                 | 1 % ~31s          
  |+                                                 | 2 % ~23s          
  |++                                                | 3 % ~21s          
  |++                                                | 4 % ~19s          
  |+++                                               | 5 % ~18s          
  |+++                                               | 6 % ~18s          
  |++++                                              | 7 % ~19s          
  |++++                                              | 8 % ~20s          
  |+++++                                             | 9 % ~19s          
  |+++++                                             | 10% ~19s          
  |++++++                                            | 11% ~18s          
  |++++++                                            | 12% ~18s          
  |+++++++                                           | 13% ~17s          
  |+++++++                                           | 14% ~17s          
  |++++++++                                          | 15% ~16s          
  |++++++++                                          | 16% ~16s          
  |+++++++++                                         | 17% ~16s          
  |+++++++++                                         | 18% ~15s          
  |++++++++++                                        | 19% ~15s          
  |++++++++++                                        | 20% ~15s          
  |+++++++++++                                       | 21% ~14s          
  |+++++++++++                                       | 22% ~14s          
  |++++++++++++                                      | 23% ~14s          
  |++++++++++++                                      | 24% ~14s          
  |+++++++++++++                                     | 25% ~13s          
  |+++++++++++++                                     | 26% ~13s          
  |++++++++++++++                                    | 27% ~13s          
  |++++++++++++++                                    | 28% ~13s          
  |+++++++++++++++                                   | 29% ~13s          
  |+++++++++++++++                                   | 30% ~12s          
  |++++++++++++++++                                  | 31% ~12s          
  |++++++++++++++++                                  | 32% ~12s          
  |+++++++++++++++++                                 | 33% ~12s          
  |+++++++++++++++++                                 | 34% ~11s          
  |++++++++++++++++++                                | 35% ~11s          
  |++++++++++++++++++                                | 36% ~11s          
  |+++++++++++++++++++                               | 37% ~11s          
  |+++++++++++++++++++                               | 38% ~11s          
  |++++++++++++++++++++                              | 39% ~11s          
  |++++++++++++++++++++                              | 40% ~10s          
  |+++++++++++++++++++++                             | 41% ~10s          
  |+++++++++++++++++++++                             | 42% ~10s          
  |++++++++++++++++++++++                            | 43% ~10s          
  |++++++++++++++++++++++                            | 44% ~10s          
  |+++++++++++++++++++++++                           | 45% ~10s          
  |+++++++++++++++++++++++                           | 46% ~09s          
  |++++++++++++++++++++++++                          | 47% ~09s          
  |++++++++++++++++++++++++                          | 48% ~09s          
  |+++++++++++++++++++++++++                         | 49% ~09s          
  |+++++++++++++++++++++++++                         | 50% ~09s          
  |++++++++++++++++++++++++++                        | 51% ~08s          
  |++++++++++++++++++++++++++                        | 52% ~08s          
  |+++++++++++++++++++++++++++                       | 53% ~08s          
  |+++++++++++++++++++++++++++                       | 54% ~08s          
  |++++++++++++++++++++++++++++                      | 55% ~08s          
  |++++++++++++++++++++++++++++                      | 56% ~08s          
  |+++++++++++++++++++++++++++++                     | 57% ~07s          
  |+++++++++++++++++++++++++++++                     | 58% ~07s          
  |++++++++++++++++++++++++++++++                    | 59% ~07s          
  |++++++++++++++++++++++++++++++                    | 60% ~07s          
  |+++++++++++++++++++++++++++++++                   | 61% ~07s          
  |+++++++++++++++++++++++++++++++                   | 62% ~06s          
  |++++++++++++++++++++++++++++++++                  | 63% ~06s          
  |++++++++++++++++++++++++++++++++                  | 64% ~06s          
  |+++++++++++++++++++++++++++++++++                 | 65% ~06s          
  |+++++++++++++++++++++++++++++++++                 | 66% ~06s          
  |++++++++++++++++++++++++++++++++++                | 67% ~06s          
  |++++++++++++++++++++++++++++++++++                | 68% ~06s          
  |+++++++++++++++++++++++++++++++++++               | 69% ~05s          
  |+++++++++++++++++++++++++++++++++++               | 70% ~05s          
  |++++++++++++++++++++++++++++++++++++              | 71% ~05s          
  |++++++++++++++++++++++++++++++++++++              | 72% ~05s          
  |+++++++++++++++++++++++++++++++++++++             | 73% ~05s          
  |+++++++++++++++++++++++++++++++++++++             | 74% ~04s          
  |++++++++++++++++++++++++++++++++++++++            | 75% ~04s          
  |++++++++++++++++++++++++++++++++++++++            | 76% ~04s          
  |+++++++++++++++++++++++++++++++++++++++           | 77% ~04s          
  |+++++++++++++++++++++++++++++++++++++++           | 78% ~04s          
  |++++++++++++++++++++++++++++++++++++++++          | 79% ~04s          
  |++++++++++++++++++++++++++++++++++++++++          | 80% ~03s          
  |+++++++++++++++++++++++++++++++++++++++++         | 81% ~03s          
  |+++++++++++++++++++++++++++++++++++++++++         | 82% ~03s          
  |++++++++++++++++++++++++++++++++++++++++++        | 83% ~03s          
  |++++++++++++++++++++++++++++++++++++++++++        | 84% ~03s          
  |+++++++++++++++++++++++++++++++++++++++++++       | 85% ~03s          
  |+++++++++++++++++++++++++++++++++++++++++++       | 86% ~02s          
  |++++++++++++++++++++++++++++++++++++++++++++      | 87% ~02s          
  |++++++++++++++++++++++++++++++++++++++++++++      | 88% ~02s          
  |+++++++++++++++++++++++++++++++++++++++++++++     | 89% ~02s          
  |+++++++++++++++++++++++++++++++++++++++++++++     | 90% ~02s          
  |++++++++++++++++++++++++++++++++++++++++++++++    | 91% ~02s          
  |++++++++++++++++++++++++++++++++++++++++++++++    | 92% ~01s          
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 93% ~01s          
  |+++++++++++++++++++++++++++++++++++++++++++++++   | 94% ~01s          
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 95% ~01s          
  |++++++++++++++++++++++++++++++++++++++++++++++++  | 96% ~01s          
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 97% ~01s          
  |+++++++++++++++++++++++++++++++++++++++++++++++++ | 98% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 99% ~00s          
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=17s  

adrf_b
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>  - inference: bootstrap
#> 

plot(adrf_b)

# ADRF in subset
adrf_m <- adrf(fit, treat = "logBLL",
               subset = Male == 1)

adrf_m
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>  - inference: unconditional
#> 

# ADRFs in subgroups
adrf_by <- adrf(fit, treat = "logBLL",
                by = ~Male)

adrf_by
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Math
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>  - by: Male
#>  - inference: unconditional
#> 
```
