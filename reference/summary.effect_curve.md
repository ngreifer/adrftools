# Test omnibus hypotheses about an effect curve

[`summary()`](https://rdrr.io/r/base/summary.html) tests an omnibus
hypothesis about an effect curve. For example, it can be used to test
that the ADRF is flat, that the contrast between two ADRFs is 0
everywhere, or that the AMEF is 0 everywhere.

## Usage

``` r
# S3 method for class 'effect_curve'
summary(
  object,
  hypothesis,
  method,
  subset = NULL,
  outcome = NULL,
  transform = TRUE,
  df = NULL,
  nsim = 1e+06,
  integrate = "uniform",
  ...
)

# S3 method for class 'summary.effect_curve'
print(x, digits = max(4L, getOption("digits") - 3L), ...)
```

## Arguments

- object:

  an object of class
  [`effect_curve`](https://ngreifer.github.io/adrftools/reference/effect_curve.md);
  the result of a call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) or
  a function that modifies it.

- hypothesis:

  the hypothesis to be tested. Allowable options include `"flat"`,
  `"linear"`, `"quadratic"`, `"cubic"`, a one-sided formula
  corresponding to a projection model, or a single number (e.g., 0). See
  Details. The default is `"flat"` for ADRFs and 0 otherwise.

- method:

  string; the method used to compute the p-value of the test. Allowable
  options include `"sim"` for simulation-based, `"imhof"` for Imhof's
  approximation, `"davies"` for Davies's approximation, `"liu"` for
  Liu's approximation, `"satterthwaite"` for Satterthwaite's
  approximation, and `"saddlepoint"` for a saddlepoint approximation.
  Default is `"imhof"` when the CompQuadForm package is installed and
  `"saddlepoint"` otherwise. See Details.

- subset:

  an optional logical expression indicating the subset of the subgroups
  for which to perform tests. Can only be used when `by` was supplied to
  the original call to
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md),
  and only to refer to variables defining subgroups.

- outcome:

  an optional character vector naming the categories of the outcome for
  which to perform tests. Can only be used when the model gives a
  separate prediction for each category of the outcome (e.g.,
  multinomial and ordinal models); see the section *Categorical
  outcomes* in Details of
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md).

- transform:

  whether to perform the test on the transformed estimates. Allowable
  options include `TRUE`, `FALSE`, a family object (e.g.,
  [`binomial()`](https://rdrr.io/r/stats/family.html)), a link object
  (e.g., the output of
  [`make.link()`](https://rdrr.io/r/stats/make.link.html)), or one of
  the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`. Ignored
  unless `object` is an ADRF object. See Details.

- df:

  the "denominator" degrees of freedom to use for the tests. Default is
  to use the residual degrees of freedom from the original model if it
  is a linear model and `Inf` otherwise, or, for models fit to multiply
  imputed data, degrees of freedom that account for the number of
  imputations (see
  [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md)).

- nsim:

  when `method` is `"sim"`, the number of iterations used to simulate
  the p-values. Higher numbers give more accurate p-values subject to
  less Monte Carlo error but are slower and require more memory. Default
  is 1,000,000.

- integrate:

  string; how the deviations of the effect curve from the hypothesized
  curve are weighted across treatment values in the test statistic.
  Allowable options include `"uniform"` (the default), which weights all
  treatment values in the range of the effect curve equally, and
  `"empirical"`, which weights them by an estimate of the density of the
  treatment. Only `"uniform"` can be used with a shift curve (see
  [`shift_curve()`](https://ngreifer.github.io/adrftools/reference/shift_curve.md)).
  See Details.

- ...:

  when `method` is `"imhof"` or `"davies"`, further arguments passed to
  [`CompQuadForm::imhof()`](https://rdrr.io/pkg/CompQuadForm/man/imhof.html)
  or
  [`CompQuadForm::davies()`](https://rdrr.io/pkg/CompQuadForm/man/davies.html)
  , respectively.

- x:

  a `summary.effect_curve` object.

- digits:

  numeric; the number of digits to print.

## Value

An object of class `"summary.effect_curve"`, which is a data.frame with
a column for the p-value and, for stratified effect curves or contrasts
thereof, additional columns identifying the subset to which the p-value
refers.

## Details

[`summary()`](https://rdrr.io/r/base/summary.html) performs an omnibus
test for an effect curve. The hypothesis tested is determined by the
argument to `"hypothesis"`. When supplied as a single number,
[`summary()`](https://rdrr.io/r/base/summary.html) tests whether all
values on the effect curve are equal to that number. When supplied as a
one-sided formula, [`summary()`](https://rdrr.io/r/base/summary.html)
tests whether the projection of the effect curve model onto the model
represented in the formula is sufficient to describe the effect curve.
The test itself tests whether the residuals around the projection are
all equal to 0, incorporating the uncertainty in estimating the effect
curve. See
[`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
for more information on how the curve projection and the uncertainty in
the residuals are computed.

`"flat"` tests whether all values on the curve are equal to each other
(i.e., whether the curve is flat), without specifying what value they
are equal to. This is equivalent to testing whether the variance around
the mean estimate is different from 0 or whether an intercept-only
projection model is sufficient. `"linear"` tests whether the curve is
linear, i.e., whether the residuals around linear projection of the
curve are all 0. `"quadratic"` and `"cubic"` test whether the curve is
quadratic or cubic, respectively, using the same method.

Rejecting the null hypothesis means that the curve is more complicated
than the specified model. For example, rejecting the null hypothesis
that the curve is linear implies that the curve is nonlinear (and,
therefore, not flat either).

The test involves computing a test statistic, specifying its
distribution under the null hypothesis, and computing the p-value using
the complement of the cumulative distribution function of the
distribution evaluated at the test statistic value. The test statistic
depends on `"hypothesis"`. For `hypothesis` equal to a constant, say,
\\c\\, the test statistic is \$\$T^\* = \int\_\mathcal{A} (\theta(a) -
c)^2 \\ da\$\$ Otherwise, the test statistic is \$\$T^\* =
\int\_\mathcal{A} (\theta(a) - \hat{\theta}\_0(a))^2 \\ da\$\$ where
\\\hat{\theta}\_0\\ is the projection of \\\theta(a)\\ onto the null
subspace specified by `hypothesis`.

Each of these can be approximated as a quadratic form,
\\T=\mathbf{\theta}'\mathbf{W}\mathbf{\theta}\\ where
\\\mathbf{\theta}\\ is a vector of estimates at evaluation points along
the curve and \\\mathbf{W}\\ is a diagonal matrix of weights
implementing a trapezoidal approximation to the integral. The null
hypothesis is that \\T=0\\, which approximates the null hypothesis that
\\T^\*=0\\. Each of the allowable options to `method` corresponds to a
method of approximating the distribution of \\T\\ under the null
hypothesis:

- `"sim"` simulates the null distribution by simulating from a
  multivariate normal distribution under the null hypothesis and
  computing the test statistic in each simulation. The p-value is the
  proportion of simulated estimates greater than the observed test
  statistic. When `df` is not `Inf`, the simulation is done from a
  multivariate t-distribution.

- `"imhof"`, `"davies"`, and `"liu"` assume the test statistic follows a
  generalized \\\chi^2\\-distribution and approximate its CDF
  numerically. `"imhof"` tends to be the most accurate and is
  recommended. Because Imhof's and Davies's methods control the absolute
  rather than the relative error of the computed probability, p-values
  below .0001, and any that these methods fail to compute, are computed
  using the saddlepoint approximation instead. `"liu"` can be inaccurate
  when `df` is not `Inf`. These methods require the
  [CompQuadForm](https://CRAN.R-project.org/package=CompQuadForm)
  package to be installed.

- `"satterthwaite"` also assumes the test statistic follows a
  generalized \\\chi^2\\-distribution, but this distribution is
  approximated using a scaled F-distribution with the same first two
  moments.

- `"saddlepoint"` also assumes the test statistic follows a generalized
  \\\chi^2\\-distribution, and this distribution is approximated using
  the saddlepoint method described by Kuonen (1999), which remains
  accurate far into the tails of the distribution. Near the center of
  the distribution, where the saddlepoint approximation is not defined,
  the Satterthwaite approximation is used instead.

In general, we recommend using `method = "imhof"`, though this requires
CompQuadForm to be installed (and is the default when it is).
`method = "saddlepoint"` is also accurate, is very fast, and doesn't
require any other packages. When using `"sim"`, increasing `nsim`
improves the accuracy of the p-value by reducing Monte Carlo error. The
default value of `1e6` ensures that the simulated p-value is within
.0005 of the true p-value with at least 97.8% confidence for p-values
less than .05. For greater precision at higher p-value thresholds,
`nsim` should be increased. The other methods exist primarily for
comparison purposes. `"satterthwaite"` is particularly fast and doesn't
require any other packages, but it can be somewhat inaccurate.

### Transform

When the effect curve is an ADRF and the outcomes are bounded (e.g., a
probability between 0 and 1), the `transform` argument can be specified,
which changes the details of the tests.

The tests above assume the estimates along the effect curve are normally
distributed (or t-distributed when `df` is not `Inf`). However, when the
outcome is bounded (e.g., a probability bounded between 0 and 1), this
assumption may not be valid for the ADRF in finite samples. How
`transform` is used depends on `hypothesis`. When `hypothesis` is a
number, the estimates and that number are transformed to be on an
unbounded scale, and the test is performed on the transformed estimates,
whose distribution is computed using the delta method (or from the
transformed bootstrap or posterior draws). When `hypothesis` is a string
or formula, the estimates are not transformed; instead, the projection
is fit with the inverse of the transformation as its inverse link, as
described in
[`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md),
and the residuals of the estimates around it are tested. By default, if
a generalized linear model is used for the outcome with a non-identity
link function, the estimates are transformed by the link function to be
on an unbounded scale. Note this is not the same as using the linear
predictor for the effect curve; this is simply a transformation of the
estimated points along the curve already computed. A value of
`hypothesis` outside the domain of the transformation (e.g.,
`hypothesis = 0` with a logit transformation) gives an error; set
`transform = FALSE` to perform the test on the original scale. A
transformation supplied to `transform` that cannot be applied to every
estimate or draw also gives an error, and when `transform = TRUE`, the
estimates are left untransformed with a warning in that case.

Tests on the transformed and untransformed ADRFs correspond to different
hypotheses; the difference is not simply a matter of the appropriate
distribution of the statistic. For example, for a binary outcome model
with a logistic transformation, testing that the transformed ADRF is
linear corresponds to testing whether the ADRF is a sigmoid function,
whereas testing that the untransformed ADRF is linear corresponds to
testing whether the ADRF is a straight line. These choices correspond to
how the projection of the ADRF is formed; see
[`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
for details. See Examples for an example with a binary outcome. In the
example, there is evidence to reject that the transformed ADRF is linear
for one of the groups, indicating that a sigmoid function is not
sufficient for describing the ADRF, but there is indeterminate evidence
(at the .05 level) to reject that the untransformed ADRF is linear for
either group, indicating that a linear function is sufficient for
describing the ADRF (at least in the treatment range examined).

### Weighting treatment values

The test statistic integrates the squared deviations of the effect curve
from the hypothesized curve over the range of treatment values at which
the effect curve was computed, and `integrate` controls how treatment
values are weighted in that integral. With `integrate = "uniform"`,
every treatment value in the range counts equally. With
`integrate = "empirical"`, the integral is weighted by the density of
the treatment, so departures from the null hypothesis count in
proportion to how common the treatment values at which they occur are;
for example, with `hypothesis = "flat"`, the statistic is the variance
of the ADRF over the distribution of the treatment. The null hypothesis
is the same under either option; only the test statistic, and therefore
the power against different alternatives, differs. Because the effect
curve is estimated least precisely where the treatment is sparse,
`integrate = "empirical"` gives the least precise estimates the least
weight, which can increase power when the treatment is unevenly
distributed; when it is evenly distributed, the two options give nearly
the same result.

The density is estimated using
[`density()`](https://rdrr.io/r/stats/density.html) with the bandwidth
given by Scott's rule (i.e.,
[`bw.nrd()`](https://rdrr.io/r/stats/bandwidth.html)) from the treatment
values of all units used to compute the effect curve (including those
outside `range`), weighted by `wts`, separately within each subgroup
defined by `by`; a contrast between subgroups uses the density of the
two subgroups pooled. The density is treated as known. Because the
integral is over the range of the effect curve, with the default `range`
of .95 in
[`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md), the
test concerns the middle 95% of the distribution of the treatment;
setting `range = 1` makes it concern the whole distribution.

## References

Kuonen, D. (1999). Saddlepoint approximations for distributions of
quadratic forms in normal variables. *Biometrika*, 86(4), 929–935.
[doi:10.1093/biomet/86.4.929](https://doi.org/10.1093/biomet/86.4.929)

## See also

- [`adrf()`](https://ngreifer.github.io/adrftools/reference/adrf.md) for
  computing the ADRF

- [`summary.curve_est()`](https://ngreifer.github.io/adrftools/reference/summary.curve_est.md)
  for performing inference on individual points on an effect curve

- [`plot.effect_curve()`](https://ngreifer.github.io/adrftools/reference/plot.effect_curve.md)
  for plotting the effect curve

- [`curve_projection()`](https://ngreifer.github.io/adrftools/reference/curve_projection.md)
  for projecting a simpler model onto an effect curve

## Examples

``` r
data("nhanes3lead")

fit <- glm(Block >= 12 ~ poly(logBLL, 3) *
             Male * (Age + Race + PIR + NICU +
                       Smoke_Pregnant),
           data = nhanes3lead,
           family = binomial)

# ADRFs of logBLL on P(Block >= 12) within
# groups defined by `Male`
adrf1 <- adrf(fit, treat = "logBLL",
              by = ~Male)

adrf1
#> An <effect_curve> object
#> 
#>  - curve type: ADRF
#>  - response: Block >= 12
#>  - treatment: logBLL
#>    + range: -0.3567 to 2.4248
#>    + grid points: 51
#>  - by: Male (2 subgroups)
#>  - sample size: 2,521
#>  - inference: unconditional (HC2)
#> 

# Test if ADRFs are flat
summary(adrf1)
#>                   Omnibus Curve Test
#> ───────────────────────────────────────────────────────
#> H₀: ADRF is flat (transformed) for values of logBLL
#>     between -0.3567 and 2.4248
#> 
#>  Male  P-value
#>     0   0.0001
#>     1 < 0.0001
#> ───────────────────────────────────────────────────────
#> Computed using the Imhof approximation
#>   * Saddlepoint approximation used for p-values below .0001

# Test if logit-transformed ADRFs are linear
# (i.e., if ADRFs have sigmoid shape)
summary(adrf1, hypothesis = "linear")
#>                   Omnibus Curve Test
#> ───────────────────────────────────────────────────────
#> H₀: ADRF is linear (transformed) for values of logBLL
#>     between -0.3567 and 2.4248
#> 
#>  Male P-value
#>     0  0.1827
#>     1  0.0242
#> ───────────────────────────────────────────────────────
#> Computed using the Imhof approximation
# summary(adrf1, hypothesis = ~logBLL) # same test

proj1 <- curve_projection(adrf1, "linear")

plot(adrf1, proj = proj1, conf_level = 0)


# Test if un-transformed ADRFs are linear
summary(adrf1, hypothesis = "linear",
        transform = FALSE)
#>                   Omnibus Curve Test
#> ───────────────────────────────────────────────────────
#> H₀: ADRF is linear for values of logBLL between
#>     -0.3567 and 2.4248
#> 
#>  Male P-value
#>     0  0.4955
#>     1  0.1046
#> ───────────────────────────────────────────────────────
#> Computed using the Imhof approximation

proj2 <- curve_projection(adrf1, "linear",
                          transform = FALSE)

plot(adrf1, proj = proj2, conf_level = 0)


## Test if ADRFs differ from each other by
## testing if the ADRF contrast is 0
curve_contrast(adrf1) |>
  summary()
#>                   Omnibus Curve Test
#> ───────────────────────────────────────────────────────
#> H₀: ADRF contrast is 0 for values of logBLL between
#>     -0.3567 and 2.4248
#> 
#>                 Contrast P-value
#>  [Male = 1] - [Male = 0]  0.0002
#> ───────────────────────────────────────────────────────
#> Computed using the Imhof approximation
```
