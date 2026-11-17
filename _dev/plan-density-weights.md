# Plan: weighting omnibus tests and projections by the treatment density (2026-11-16)

Implements the recommendation of `_dev/omnibus-test-choice/README.md`: an option to weight the grid in omnibus tests and projections by a kernel density estimate of the treatment, in place of the trapezoidal (uniform) weights.

**Decided and implemented (2026-11-16):** one grid (no separate internal range); the argument is `integrate`, with values `"uniform"` (the default) and `"empirical"`. The tests are in `tests/testthat/test-integrate.R`.

## Where things stand

- The grid weights are computed in two places, both from the grid alone: `summary.effect_curve()` (`w <- get_trapezoidal_w(.values)`, one vector for every subgroup) and `curve_projection()` (the same, repeated for each subgroup). `anova()` compares projections, which carry their weights implicitly.
- The `effect_curve` object stores the grid (`.values`), estimates, covariance or draws, subgroups (`.by_grid`), and contrasts, but nothing about the distribution of the treatment. `adrf()` has the treatment values (`treat_var`, stacked across imputations), the target weights (`wts.list`), `subset`, and the subgroup codes, and discards them after fitting.
- `range` is applied once, in `adrf()`, to build the grid; every later function treats the grid as the domain of the curve (evaluation and the rule against values outside the grid, plots, `summary.curve_est()`, tests, and projections).

## The range question

The proposal is that the density weighting be over the original distribution of the treatment, which would mean estimating the curve over the full range (`range = 1`) internally, even when it is displayed over a narrower range such as the default `.95`.

**Recommendation: keep one grid, compute the density from all of the original treatment values, and do not separate an internal range from a displayed one.** The density would be estimated from every unit's treatment value (in the subgroup, weighted by `wts`), not only those within the range, so the weights are the original density evaluated on the grid. The test then integrates over the original distribution restricted to the range (the middle 95% of its mass by default); a user who wants it over the full distribution sets `range = 1`.

Reasons:

- **Two domains would give tests two meanings.** With a full internal grid, uniformly weighted tests would either cover the full range, where the sparse extremes dominate (the reason `range = .95` is the default), or the displayed range, so the domain of a test would depend on the weighting. The same would apply to projections, `summary.curve_est()`, and `amef()`.
- **The extremes are what the default range excludes on purpose.** Beyond the .025 and .975 quantiles, the ADRF depends on few units and, with confounding, on extrapolation over the covariates. Density weights make their contribution to the statistic small, but their estimates still come from the model's extrapolation.
- **The difference is small.** With density weights, the excluded tails carry 5% of the weight by construction, and the curve there is the least precise.
- **Resolution.** `n` grid points spread over the full range leave fewer in the displayed range (about 30% fewer for `logBLL`), which matters for interpolation and plots.
- **The change would touch every function.** Evaluation, the outside-the-grid rule, plotting, printing, `summary.curve_est()`, simultaneous bands, reference curves, AMEFs, contrasts, and projections all use `.values` as the domain.

If the full range is wanted with a narrower display, a lighter alternative is a display argument to `plot()` (and perhaps `print()`), so `adrf(range = 1)` can be plotted over the middle 95%.

## Design

**Argument.** `summary.effect_curve()` and `curve_projection()` gain `weights = "uniform"` (or another name, since "weights" may be read as sampling weights; alternatives are `grid_weights` or `integrate`). Allowable values: `"uniform"` (the trapezoidal weights, as now) and `"density"` (trapezoidal weights times the estimated treatment density at each grid point, normalized), matched with `arg::match_arg()` so that further options (e.g., the GPS-based weights used in *discoverHTE*, or a user-supplied function of the treatment values) can be added later without changing the interface.

**What `adrf()` stores.** A new attribute, `.treat_density`: the kernel density estimate of the treatment for each subgroup, evaluated on a fine grid (e.g., 512 points from `stats::density()`, covering the range of the treatment), so that weights can be computed at the grid points and at any added point (e.g., the reference value added by `reference_curve()`) by interpolation, without storing the treatment values.

- Estimated with `stats::density()`, using all treatment values of the units in `subset` and in the subgroup (not only those within `range`), weighted by the target weights `wts` (normalized), with the bandwidth from Scott's rule (`bw.nrd()`, computed from `bw.nrd0()`; changed from `bw.nrd0()` on 2026-11-17).
- With multiply imputed data, from the treatment values of all imputed datasets, each weighted by $1/m$.
- Each subgroup's total target weight is stored too, for contrasts.
- With bootstrapping or Bayesian models, the density is from the original sample and treated as fixed, as the grid already is.

**Propagation.** `.make_fun()` copies `.treat_density` like the other attributes, so `amef()`, `reference_curve()`, `curve_projection()`, and the subgroup `subset` argument inherit it. For `curve_contrast()`, the density of each contrast is the pooled density of the two subgroups, weighted by their total target weights.

**Use.** A helper `get_grid_w(values, weights, density, b)` returns the trapezoidal or density weights for subgroup (or contrast) `b`. `summary.effect_curve()` computes them per subgroup (currently one vector for all subgroups). `curve_projection()` uses them in the projection fit and its influence (`.nls_fit()`, `.proj_influence()`) and stores the setting, so `anova()` can refuse projections with different weights and `print()` and `summary()` can report it.

**Objects created before this change** have no `.treat_density`; `weights = "density"` then gives an error asking to rerun `adrf()`.

## Documentation

- `?summary.effect_curve` and `?curve_projection`: the argument, and a Details section explaining that the null hypothesis is the same under either weighting but the statistic measures departures over the population's treatment values (the variance of the ADRF over the treatment distribution) rather than uniformly over the range, that sparse regions, whose estimates are least precise, get less weight, that for projections the weights define the projection itself, and the recommendation from the study, with the alternatives' power.
- `?adrf`: that the density is stored, and that `range = 1` with density weights tests over the full distribution.
- Vignette: a short subsection.
- NEWS.

## Tests

- `"uniform"` reproduces current results exactly.
- `"density"` weights match a hand computation with `stats::density()` (weighted, by subgroup) and the trapezoidal widths, and the test matches a hand-computed quadratic form with Imhof.
- With a uniform treatment, the two weightings nearly agree.
- Subgroups get their own densities; a contrast gets the pooled density.
- A reference curve's added grid point gets an interpolated density.
- A density-weighted projection equals weighted least squares on the grid with those weights (and the NLS fit with a transformation); its standard errors match a numerical Jacobian.
- `anova()` refuses projections with different weights.
- Multiply imputed data pool the densities; bootstrap and Bayesian paths run.
- An object without `.treat_density` gives the error.

Afterward, rerun a few cells of the test-choice study with the package's implementation (kernel density in place of the bin counts used there) to confirm the power is unchanged.

## Decisions needed

1. Keep one grid as recommended, or separate an internal full range from the displayed range.
2. The argument's name and its value names (`"uniform"`, `"density"`, `"empirical"`).
3. The default: `"uniform"` keeps current results (including projection coefficients, which the vignette reports); `"density"` follows the study's recommendation for tests but also changes the projections, which `summary()` uses for formula and named hypotheses, so the two should share a default.
