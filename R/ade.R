#' Compute the average derivative effect
#'
#' @description
#' `ade()` computes the average derivative effect (ADE) from a shift curve, which is the derivative of the shift curve at a shift of 0.
#'
#' @param x a `shift_curve` object; the output of a call to [shift_curve()], or of [reference_curve()] or [curve_contrast()] applied to one.
#' @param subset an optional logical expression indicating the subset of the subgroups for which to compute the ADE. Can only be used when `by` was supplied to the original call to [shift_curve()], and only to refer to variables defining subgroups.
#' @param eps numeric; the step size to use when calculating the numerical derivative, relative to the range of shifts at which the shift curve was computed: the step is `eps` times a quarter of that range. Default is `1e-5` (.00001).
#'
#' @returns
#' An object of class `ade_est`, which is like a [`curve_est`] object but with one estimate per subgroup (or per contrast) and no column for the shift. `summary()`, `coef()`, and `vcov()` work as they do for `curve_est` objects; see [summary.curve_est()].
#'
#' @details
#' With an additive shift, the ADE is \eqn{E[\partial \mu(A, X) / \partial a]}, where \eqn{\mu(a, x)} is the expected outcome given treatment \eqn{a} and covariates \eqn{x}: the slope of the outcome with respect to the treatment for each unit at the treatment it actually received, averaged over units. It is the change in the mean outcome per unit of an infinitesimal shift of every unit's treatment (Hines, Diaz-Ordaz, & Vansteelandt, 2026). Because it involves only treatment values near those actually received, it requires positivity only locally, whereas the ADRF and its derivative, the AMEF, require that every unit could have received each treatment value in the range.
#'
#' The ADE differs from the AMEF averaged over the treatment. The AMEF at a treatment value \eqn{a} averages the slope at \eqn{a} over the distribution of the covariates in the whole sample, whereas the ADE averages each unit's slope at its own treatment value. The two coincide when the treatment is unrelated to the covariates (e.g., in a randomized trial) but differ otherwise. A single-number summary of the ADRF itself is the slope of its linear projection; see [curve_projection()].
#'
#' With `shift_fun = "multiplicative"` in [shift_curve()], the ADE is \eqn{E[A \, \partial \mu(A, X) / \partial a]}, the change in the mean outcome per unit proportional increase in every unit's treatment; multiplying it by .01 gives the change for a 1% increase. More generally, for a function \eqn{g(a, \delta)} supplied to `shift_fun`, the ADE averages each unit's slope weighted by the derivative of \eqn{g} with respect to \eqn{\delta} at 0. `ade()` can only be used when \eqn{g(a, 0) = a} for every unit, i.e., when a shift of 0 leaves every unit's treatment unchanged, and when the range of shifts in `x` includes 0.
#'
#' The ADE is computed as the derivative of the interpolated shift curve at 0, using the same method as [amef()], and its standard error is computed from that of the shift curve using the delta method, or from the bootstrap replicates or posterior draws of the shift curve when those are used. For smooth outcome models, it agrees closely with the ADE computed directly from the predictions (e.g., by [marginaleffects::avg_slopes()]) and does not depend on the range of shifts. For outcome models whose predictions are step functions of the treatment (e.g., Bayesian additive regression trees), the derivative of each unit's prediction is 0 almost everywhere, and the ADE computed from the shift curve reflects the change in the mean outcome over the spacing of the grid of shifts, which depends on `shift` and `n` in the call to [shift_curve()].
#'
#' The default unconditional standard errors account for sampling from the population, which matters when the slope varies across units. The ADE of a contrast between subgroups (i.e., when `x` is the output of [curve_contrast()]) is the difference between the ADEs in those subgroups.
#'
#' @seealso
#' * [shift_curve()] for estimating a shift curve
#' * [amef()] for computing the derivative of an effect curve
#' * [summary.curve_est()] for performing inference on the ADE
#' * [curve_projection()] for summarizing an ADRF with the slope of its linear projection
#' * [marginaleffects::avg_slopes()] for computing the ADE with conditional inference
#'
#' @references
#' Hines, O. J., Diaz-Ordaz, K., & Vansteelandt, S. (2026). Parameterizing the effect of a continuous treatment using average derivative effects. *Biometrika*, 113(2), asag012. \doi{10.1093/biomet/asag012}
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             Male * (Age + Race + PIR +
#'                       Enough_Food),
#'           data = nhanes3lead)
#'
#' # Shift curve of logBLL on Math, unconditional
#' # inference
#' sc <- shift_curve(fit, treat = "logBLL")
#'
#' # ADE of logBLL on Math
#' ade(sc) |>
#'   summary()
#'
#' # ADEs within subgroups and their difference
#' sc_by <- shift_curve(fit, treat = "logBLL",
#'                      by = ~Male)
#'
#' ade(sc_by) |>
#'   summary()
#'
#' curve_contrast(sc_by) |>
#'   ade() |>
#'   summary()

#' @export
ade <- function(x, subset = NULL, eps = 1e-5) {
  arg::arg_supplied(x)

  .shift <- .attr(x, ".shift")

  if (is_null(.shift) || !inherits(x, "effect_curve")) {
    arg::err("{.arg x} must be a shift curve, the output of {.fun shift_curve}; for a single-number summary of an ADRF, see {.fun curve_projection} with {.code integrate = \"empirical\"}")
  }

  check_effect_curve(x, amef_ok = FALSE, projection_ok = FALSE)

  if (!isTRUE(.shift[["identity0"]])) {
    arg::err("the ADE can only be computed when a shift of 0 leaves every unit's treatment unchanged, which is not the case for the {.arg shift_fun} used to compute {.arg x}")
  }

  .values <- .attr(x, ".values")

  if (min(.values) > 0 || max(.values) < 0) {
    arg::err("the ADE can only be computed when the range of shifts in {.arg x} includes 0")
  }

  d <- amef(x, eps = eps)

  subset_substitute <- substitute(subset)

  res <- {
    if (anyv(rlang::fn_fmls_names(d), "subset")) {
      rlang::inject(d(shift = 0, subset = !!subset_substitute))
    }
    else if (is_null(subset_substitute)) {
      d(shift = 0)
    }
    else {
      arg::err("{.arg subset} can only be used when {.arg by} was supplied to {.fun shift_curve} and {.arg x} is not a curve contrast")
    }
  }

  res[["shift"]] <- NULL

  attr(res, "values") <- NULL
  attr(res, ".est_names") <- "ADE"

  class(res) <- unique(c("ade_est", class(res)))

  res
}
