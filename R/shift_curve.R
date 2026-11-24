#' Estimate a shift curve
#'
#' @description
#' Estimates the shift curve for a fitted model object: the mean outcome that would be observed if every unit's treatment were shifted from its observed value, as a function of the size of the shift.
#'
#' @inheritParams adrf
#' @param shift numeric; the sizes of the shift at which to compute the shift curve. Can be two numbers corresponding to the lower and upper bounds of the shift or a single number `s` corresponding to the range from `-s` to `s`. When `shift_fun = "additive"`, a single number is in standard deviations of the treatment; otherwise, it is in the units of the shift. Default is .5 when `shift_fun = "additive"`, i.e., shifts of up to half a standard deviation of the treatment in either direction; `shift` must be supplied otherwise. See Details.
#' @param shift_fun how each unit's treatment is to be shifted. Allowable options include `"additive"` (the default), which adds the shift to the treatment, `"multiplicative"`, which multiplies the treatment by 1 plus the shift (e.g., a shift of .1 is a 10% increase), and a function of two arguments, the vector of observed treatment values and a single value of the shift, that returns the vector of shifted treatment values. `"multiplicative"` can only be used when the treatment is nonnegative, and the shifts cannot be below -1. See Details.
#'
#' @returns
#' An object of class `shift_curve`, which inherits from `adrf_curve` and [`effect_curve`]. Its first argument is named `shift` rather than after the treatment, so estimates are requested as, e.g., `f(shift = c(-.2, 0, .2))`. See [`effect_curve-class`] for details on this function and its outputs.
#'
#' @details
#' `shift_curve()` estimates \eqn{\psi(\delta) = E[Y(g(A, \delta))]}, the mean potential outcome when each unit's treatment \eqn{A} is replaced by \eqn{g(A, \delta)}, where \eqn{g} is the function supplied to `shift_fun` and \eqn{\delta} is the shift. With `shift_fun = "additive"`, \eqn{g(a, \delta) = a + \delta}, so \eqn{\psi(\delta)} is the mean outcome if every unit's treatment were increased by \eqn{\delta} (or decreased, when \eqn{\delta} is negative); \eqn{\psi(0)} is the mean outcome under the treatment actually received. The estimate at each shift is the average of the predicted outcomes computed with each unit's treatment set to its shifted value, in the same way that [adrf()] averages predicted outcomes with every unit's treatment set to the same value. Such interventions, which depend on the treatment each unit would have received, are known as modified treatment policies or shift interventions (Díaz Muñoz & van der Laan, 2012; Haneuse & Rotnitzky, 2013).
#'
#' Whereas the ADRF at a given treatment value requires that every unit could have received that value, a shift curve only requires that units could have received treatment values near those they actually received, so it involves less extrapolation from the outcome model when the treatment is strongly related to the covariates (Lundberg & Brand, 2026). This holds only for small shifts; a large shift moves some units to treatment values that units like them rarely receive, in which case the shift curve depends on extrapolation just as the ADRF does. The derivative of the shift curve at 0 is the average derivative effect (ADE), which can be computed using [ade()].
#'
#' ## `shift` and `shift_fun`
#'
#' When `shift` is a single number and `shift_fun = "additive"`, the shifts are in standard deviations of the treatment, computed over the units in `subset`, weighted by `wts`, and pooled over the subgroups defined by `by` (and over the imputed datasets for multiply imputed data), so that every subgroup has the same grid of shifts and their shift curves can be contrasted. The grid itself is in the units of the treatment and is treated as fixed in computing standard errors, as the range of the ADRF is in [adrf()]. With the default `n` (odd) and a symmetric range, a shift of 0 is one of the grid points.
#'
#' A function supplied to `shift_fun` is called once for each grid value of the shift with the observed treatment values of all units in `subset` (in all subgroups defined by `by` together), so it can be used to define other modified treatment policies. For example, `shift_fun = function(a, delta) pmin(a, delta)` caps every unit's treatment at `delta`, so the shift curve gives the mean outcome as a function of the cap, and `shift` gives the range of caps to consider. The function must return a finite treatment value for every unit at every shift, and the outcome model must be able to make predictions at those values. The function must also be deterministic (i.e., it must give the same result each time it is called with the same inputs). When it depends on the distribution of the observed treatment values (e.g., shifting each unit to a quantile of the treatment), that dependence is treated as fixed in computing standard errors. [ade()] can only be used when a shift of 0 leaves every unit's treatment unchanged, which is not the case for a cap.
#'
#' ## Inference
#'
#' The options for `vcov` are the same as for [adrf()]. Unconditional standard errors account for the sampling of both the covariates and the treatment, since the shifted treatment of each unit depends on the treatment it received.
#'
#' Once computed, a shift curve can be used with the other functions that take an [`effect_curve`] object. [reference_curve()] with `reference = 0` gives the differences between the mean outcome under each shift and under the treatment actually received (i.e., the effect of each shift), [summary.effect_curve()] tests whether the shift curve is flat (i.e., whether no shift in the range changes the mean outcome), and [curve_contrast()] compares shift curves computed within subgroups.
#'
#' @seealso
#' * [ade()] for computing the average derivative effect from a shift curve
#' * [adrf()] for estimating the average dose-response function
#' * [reference_curve()] for computing the effect of each shift relative to no shift
#' * [plot.effect_curve()] for plotting the shift curve
#' * [summary.effect_curve()] for testing hypotheses about the shift curve
#' * [marginaleffects::avg_comparisons()] for computing the effect of a single shift with conditional inference
#'
#' @references
#' Díaz Muñoz, I., & van der Laan, M. (2012). Population intervention causal effects based on stochastic interventions. *Biometrics*, 68(2), 541–549. \doi{10.1111/j.1541-0420.2011.01685.x}
#'
#' Haneuse, S., & Rotnitzky, A. (2013). Estimation of the effect of interventions that modify the received treatment. *Statistics in Medicine*, 32(30), 5260–5277. \doi{10.1002/sim.5907}
#'
#' Lundberg, I., & Brand, J. E. (2026). Causal inference with a continuous treatment: Addressing positivity constraints, nonlinearity, and effect heterogeneity. *Sociological Methodology*, 56(2), 206–236. \doi{10.1177/00811750261459364}
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             Male * (Age + Race + PIR +
#'                       Enough_Food),
#'           data = nhanes3lead)
#'
#' # Mean Math score when every unit's logBLL is
#' # shifted by up to half a standard deviation
#' sc <- shift_curve(fit, treat = "logBLL")
#'
#' sc
#'
#' plot(sc)
#'
#' # Effects of shifts relative to no shift
#' reference_curve(sc, reference = 0)(shift = c(-.3, .3)) |>
#'   summary()
#'
#' # Average derivative effect
#' ade(sc) |>
#'   summary()
#'
#' # Mean Math score when logBLL is capped at
#' # each value from 0 to 2
#' sc_cap <- shift_curve(fit, treat = "logBLL",
#'                       shift = c(0, 2),
#'                       shift_fun = function(a, delta) pmin(a, delta))
#'
#' plot(sc_cap)

#' @export
shift_curve <- function(x, ...) {
  UseMethod("shift_curve")
}

#' @exportS3Method shift_curve default
#' @rdname shift_curve
shift_curve.default <- function(x, treat, shift = NULL, shift_fun = "additive",
                                vcov = "unconditional", cluster = NULL, type = NULL,
                                outcome = NULL, data = NULL, subset = NULL, by = NULL, wts = NULL,
                                n = 51, fwb.args = list(), ...) {

  arg::arg_supplied(treat)

  .check_misplaced_args(...names(), "shift_curve")

  subset_substitute <- substitute(subset)

  shift_fun <- process_shift_fun(shift_fun,
                                 label = .shift_fun_label(substitute(shift_fun)))

  .effect_curve_internal(
    x = x, treat = treat, vcov = vcov,
    range = NULL, n = n,
    data = data, subset = subset_substitute,
    by = by, wts = wts, cluster = cluster, fwb.args = fwb.args,
    type = type, outcome = outcome, shift_fun = shift_fun, shift = shift, ...
  )
}

#' @exportS3Method shift_curve mira
shift_curve.mira <- function(x, treat, shift = NULL, shift_fun = "additive",
                             vcov = "unconditional", cluster = NULL, type = NULL,
                             outcome = NULL, data = NULL, subset = NULL, by = NULL, wts = NULL,
                             n = 51, ...) {

  arg::arg_supplied(treat)

  .check_misplaced_args(...names(), "shift_curve")

  subset_substitute <- substitute(subset)

  shift_fun <- process_shift_fun(shift_fun,
                                 label = .shift_fun_label(substitute(shift_fun)))

  .effect_curve_internal(
    x = x, treat = treat, vcov = vcov,
    range = NULL, n = n,
    data = data, subset = subset_substitute,
    by = by, wts = wts, cluster = cluster,
    type = type, outcome = outcome, shift_fun = shift_fun, shift = shift, ...
  )
}

# How a user-supplied shift function is described in printed output: its name
# when it was supplied by name, its code when that is short, and "custom"
# otherwise
.shift_fun_label <- function(expr) {
  if (is.symbol(expr)) {
    return(as.character(expr))
  }

  if (rlang::is_call(expr, "function")) {
    label <- deparse1(expr, collapse = " ")

    if (nchar(label) <= 50L) {
      return(label)
    }
  }

  "custom"
}
