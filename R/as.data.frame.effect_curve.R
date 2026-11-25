#' Extract estimates and confidence bands of an effect curve
#'
#' @description
#' `as.data.frame()` returns the estimates of an effect curve at its grid points and the limits of their confidence bands as a data frame. These are the values `plot()` displays, so the data frame can be used to build a custom plot of the effect curve.
#'
#' @inheritParams plot.effect_curve
#' @param row.names `NULL` or a character vector of row names for the returned data frame. Default is `NULL` to use the row numbers.
#' @param optional ignored.
#'
#' @returns
#' A data frame with one row for each grid point of the effect curve within each subgroup or contrast, ordered by subgroup or contrast and then by treatment value. It has the following columns:
#' * the variables supplied to `by` in the call to [adrf()] (when supplied), with one value for each subgroup, or `contrast`, the label of each contrast, when `x` is the output of a call to [curve_contrast()]. For models that give a separate prediction for each category of the outcome, a column (usually named `outcome`) contains the category.
#' * a column named for the treatment (or `shift` for a shift curve) containing the grid points.
#' * `estimate`: the estimate of the effect curve at each grid point.
#' * `conf.low` and `conf.high`: the lower and upper limits of the confidence band at each grid point. These are omitted when `vcov = "none"` in the call to [adrf()].
#' * `proj_estimate`: the estimate of the projection at each grid point, when `proj` is supplied.
#'
#' @details
#' The confidence bands are computed exactly as they are by [plot.effect_curve()], with the same arguments, so a call to `as.data.frame()` and a call to `plot()` with the same arguments give the same band (with Wald inference and `simultaneous = TRUE`, the critical value is computed using a randomized algorithm, so [set.seed()] should be called before each for the bands to be identical). To compute estimates and confidence intervals at specific treatment values rather than at the grid points, call the effect curve with those values and use [summary.curve_est()] on the result.
#'
#' @seealso
#' * [plot.effect_curve()] for plotting an effect curve
#' * [summary.curve_est()] for confidence intervals at specific treatment values
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 3) * Male + Age + Race +
#'             PIR + NICU + Smoke_Pregnant,
#'           data = nhanes3lead)
#'
#' # ADRFs of logBLL on Math within subgroups
#' adrf1 <- adrf(fit, treat = "logBLL", by = ~Male,
#'               n = 21)
#'
#' # Estimates and simultaneous 95% confidence bands
#' # at the grid points
#' set.seed(123)
#' bands <- as.data.frame(adrf1)
#'
#' head(bands)
#'
#' # Pointwise bands for one subgroup, with the linear
#' # projection of the ADRF
#' as.data.frame(adrf1, simultaneous = FALSE,
#'               subset = Male == 1,
#'               proj = curve_projection(adrf1, "linear")) |>
#'   head()
#'
#' @exportS3Method as.data.frame effect_curve
as.data.frame.effect_curve <- function(x, row.names = NULL, optional = FALSE,
                                       conf_level = 0.95, simultaneous = TRUE,
                                       subset = NULL, outcome = NULL, proj = NULL,
                                       transform = TRUE, ci.type = "perc",
                                       df = NULL, ...) {
  check_proj(x, proj)

  .values <- .attr(x, ".values")
  .treat <- .attr(x, ".treat")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")

  band <- .curve_band_data(x, conf_level = conf_level, simultaneous = simultaneous,
                           subset = substitute(subset), outcome = outcome,
                           proj = proj, transform = transform, ci.type = ci.type,
                           df = df)

  res <- band[["res"]]
  .subset <- band[["subset"]]

  # One row per treatment value in each selected subgroup or contrast, with the
  # variables defining the subgroups (and the category of the outcome) or the
  # contrast first
  rows <- rep(which(.subset), each = length(.values))

  id <- NULL

  if (!is_null(.contrast)) {
    id <- data.frame(contrast = .contrast[rows])

    # The category of the outcome of each contrast
    if (!is_null(.by_grid)) {
      id <- cbind(id, ss(.by_grid, rows))
    }
  }
  else if (!is_null(.by_grid)) {
    id <- ss(.by_grid, rows)
  }

  out <- {
    if (is_null(id)) res
    else cbind(id, res)
  }

  rownames(out) <- row.names

  out
}
