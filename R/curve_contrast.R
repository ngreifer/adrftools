#' Contrast multiple subgroup effect curves
#'
#' @description
#' `curve_contrast()` computes the difference between effect curves across levels of a subgrouping variable supplied to `by` in the original call to [adrf()].
#'
#' @param x an [`effect_curve`] object; the output of a call to [adrf()] with `by` supplied.
#'
#' @returns
#' An object of class `contrast_curve`, which inherits from [`effect_curve`], with additional information about the groups being contrasted.
#'
#' @details
#' `curve_contrast()` creates a new effect curve corresponding to the difference between effect curves in two groups. When multiple subgroups are specified by `by` in the original call to [adrf()], all pairwise comparisons are included. Use the `subset` argument in the original function call to restrict comparisons to fewer groups.
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [reference_curve()] for comparing effect curves to a point along the curve
#' * [plot.effect_curve()] for plotting the effect curve contrasts
#' * [summary.curve_est()] for performing tests of effect curve contrasts at specific points
#' * [summary.effect_curve()] for performing omnibus tests of effect curve contrasts (e.g., whether the contrast curve differs from 0)
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             Male * Smoke_in_Home *
#'             (Age + Race + PIR),
#'           data = nhanes3lead)
#'
#' # ADRFs in Race subgroups, excluding Other
#' adrf_by <- adrf(fit, treat = "logBLL",
#'                 by = ~Race,
#'                 subset = Race != "Other")
#'
#' adrf_by
#'
#' # Contrast subgroup ADRFs
#' adrf_contrast <- curve_contrast(adrf_by)
#'
#' adrf_contrast
#'
#' # Plot contrast ADRFs
#' plot(adrf_contrast, simultaneous = FALSE)
#'
#' # Compute and test difference between subgroup
#' # ADRFs at specific points
#' adrf_contrast(logBLL = c(0, 2)) |>
#'   summary()
#'
#' # Test if ADRF differences are present
#' summary(adrf_contrast)

#' @export
curve_contrast <- function(x) {
  arg_not_missing(x)

  check_effect_curve(x, contrast_ok = FALSE, projection_ok = FALSE)

  .by_grid <- .attr(x, ".by_grid")

  if (is_null(.by_grid)) {
    fn_name <- rlang::current_call() |> rlang::call_name()
    .err("{.fun adrf} must have been called with {.arg by} specified to use {.fun {fn_name}}")
  }

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")

  # All combinations of 'by' levels
  combos <- utils::combn(seq_row(.by_grid), 2L, simplify = FALSE)

  comp_grid <- expand.grid(val_id = seq_along(.values),
                           combo = seq_along(combos),
                           KEEP.OUT.ATTRS = FALSE)

  m <- matrix(0, nrow = nrow(comp_grid), ncol = length(.est))

  by_id <- rep(seq_row(.by_grid), each = length(.values))

  for (co in seq_along(combos)) {
    diag(m[comp_grid$combo == co, by_id == combos[[co]][2L]]) <- 1
    diag(m[comp_grid$combo == co, by_id == combos[[co]][1L]]) <- -1
  }

  .contrast <- vapply(combos, function(combo) {
    sprintf("[%s] - [%s]",
            get_by_grid_labels(ss(.by_grid, combo[2L])),
            get_by_grid_labels(ss(.by_grid, combo[1L]))
    )
  }, character(1L))

  .make_fun(x,
            .contrast_mat = m,
            .contrast = .contrast)
}
