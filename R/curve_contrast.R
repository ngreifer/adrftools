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
  arg::arg_supplied(x)

  check_effect_curve(x, contrast_ok = FALSE, projection_ok = FALSE)

  .by_grid <- .attr(x, ".by_grid")
  .outcome <- .attr(x, ".outcome")

  by_vars <- setdiff(names(.by_grid), .outcome)

  if (is_null(by_vars)) {
    fn_name <- rlang::current_call() |> rlang::call_name()
    arg::err("{.fun {(.curve_fun(x))}} must have been called with {.arg by} specified to use {.fun {fn_name}}")
  }

  # The rows of `.by_grid` for each category of the outcome (or all of them);
  # subgroups are contrasted only within a category
  blocks <- {
    if (is_null(.outcome)) list(seq_row(.by_grid))
    else gsplit(seq_row(.by_grid), .by_grid[[.outcome]], use.g.names = FALSE)
  }

  n_sub <- length(blocks[[1L]])

  if (n_sub < 2L) {
    fn_name <- rlang::current_call() |> rlang::call_name()
    arg::err("{.fun {fn_name}} requires at least two subgroups, but {.arg x} has only one")
  }

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")
  n <- length(.values)

  # All combinations of 'by' levels, within each category of the outcome, as
  # rows of `.by_grid`
  combos0 <- utils::combn(seq_len(n_sub), 2L, simplify = FALSE)

  combos <- lapply(blocks, function(b) {
    lapply(combos0, function(combo) b[combo])
  }) |>
    unlist(recursive = FALSE)

  m <- matrix(0, nrow = n * length(combos), ncol = length(.est))

  for (co in seq_along(combos)) {
    rows <- (co - 1L) * n + seq_len(n)

    m[cbind(rows, (combos[[co]][2L] - 1L) * n + seq_len(n))] <- 1
    m[cbind(rows, (combos[[co]][1L] - 1L) * n + seq_len(n))] <- -1
  }

  by_grid0 <- ss(.by_grid, j = by_vars)

  .contrast <- vapply(combos, function(combo) {
    sprintf("[%s] - [%s]",
            get_by_grid_labels(ss(by_grid0, combo[2L])),
            get_by_grid_labels(ss(by_grid0, combo[1L]))
    )
  }, character(1L))

  # With a model that gives one prediction per category of the outcome, the
  # category of each contrast is kept in `.by_grid`
  .by_grid <- {
    if (is_null(.outcome)) NULL
    else ss(.by_grid, vapply(combos, `[`, integer(1L), 1L), .outcome)
  }

  # Each contrast is weighted by the pooled distribution of the treatment in
  # its two subgroups
  .treat_density <- .attr(x, ".treat_density")

  if (!is_null(.treat_density)) {
    .treat_density[c("y", "tot")] <- list(
      vapply(combos, function(combo) {
        drop(.treat_density[["y"]][, combo, drop = FALSE] %*% .treat_density[["tot"]][combo]) /
          sum(.treat_density[["tot"]][combo])
      }, numeric(length(.treat_density[["x"]]))) |>
        matrix(nrow = length(.treat_density[["x"]])),
      vapply(combos, function(combo) sum(.treat_density[["tot"]][combo]), numeric(1L))
    )
  }

  .make_fun(x,
            .contrast_mat = m,
            .contrast = .contrast,
            .by_grid = .by_grid,
            .treat_density = .treat_density)
}
