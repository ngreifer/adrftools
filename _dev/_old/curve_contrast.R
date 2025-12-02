#' Contrast Multiple ADRFs or AMEFs
#'
#' @description
#' `curve_contrast()` computes the pointwise difference between ADRFs or between AMEFs across levels of a stratifying variable supplied to `by` in the original call to [adrf()] or [amef()].
#'
#' @param x an `adrf_est` or `amef_est` object; the output of a call to [adrf()] or [amef()] with `by` supplied.
#' @param \dots ignored.
#'
#' @return An object of class `"curve_contrast"` with elements:
#' \describe{
#'   \item{est}{Estimated differences between ADRFs or AMEFs. These should be examined using [`summary()`][summary.curve_contrast].}
#'   \item{values}{The treatment values used.}
#'   \item{treat}{The treatment variable name.}
#'   \item{vcov}{The estimated variance-covariance matrix of the estimated differences.}
#'   \item{vcov_type}{Type of variance estimation used.}
#'   \item{contrast}{A character vector denoting the contrasts between groups.}
#'   \item{call}{The function call.}
#' }
#'
#' It will additionally inherit from the class of `x`.
#'
#' @details
#' `curve_contrast()` is a helper function for computing the difference between ADRFs or AMEFs computed in subgroups to perform moderation analysis. When multiple subgroups are specified by `by` in the original call to `adrf()` or `amef()`, all pairwise comparisons are included. Use the `subset` argument in the original function call to restrict comparisons to fewer groups.
#'
#'
#' @seealso
#' * [adrf()] and [amef()] for the average dose-response function (ADRF) and average marginal effect function (AMEF), respectively
#' * [summary.curve_contrast()] and [plot.curve_contrast()] to summarize and plot the difference between curves
#' * [curve_test()] for omnibus tests of the difference between curves (e.g., whether the curves differ)
#' * [point_contrast()] for comparing estimates along effect curves
#'
#' @examples
#'

#' @export
curve_contrast <- function(x, ...) {
  .chk_is(x, c("adrf_est", "amef_est"))

  if (inherits(x, "curve_contrast")) {
    .err("`curve_contrast()` cannot be used on a `curve_contrast` object")
  }

  if (is_null(x$by_grid)) {
    .err(sprintf("`%s()` must have been called with `by` specified to use `curve_contrast()`",
                 switch(class(x)[1L], "adrf_est" = "adrf", "amef_est" = "amef")))
  }

  combos <- combn(seq_row(x$by_grid), 2L, simplify = FALSE)

  comp_grid <- expand.grid(val_id = seq_along(x$values),
                           combo = seq_along(combos),
                           KEEP.OUT.ATTRS = FALSE)

  m <- matrix(0, nrow = nrow(comp_grid), ncol = length(x$est))

  by_id <- rep(seq_row(x$by_grid), each = length(x$values))

  for (co in seq_along(combos)) {
    m[comp_grid$combo == co, by_id == combos[[co]][2L]] <- diag(length(x$values))
    m[comp_grid$combo == co, by_id == combos[[co]][1L]] <- -diag(length(x$values))
  }

  est <- tcrossprod(x$est, m)

  vcov <- {
    if (is_null(x$vcov)) NULL
    else quad_mult(m, x$vcov)
  }

  contrast <- vapply(combos, function(combo) {
    sprintf("[%s] - [%s]",
            do.call(function(...) paste(..., sep = ", "),
                    lapply(names(x$by_grid), function(i) {
                      sprintf("%s = %s",
                              i,
                              add_quotes(x$by_grid[[i]][combo[2L]],
                                         chk::vld_character_or_factor(x$by_grid[[i]])))
                    })),
            do.call(function(...) paste(..., sep = ", "),
                    lapply(names(x$by_grid), function(i) {
                      sprintf("%s = %s",
                              i,
                              add_quotes(x$by_grid[[i]][combo[1L]],
                                         chk::vld_character_or_factor(x$by_grid[[i]])))
                    }))
    )
  }, character(1L))

  x[["est"]] <- unname(drop(est))
  x[["vcov"]] <- unname(vcov)
  x[["contrast"]] <- contrast

  x[["by_grid"]] <- NULL
  attr(x, "by_id") <- NULL

  # out <- list(est = unname(drop(est)),
  #             values = x$values,
  #             treat = x$treat,
  #             vcov = unname(vcov),
  #             vcov_type = x$vcov_type,
  #             contrast = contrast,
  #             call = x$call)
  #
  # attr(out, "treat_var") <- .attr(x, "treat_var")

  class(x) <- c("curve_contrast", class(x))

  x
}
