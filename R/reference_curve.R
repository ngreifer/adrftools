#' Compare an Effect Curve to a Single Point
#'
#' @description
#' `reference_curve()` creates a new effect curve as the contrast between each point on a given effect curve and a specified point along that curve. The new curve is called a "reference effect curve".
#'
#' @param x an [`effect_curve`] object; the output of a call to [adrf()] or a function that modifies it.
#' @param reference numeric; the value of the treatment to use as the reference value.
#'
#' @returns
#' An object of class `reference_curve`, which inherits from [`effect_curve`], with the value supplied to `reference` as an additional attribute.
#'
#' @details
#' The value supplied to `reference` is added as a grid point on the reference effect curve using the interpolation method described in [`effect_curve`]. The delta method is used to compute the variance of the difference between each point along the effect curve and the reference point.
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [plot.effect_curve()] for plotting the reference effect curve
#' * [summary.effect_curve()] for testing hypotheses about the reference effect curve
#' * [summary.curve_est()] for performing inference on individual points on an effect curve, including a reference effect curve
#' * [point_contrast()] for effect curve estimates to each other (rather than to a single point)
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             (Male + Age + Race + PIR +
#'                Enough_Food),
#'           data = nhanes3lead)
#'
#' # ADRF of logBLL on Math, unconditional
#' # inference
#' adrf1 <- adrf(fit, treat = "logBLL")
#'
#' # Differences between ADRF estimates and estimate
#' # at `logBLL = 0`
#' ref1 <- reference_curve(adrf1, reference = 0)
#'
#' ref1
#'
#' # Plot the reference effect curve
#' plot(ref1)
#'
#' # Reference effect curve estimates at given points
#' ref1(logBLL = c(0, 1, 2)) |>
#'   summary()
#'
#' # Test if reference effect curve is 0 (equivalent
#' # to testing if ADRF is flat)
#' summary(ref1)

#' @export
reference_curve <- function(x, reference) {
  chk::chk_not_missing(x, "`x`")
  check_effect_curve(x, reference_ok = FALSE, projection_ok = FALSE)

  chk::chk_not_missing(reference, "`reference`")
  check_reference(reference, values = .attr(x, ".values"))

  x <- .add_values(x, reference)

  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")
  .values <- .attr(x, ".values")

  n <- length(.values)
  ref_ind <- match(reference, .values)

  # contrast_mat0 encodes contrasts between estimates at values
  contrast_mat0 <- diag(1, n)
  contrast_mat0[, ref_ind] <- contrast_mat0[, ref_ind] - 1

  n_by <- get_n_by(.contrast, .by_grid)

  contrast_mat <- block_diag(contrast_mat0, n_by)

  .make_fun(x,
            .contrast_mat = contrast_mat,
            .reference = reference)
}

.add_values <- function(x, new_values) {
  .values <- .attr(x, ".values")

  if (all(new_values %in% .values)) {
    return(x)
  }

  n_old <- length(.values)

  .est <- .attr(x, ".est")
  .vcov_type <- .attr(x, ".vcov_type")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")

  values <- funique(c(.values, new_values), sort = TRUE)

  n_new <- length(values)

  val_mat <- matrix(0, nrow = n_new, ncol = n_old)

  matches <- max.col(-abs(outer(values, .values, "-")), "first")

  to_interp <- !check_if_zero(values - .values[matches])

  if (any(to_interp)) {
    val_mat[to_interp, ] <- get_locpoly_w(x = values[to_interp],
                                          v = .values)
  }

  if (!all(to_interp)) {
    val_mat[cbind(which(!to_interp), matches[!to_interp])] <- 1
  }

  n_by <- get_n_by(.contrast, .by_grid)

  val_mat <- block_diag(val_mat, n_by)

  # est0 is the estimates at values, including interp'd
  new_est <- drop(val_mat %*% .est)

  if (.vcov_type == "none") {
    return(.make_fun(x,
                     .est = new_est,
                     .values = values))
  }

  if (.vcov_type == "bootstrap") {
    new_boot <- .attr(x, ".boot")

    new_boot[["t"]] <- new_boot[["t"]] |>
      ss(j = seq_along(.est)) |>
      tcrossprod(val_mat)

    new_boot[["t0"]] <- new_boot[["t0"]][seq_along(.est)] |>
      tcrossprod(val_mat) |>
      drop()

    return(.make_fun(x,
                     .est = new_est,
                     .values = values,
                     .vcov = stats::vcov(new_boot),
                     .boot = new_boot))
  }

  if (.vcov_type == "posterior") {
    new_draws <- .attr(x, ".draws") |>
      ss(j = seq_along(.est)) |>
      tcrossprod(val_mat)

    return(.make_fun(x,
                     .est = new_est,
                     .values = values,
                     .vcov = stats::cov(new_draws),
                     .draws = new_draws))
  }

  .vcov <- .attr(x, ".vcov")

  .make_fun(x,
            .est = new_est,
            .values = values,
            .vcov = quad_mult(val_mat, .vcov))
}
