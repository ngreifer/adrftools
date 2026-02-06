#' Estimate the average marginal effect function (AMEF)
#'
#' @description
#' `amef()` computes the average marginal effect function (AMEF), the derivative of the average dose-response function (ADRF). This computed from an `adrf_curve` object or from a fitting outcome model directly.
#'
#' @param x an `adrf_curve` object; the output of a call to [adrf()].
#' @param eps numeric; the step size to use when calculating numerical derivatives. Default is `1e-5` (.00001). See Details.
#'
#' @returns
#' An object of class `amef_curve`, which inherits from [`effect_curve`].
#'
#' @details
#' The AMEF is calculated numerically using the central finite derivative formula:
#'
#'  \deqn{\frac{df(x)}\frac{dx} \approx \frac{f(x + e) - f(x - e)}{2e}}
#'
#' The values of the ADRF at the evaluation points are computed using a local polynomial regression as described at [`effect_curve`]. At the boundaries of the ADRF, one-sided derivatives are used.
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [plot.effect_curve()] for plotting the AMEF
#' * [summary.effect_curve()] for testing hypotheses about the AMEF
#' * [`effect_curve`] for computing point estimates along the AMEF
#' * [curve_projection()] for projecting a simpler model onto the AMEF
#' * [reference_curve()] for computing the difference between each point on the AMEF and a specific reference point
#' * [curve_contrast()] for contrasting AMEFs computed within subgroups
#' * [marginaleffects::avg_slopes()] for computing average adjusted slopes for fitted models (similar to the AMEF)
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             Male * (Age + Race + PIR +
#'                       Enough_Food),
#'           data = nhanes3lead)
#'
#' # ADRF of logBLL on Math
#' adrf1 <- adrf(fit, treat = "logBLL")
#'
#' # AMEF of logBLL on Math
#' amef1 <- amef(adrf1)
#'
#' amef1
#'
#' # Plot the AMEF
#' plot(amef1)
#'
#' # AMEF estimates at given points
#' amef1(logBLL = c(0, 1, 2)) |>
#'   summary()

#' @export
amef <- function(x, eps = 1e-5) {
  arg_not_missing(x)

  check_effect_curve(x, amef_ok = FALSE, projection_ok = FALSE)

  .values <- .attr(x, ".values")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")

  eps <- process_eps(eps, .values)

  n_old <- length(.values)

  i0 <- 2 * seq_len(n_old) - 1L
  i1 <- 2 * seq_len(n_old)

  n_new <- 2 * n_old

  values <- numeric(n_new)

  values[i0] <- .values - eps
  values[i1] <- .values + eps

  values[c(1L, 2L)] <- values[c(1L, 2L)] + eps
  values[c(n_new - 1L, n_new)] <- values[c(n_new - 1L, n_new)] - eps

  val_mat <- get_locpoly_w(x = values, v = .values)

  contrast_mat <- matrix(0, nrow = n_old, ncol = n_new)

  diag(contrast_mat[, i0]) <- -.5 / eps
  diag(contrast_mat[, i1]) <- .5 / eps

  n_by <- get_n_by(.contrast, .by_grid)

  contrast_mat <- block_diag(contrast_mat %*% val_mat, n_by)

  .make_fun(x, .contrast_mat = contrast_mat,
            .curve_type = "AMEF")
}
