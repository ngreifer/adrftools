#' Effect curve objects
#' @name effect_curve
#' @aliases effect_curve-class curve_est curve_est-class
#'
#' @description
#' An `effect_curve` object is a function that takes in values of the treatment and produces estimates of the effect curve at those values. `effect_curve` objects are produces by [adrf()] and functions that modify effect curves, such as [amef()], [curve_contrast()], [reference_curve()], and [curve_projection()]. The output of an `effect_curve` object is a `curve_est` object containing the effect curve estimates. This page describes `effect_curve` and `curve_est` objects.
#'
#' ## Usage
#' \preformatted{
#' f <- adrf(x, ...)
#'
#' f({treat}, subset = NULL)
#' }
#'
#' @param \{treat\} the values of the treatment at which to evaluate the effect curve.
#' @param subset an optional logical expression indicating the subset of the subgroups for which to compute estimates. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param x a `curve_est` object; the output of an `effect_curve` object call.
#' @param digits the number of digits to display.
#' @param \dots arguments passed to [print.data.frame()].
#'
#' @returns
#' A call to an `effect_curve` object returns a `curve_est` object, which is a `data.frame` containing a column for the treatment and a column for the effect curve estimates. `curve_est` objects have `print()`, [`summary()`][summary.curve_est], [`ceof()`][coef.curve_est], and [`vcov()`][vcov.curve_est] methods.
#'
#' @details
#' An `effect_curve` object contains a set of grid points on which the effect curve is initially evaluated. The effect curve estimates produced by a call to the `effect_curve` object are interpolated using 3rd-degree local polynomial regression with a Gaussian kernel and bandwidth equal to half the distance between grid points, unless they coincide with the grid points; this means the produced estimates are linear combinations of the grid point estimates.
#'
#' @seealso
#' * [adrf()] for generating an effect curve
#' * [summary.curve_est()] for performing inference on effect curve estimates
#' * [plot.effect_curve()] for plotting the effect curve
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
#' adrf1
#'
#' # Compute estimates along the ADRF
#' adrf1(logBLL = c(0, 1, 2))
#'
#' # Perform inference on the estimates
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary()
#'
#' # ADRF within groups defined by `Male`
#' adrf2 <- adrf(fit, treat = "logBLL",
#'               by = ~Male)
#'
#' adrf2
#'
#' # Estimates in both groups
#' adrf2(logBLL = c(0, 1, 2))
#'
#' # Estimates in one group
#' adrf2(logBLL = c(0, 1, 2), subset = Male == 1)
NULL

.make_fun <- function(fn = NULL, .contrast_mat = NULL,
                      .est, .vcov, .values, .treat, .vcov_type, .boot, .draws,
                      .curve_type, .family, .df, .response, .call,
                      .by_grid, .contrast, .reference) {

  if (missing(.curve_type) && is_not_null(fn)) {
    if (inherits(fn, "adrf_curve")) {
      .curve_type <- "ADRF"
    }
    else if (inherits(fn, "amef_curve")) {
      .curve_type <- "AMEF"
    }
  }

  for (i in rlang::fn_fmls_names()[-(1:2)]) {
    if (rlang::is_missing(environment()[[i]])) {
      if (is_null(fn)) {
        assign(i, NULL)
      }
      else {
        assign(i, .attr(fn, i))
      }
    }
  }

  if (is_not_null(.contrast_mat)) {
    # est is the contrast estimates
    .est <- drop(tcrossprod(.est, .contrast_mat))

    if (.vcov_type == "bootstrap") {
      .boot[["t"]] <- .boot[["t"]] |>
        ss(j = seq_along(.est)) |>
        tcrossprod(.contrast_mat)

      .boot[["t0"]] <- .boot[["t0"]][seq_along(.est)] |>
        tcrossprod(.contrast_mat) |>
        drop()

      .vcov <- stats::vcov(.boot)
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = seq_along(.est)) |>
        tcrossprod(.contrast_mat)

      .vcov <- stats::cov(.draws)
    }
    else if (.vcov_type != "none") {
      .vcov <- quad_mult(.contrast_mat, .vcov)
    }
  }

  fn <- function(x, subset = NULL) {
    values <- get(.treat, inherits = FALSE)

    chk::chk_numeric(values, add_quotes(.treat, "`"))

    # subset
    if (any(rlang::fn_fmls_names() == "subset") && is_null(.contrast)) {
      .subset <- process_subset_by_grid(substitute(subset),
                                        .by_grid = .by_grid,
                                        .contrast = .contrast)

      .s <- which(rep(.subset, each = length(.est) / nrow(.by_grid)))

      .by_grid <- ss(.by_grid, .subset)

      .est <- .est[.s]
    }
    else {
      .s <- seq_along(.est)
    }

    n <- length(values)

    val_mat <- matrix(0, nrow = n, ncol = length(.values))

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
    .est0 <- drop(val_mat %*% .est)

    res_names <- c(.treat, "estimate")
    res <- make_df(res_names, length(.est0))

    res$estimate <- .est0

    res <- add_est_labels(res, .contrast, .by_grid, values, .treat)

    if (.vcov_type == "bootstrap") {
      .boot[["t"]] <- .boot[["t"]] |>
        ss(j = .s) |>
        tcrossprod(val_mat)

      .boot[["t0"]] <- .boot[["t0"]][.s] |>
        tcrossprod(val_mat) |>
        drop()

      attr(res, ".vcov") <- stats::vcov(.boot)
      attr(res, ".boot") <- .boot
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = .s) |>
        tcrossprod(val_mat)

      attr(res, ".vcov") <- stats::cov(.draws)
      attr(res, ".draws") <- .draws
    }
    else if (.vcov_type != "none") {
      attr(res, ".vcov") <- quad_mult(val_mat, ss(.vcov, .s, .s))
    }

    attr(res, "values") <- values

    attr(res, ".treat") <- .treat
    attr(res, ".vcov_type") <- .vcov_type
    attr(res, ".curve_type") <- .curve_type
    attr(res, ".by_grid") <- .by_grid
    attr(res, ".contrast") <- .contrast
    attr(res, ".family") <- .family
    attr(res, ".df") <- .df
    attr(res, ".reference") <- .reference

    class(res) <- unique(c("curve_est", class(res)))

    res
  }

  rlang::fn_fmls_names(fn)[1L] <- .treat

  if (is_null(.by_grid) || is_not_null(.contrast)) {
    rlang::fn_fmls(fn)[-seq_along(.treat)] <- NULL
  }

  attr(fn, ".est") <- .est
  attr(fn, ".vcov") <- .vcov
  attr(fn, ".values") <- .values
  attr(fn, ".treat") <- .treat
  attr(fn, ".vcov_type") <- .vcov_type
  attr(fn, ".boot") <- .boot
  attr(fn, ".draws") <- .draws
  attr(fn, ".reference") <- .reference
  attr(fn, ".by_grid") <- .by_grid
  attr(fn, ".contrast") <- .contrast
  attr(fn, ".family") <- .family
  attr(fn, ".df") <- .df
  attr(fn, ".family") <- .family
  attr(fn, ".response") <- .response
  attr(fn, ".call") <- .call

  classes <- c("effect_curve", class(fn))

  if (.curve_type == "ADRF") {
    classes <- c("adrf_curve", classes)
  }

  if (.curve_type == "AMEF") {
    classes <- c("amef_curve", classes)
  }

  if (is_not_null(.reference)) {
    classes <- c("reference_curve", classes)
  }

  if (is_not_null(.contrast)) {
    classes <- c("contrast_curve", classes)
  }

  class(fn) <- unique(classes)

  fn
}

#' @exportS3Method print curve_est
print.curve_est <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  summary(x, conf_level = 0, null = NA) |>
    print(digits = digits, ...)
}
