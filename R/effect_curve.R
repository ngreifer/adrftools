#' Effect curve objects
#' @name effect_curve
#' @aliases effect_curve-class curve_est curve_est-class
#'
#' @description
#' An `effect_curve` object is a function that takes in values of the treatment and produces estimates of the effect curve at those values. `effect_curve` objects are produced by [adrf()] and [shift_curve()] and by functions that modify effect curves, such as [amef()], [curve_contrast()], [reference_curve()], and [curve_projection()]. The output of an `effect_curve` object is a `curve_est` object containing the effect curve estimates. This page describes `effect_curve` and `curve_est` objects.
#'
#' ## Usage
#' \preformatted{f <- adrf(x, ...)
#'
#' f({treat}, subset = NULL, outcome = NULL)}
#'
#' @param \{treat\} the values of the treatment at which to evaluate the effect curve. These must be within the range of the grid points. See Details.
#' @param subset an optional logical expression indicating the subset of the subgroups for which to compute estimates. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param outcome an optional character vector naming the categories of the outcome for which to compute estimates. Can only be used when the model gives a separate prediction for each category of the outcome (e.g., multinomial and ordinal models); see the section *Categorical outcomes* in Details of [adrf()].
#' @param x a `curve_est` object; the output of an `effect_curve` object call.
#' @param digits the number of digits to display.
#' @param \dots arguments passed to [print.data.frame()].
#'
#' @returns
#' A call to an `effect_curve` object returns a `curve_est` object, which is a `data.frame` containing a column for the treatment and a column for the effect curve estimates. `curve_est` objects have `print()`, [`summary()`][summary.curve_est], [`coef()`][coef.curve_est], and [`vcov()`][vcov.curve_est] methods.
#'
#' @details
#' An `effect_curve` object contains a set of grid points on which the effect curve is initially evaluated. The effect curve estimates produced by a call to the `effect_curve` object are interpolated using 3rd-degree local polynomial regression with a Gaussian kernel and bandwidth equal to half the distance between grid points, unless they coincide with the grid points; this means the produced estimates are linear combinations of the grid point estimates.
#'
#' The effect curve can only be evaluated within the range of the grid points, which is set by the `range` argument to [adrf()]; local polynomial regression cannot reliably extrapolate beyond it. Values outside the range of the grid points by no more than 0.1% of its width (e.g., an end point supplied as it was printed, after rounding) are treated as the nearest end point, and values farther outside give an error. To estimate the effect curve at more extreme values of the treatment, supply a wider range to `range` in the call to [adrf()].
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
                      .by_grid, .contrast, .reference, .treat_density, .shift,
                      .outcome) {

  if (missing(.curve_type) && !is_null(fn)) {
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

  if (!is_null(.contrast_mat)) {
    # The contrast matrix has one column per original estimate; contrasts (as
    # in curve_contrast()) can have fewer rows than columns
    est_ind <- seq_len(fncol(.contrast_mat))

    # est is the contrast estimates
    .est <- drop(tcrossprod(.est, .contrast_mat))

    if (.vcov_type == "bootstrap") {
      .boot[["t"]] <- .boot[["t"]] |>
        ss(j = est_ind) |>
        tcrossprod(.contrast_mat)

      .boot[["t0"]] <- .boot[["t0"]][est_ind] |>
        tcrossprod(.contrast_mat) |>
        drop()

      # The map is linear, so the covariance of the transformed draws is the
      # transformed covariance, which is much faster to compute
      .vcov <- quad_mult(.contrast_mat, ss(.vcov, est_ind, est_ind))
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = est_ind) |>
        tcrossprod(.contrast_mat)

      .vcov <- quad_mult(.contrast_mat, ss(.vcov, est_ind, est_ind))
    }
    else if (.vcov_type != "none") {
      .vcov <- quad_mult(.contrast_mat, .vcov)
    }
  }

  fn <- function(x, subset = NULL, outcome = NULL) {
    values <- get(.treat, inherits = FALSE)

    arg::arg_numeric(values, .treat)

    # subset and outcome, which are formals only when they can be used
    fmls <- rlang::fn_fmls_names()

    .subset <- process_subset_by_grid(if (anyv(fmls, "subset")) substitute(subset),
                                      .by_grid = .by_grid,
                                      .contrast = .contrast,
                                      outcome = if (anyv(fmls, "outcome")) outcome,
                                      .outcome = .outcome)

    if (all(.subset)) {
      .s <- seq_along(.est)
    }
    else {
      .s <- which(rep(.subset, each = length(.est) / length(.subset)))

      .by_grid <- ss(.by_grid, .subset)

      if (!is_null(.contrast)) {
        .contrast <- .contrast[.subset]
      }

      .est <- .est[.s]
    }

    n <- length(values)

    # Values just outside the grid are moved to its ends; the estimates keep
    # the values requested as labels
    values_in <- .snap_to_grid(values, .values, .treat, shift = !is_null(.shift))

    val_mat <- matrix(0, nrow = n, ncol = length(.values))

    matches <- max.col(-abs(outer(values_in, .values, "-")), "first")

    to_interp <- abs(values_in - .values[matches]) > .grid_tol(.values)

    if (any(to_interp)) {
      val_mat[to_interp, ] <- get_locpoly_w(x = values_in[to_interp],
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

      # The map is linear, so the covariance of the transformed draws is the
      # transformed covariance
      attr(res, ".vcov") <- quad_mult(val_mat, ss(.vcov, .s, .s))
      attr(res, ".boot") <- .boot
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = .s) |>
        tcrossprod(val_mat)

      attr(res, ".vcov") <- quad_mult(val_mat, ss(.vcov, .s, .s))
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
    attr(res, ".outcome") <- .outcome
    attr(res, ".family") <- .family
    attr(res, ".df") <- .df
    attr(res, ".reference") <- .reference
    attr(res, ".shift") <- .shift

    class(res) <- unique(c("curve_est", class(res)))

    res
  }

  rlang::fn_fmls_names(fn)[1L] <- .treat

  # `subset` refers to the variables in `by`, and `outcome` to the categories
  # of the outcome
  if (!is_null(.contrast) || is_null(setdiff(names(.by_grid), .outcome))) {
    rlang::fn_fmls(fn)[["subset"]] <- NULL
  }

  if (is_null(.outcome)) {
    rlang::fn_fmls(fn)[["outcome"]] <- NULL
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
  attr(fn, ".outcome") <- .outcome
  attr(fn, ".family") <- .family
  attr(fn, ".df") <- .df
  attr(fn, ".response") <- .response
  attr(fn, ".call") <- .call
  attr(fn, ".treat_density") <- .treat_density
  attr(fn, ".shift") <- .shift

  classes <- c("effect_curve", class(fn))

  if (!is_null(.shift)) {
    classes <- c("shift_curve", classes)
  }

  if (.curve_type == "ADRF") {
    classes <- c("adrf_curve", classes)
  }

  if (.curve_type == "AMEF") {
    classes <- c("amef_curve", classes)
  }

  if (!is_null(.reference)) {
    classes <- c("reference_curve", classes)
  }

  if (!is_null(.contrast)) {
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
