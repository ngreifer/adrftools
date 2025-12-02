#' Plot and Summarize Effect Curves
#' @name summary.adrf_est

#' @description
#' `summary()` computes point estimates and confidence intervals for the ADRF or AMEF. `plot()` plots the ADRF or AMEF and its confidence band.
#'
#' @param object,x an object of class `adrf_est`, `amef_est`, or `curve_contrast`; the result of a call to [adrf()], [amef()], or [curve_contrast()].
#' @param conf_level the desired confidence level. Default is .95 for 95% confidence bands.
#' @param simultaneous `logical`; whether the computed confidence bands and p-values should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known as uniform) bands cover the full line at the desired confidence level, whereas pointwise confidence bands only cover each point at the desired level. Default is `TRUE`. See Details.
#' @param null the null value for the hypothesis tests. `NULL` means no p-values or test statistics will be produced. Default is `NULL` when applied to `adrf_est` objects and 0 otherwise. For `plot()`, a horizontal line will be displayed at this value.
#' @param transform whether to compute intervals and perform tests on the transformed estimates. Allowable options include `TRUE`, `FALSE`, or a function specifying a transformation. Ignored unless `object` is an `adrf_est` object. See Details.
#' @param df the "denominator" degrees of freedom to use for the tests and critical test statistics for confidence bands. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise.
#' @param ci.type string; when bootstrapping is used in the original effect curve, what type of confidence interval is to be computed. Allowable options include `"perc"` for percentile intervals, `"wald"` for Wald intervals, and other options allowed by \pkgfun{fwb}{summary.fwb()}. When `simultaneous = TRUE`, only `"perc"` and `"wald"` are allowed. Default is `"perc"`. Ignored when bootstrapping is not used.
#' @param proj an optional `curve_proj` object, the output of a call to [curve_proj()]. Supplying this adds the projection curve to the effect curve plot. See [plot.curve_proj()] to plot the projection curve alone.
#' @param values an optional vector of values of the treatment for which to display the effect curve estimates. If `NULL`, will produce estimates for all points specified in the original effect curve object.
#' @param ... additional arguments passed to other methods. For `plot()`, additional arguments passed to `summary()`.
#'
#' @returns
#' `summary()` returns an object of class `"summary.adrf_est"` or `"summary.amef_est"`, which is a data frame containing the estimates. `plot()` returns a `ggplot` object.
#'
#' @details
#' `summary()` is the preferred method to extract and compute estimates, standard errors, and confidence intervals for estimates of the effect curve (the ADRF or AMEF). For `adrf_est` objects, `summary()` produces the estimated expected potential outcomes at each point on the ADRF along with its standard error and confidence interval. For `amef_est` objects, `summary()` produces the estimated average derivative effect of the treatment at each point on the AMEF along with its standard error, confidence interval, and test statistic and p-value for the test that the estimate is 0. For `curve_contrast` objects, `summary()` produces the difference between the effect curves at each point on the curve along with its standard error, confidence interval, and test statistic and p-value for the test that the difference between estimates is 0. When `values` is supplied, only estimates that specified values will be computed. If those values are not present in the original effect curve, interpolation will be used.
#'
#' `plot()` displays the effect curve in a plot. The solid line represents the effect curve and the ribbon around it represents the confidence bands. When `null` is not `NULL`, an additional flat line at `null` is displayed. When `proj` is supplied, a dashed line corresponding to the projection is added.
#'
#' ## Transform
#'
#' The usual confidence intervals and tests assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. `transform` transforms the estimates to ones that are unbounded and computes the corresponding distribution of transformed estimates using the delta method. By default, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simple a transformation of the estimated points along the curve already computed. Confidence intervals are computed using the transformed estimates before being back-transformed to ensure they are within the bounds of the outcome. When `transform` is specified, standard errors are not reported.
#'
#' ## Simultaneous confidence bands
#'
#' Simultaneous confidence bands ensure the whole effect curve, not just a given individual point, is contained within the band at the given confidence level. These are wider than pointwise bands to reflect that they are covering multiple estimates, which otherwise would decrease the true coverage rate from that specified. `summary()` uses the "sup-t" simultaneous confidence band, which is the smallest one-parameter band that covers the whole effect curve at the desired rate. It is computed similarly to a typical Wald confidence band but with a different critical test statistic. This statistic comes from a complicated distribution and therefore must be simulated, so it can be subject to Monte Carlo error.
#'
#' @seealso
#' * [adrf()], [amef()], and [curve_contrast()] for estimating effect curves
#' * [curve_test()] for testing omnibus null hypotheses about effect curves, such as whether they are flat
#' * [curve_proj()] for projecting the effect curve onto a simpler linear model
#' * [point_contrast()] for comparing estimates along effect curves
#' * [ggplot2::geom_line()] and [ggplot2::geom_ribbon()] for the plotting functions
#'
#' @examples
#'

#' @exportS3Method summary adrf_est
summary.adrf_est <- function(object, conf_level = 0.95, simultaneous = TRUE, null = NULL,
                             transform = TRUE, df = NULL, ci.type = "perc",
                             values = NULL, ...) {

  vcov_type <- object[["vcov_type"]]

  if (is_null(vcov_type)) {
    conf_level <- 0
    null <- NULL
    transform <- FALSE
    simultaneous <- FALSE
  }

  chk::chk_number(conf_level)
  chk::chk_range(conf_level, c(0, 1))

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  if (do_inference) {
    # Process simultaneous
    if (length(coef(object)) > 1L) {
      chk::chk_flag(simultaneous)
    }
    else {
      simultaneous <- FALSE
    }

    if (vcov_type == "bootstrap") {
      chk::chk_string(ci.type)
      if (simultaneous) {
        ci.type <- match_arg(ci.type, c("perc", "wald"))
      }
    }
    else {
      ci.type <- "wald"
    }

    # Process df
    if (is_null(df)) {
      df <- get_df(.attr(object, "model"))
    }
    else {
      chk::chk_number(df)
      chk::chk_gt(df, 2)
    }
  }

  if (is_not_null(values)) {
    values <- process_values(values, n = NULL, .attr(object, "treat_var"),
                             strict = TRUE)
  }

  family <- .attr(object, "model")$family
  transform_list <- process_transform(transform, family)

  res <- .summary_internal(object, conf_level, null = null,
                           simultaneous = simultaneous,
                           transform_list = transform_list,
                           df = df,
                           ci.type = ci.type,
                           values = values)

  class(res) <- c("summary.adrf_est", class(res))

  res
}

#' @exportS3Method summary amef_est
#' @rdname summary.adrf_est
summary.amef_est <- function(object, conf_level = 0.95, simultaneous = TRUE, null = 0,
                             df = NULL, ci.type = "perc",
                             values = NULL, ...) {

  vcov_type <- object[["vcov_type"]]

  if (is_null(vcov_type)) {
    conf_level <- 0
    null <- NULL
    transform <- FALSE
    simultaneous <- FALSE
  }

  chk::chk_number(conf_level)
  chk::chk_range(conf_level, c(0, 1))

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  if (do_inference) {
    # Process simultaneous
    if (length(coef(object)) > 1L) {
      chk::chk_flag(simultaneous)
    }
    else {
      simultaneous <- FALSE
    }

    if (vcov_type == "bootstrap") {
      chk::chk_string(ci.type)
      if (simultaneous) {
        ci.type <- match_arg(ci.type, c("perc", "wald"))
      }
    }
    else {
      ci.type <- "wald"
    }

    # Process df
    if (is_null(df)) {
      df <- get_df(.attr(object, "model"))
    }
    else {
      chk::chk_number(df)
      chk::chk_gt(df, 2)
    }
  }

  if (is_not_null(values)) {
    values <- process_values(values, n = NULL, .attr(object, "treat_var"),
                             strict = TRUE)
  }

  res <- .summary_internal(object, conf_level, null = null,
                           simultaneous = simultaneous,
                           df = df,
                           ci.type = ci.type,
                           values = values)

  class(res) <- c("summary.amef_est", class(res))

  res
}

#' @exportS3Method summary curve_contrast
#' @rdname summary.adrf_est
summary.curve_contrast <- function(object, conf_level = 0.95, simultaneous = TRUE, null = 0,
                                   df = NULL, ci.type = "perc",
                                   values = NULL, ...) {

  vcov_type <- object[["vcov_type"]]

  if (is_null(vcov_type)) {
    conf_level <- 0
    null <- NULL
    transform <- FALSE
    simultaneous <- FALSE
  }

  chk::chk_number(conf_level)
  chk::chk_range(conf_level, c(0, 1))

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  if (do_inference) {
    # Process simultaneous
    if (length(coef(object)) > 1L) {
      chk::chk_flag(simultaneous)
    }
    else {
      simultaneous <- FALSE
    }

    if (vcov_type == "bootstrap") {
      chk::chk_string(ci.type)
      if (simultaneous) {
        ci.type <- match_arg(ci.type, c("perc", "wald"))
      }
    }
    else {
      ci.type <- "wald"
    }

    # Process df
    if (is_null(df)) {
      df <- get_df(.attr(object, "model"))
    }
    else {
      chk::chk_number(df)
      chk::chk_gt(df, 2)
    }
  }

  if (is_not_null(values)) {
    values <- process_values(values, n = NULL, .attr(object, "treat_var"),
                             strict = TRUE)
  }

  res <- .summary_internal(object, conf_level, null = null,
                           simultaneous = simultaneous,
                           df = df,
                           ci.type = ci.type,
                           values = values)

  class(res) <- c("summary.curve_contrast", "summary.adrf_est",
                  class(res))

  res
}

#' @exportS3Method plot adrf_est
#' @rdname summary.adrf_est
plot.adrf_est <- function(x, conf_level = 0.95, simultaneous = TRUE, null = NULL,
                          transform = TRUE, df = NULL, ci.type = "perc", proj = NULL, ...) {
  if (is_not_null(...get("values"))) {
    .err("`values` cannot be supplied to `plot()`")
  }

  summary.adrf_est(object = x, conf_level = conf_level,
                   simultaneous = simultaneous,
                   null = NULL, transform = transform, df = df,
                   ci.type = ci.type, values = NULL, ...) |>
    plot(null = null, proj = proj, ...)
}

#' @exportS3Method plot amef_est
#' @rdname summary.adrf_est
plot.amef_est <- function(x, conf_level = 0.95, simultaneous = TRUE, null = 0,
                          df = NULL, ci.type = "perc", proj = NULL, ...) {
  if (is_not_null(...get("values"))) {
    .err("`values` cannot be supplied to `plot()`")
  }

  summary.amef_est(object = x, conf_level = conf_level,
                   simultaneous = simultaneous,
                   null = NULL, df = df,
                   ci.type = ci.type, values = NULL, ...) |>
    plot(null = null, proj = proj, ...)
}

#' @exportS3Method plot curve_contrast
#' @rdname summary.adrf_est
plot.curve_contrast <- function(x, conf_level = 0.95, simultaneous = TRUE, null = 0,
                                df = NULL, ci.type = "perc", proj = NULL, ...) {
  if (is_not_null(...get("values"))) {
    .err("`values` cannot be supplied to `plot()`")
  }

  summary.curve_contrast(object = x, conf_level = conf_level,
                         simultaneous = simultaneous,
                         null = NULL, df = df,
                         ci.type = ci.type, values = NULL, ...) |>
    plot(null = null, proj = proj, ...)
}

#' @exportS3Method plot summary.adrf_est
plot.summary.adrf_est <- function(x, null = NULL, proj = NULL, ...) {
  chk_proj(x, proj)

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  .plot_internal(x, ylab = "E[Y]", null = null, proj = proj)
}

#' @exportS3Method plot summary.amef_est
plot.summary.amef_est <- function(x, null = 0, proj = NULL, ...) {
  chk_proj(x, proj)

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  .plot_internal(x, ylab = "dE[Y]/da", null = null, proj = proj)
}

#' @exportS3Method plot summary.curve_contrast
plot.summary.curve_contrast <- function(x, null = 0, proj = NULL, ...) {
  chk_proj(x, proj)

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  .plot_internal(x, ylab = "E[Y] Difference", null = null, proj = proj)
}

#' @exportS3Method print summary.adrf_est
#' @rdname summary.adrf_est
print.summary.adrf_est <- function(x, digits = max(3L, getOption("digits") - 3L), topn = 5, ...) {
  chk::chk_whole_number(digits)
  chk::chk_count(topn)

  .print_summary_internal(x = x, digits = digits, topn = topn,
                          ...)

  invisible(x)
}

#' @exportS3Method print summary.amef_est
print.summary.amef_est <- print.summary.adrf_est
