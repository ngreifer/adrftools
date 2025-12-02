#' Average Dose Response Function Estimation
#'
#' @description
#' Estimates the average dose response function (ADRF) for a fitted model object,
#' evaluating expected potential outcomes over a range of treatment values.
#'
#' @param x a fitted model object (e.g., from `lm()`, `glm()`, etc.). For `print()`, an `adrf_est` object.
#' @param treat a string specifying the name of the treatment variable.
#' @param vcov how the covariance matrix of the estimates should be computed. If `"unconditional"` (default), use the sandwich estimator including sampling uncertainty. If `"boot"` or `"fwb"`, use the traditional or fractional weighted bootstrap, respectively (both of which require the \CRANpkg{fwb} package to be installed). Otherwise, may be a covariance matrix or other allowed input to the `vcov` argument of [marginaleffects::get_vcov()]. Can also be `"none"` to avoid computing the uncertainty.
#' @param values an optional numeric vector of treatment values at which to compute the ADRF.
#'   If `NULL`, a sequence of `n` values spanning the treatment range is used.
#' @param n integer specifying the number of equally spaced values for which to compute the ADRF when `values` is `NULL`.
#' @param type character string indicating the type of prediction. Passed to [marginaleffects::get_predict()]. Default is `"response"` for predictions on the scale of the outcome variable. Other options might include `"link"` for the linear predictor.
#' @param subset an optional logical expression indicating the subset of data to use for estimation. Will be evaluated in the environment of the original dataset supplied to the model fitting function.
#' @param by optional variable(s) over which to group the estimation. Can be a character vector or one-sided formula.
#' @param wts optional numeric vector of weights to generalize the ADRF to a weighted population.
#' @param digits the number of digits to print.
#' @param topn half the number of estimates to print.
#' @param ... further arguments passed to [marginaleffects::get_predict()] or \pkgfun{fwb}{fwb}.
#'
#' @returns
#' An object of class `"adrf_est"` with elements:
#' \describe{
#'   \item{est}{Estimated average treatment response values. These should be examined using [summary.adrf_est()].}
#'   \item{values}{The treatment values used.}
#'   \item{treat}{The treatment variable name.}
#'   \item{vcov}{The estimated variance-covariance matrix of the estimated values.}
#'   \item{vcov_type}{Type of variance estimation used.}
#'   \item{by_grid}{Grid of combinations for `by` variables, if used.}
#'   \item{call}{The function call.}
#' }
#'
#' @details
#' `ardf()` estimates the ADRF by computing average
#' predicted outcomes in the sample for counterfactual treatment values, optionally stratified
#' by grouping variables and accounting for estimation uncertainty via
#' unconditional or conditional variance estimation or bootstrapping. Unconditional variance estimation and bootstrapping treat the sample as random. Unconditional variance estimation requires [sandwich::estfun()] and [sandwich::bread()] methods for the supplied object to be available.
#'
#' For each value of the treatment supplied to `values` (or generated automatically by `n`), the treatment variable for all units is set to that value, and the average of the predictions from the model (weighted when `wts` is specified) is used as the estimated expected potential outcome under that treatment value. The variance is computed using the formula in Hansen et al. (2024) when `vcov = "unconditional"`, which involves augmenting the influence function with a term to account for sampling from the superpopulation, or the delta method based on the original parameter estimates and their covariance otherwise.
#'
#' ## Specific methods
#'
#' * for `lm` objects, `type` is not used; predictions are automatically on the response scale.
#' * for `svyglm` objects from \pkg{survey}, `wts` is automatically set to `TRUE` and cannot be changed.
#' * for `mira` objects from \pkg{mice} and `mimira` objects from \pkg{MatchThem}, analyses are applied to each imputed dataset and pooled using Rubin's rules. Bootstrapping is not allowed with with objects.
#'
#' @seealso
#' * [amef()] for the average marginal effect function (AMEF), the derivative of the ADRF
#' * [summary.adrf_est()] and [plot.adrf_est()] to summarize and plot the ADRF
#' * [curve_contrast()] to compare multiple ADRFs when `by` is specified
#' * [curve_test()] for omnibus tests of the ADRF
#' * [point_contrast()] for comparing estimates along an ADRF
#' * [marginaleffects::avg_predictions()] for a more general implementation of computing average adjusted predictions
#'
#' @examples
#'

#' @export
adrf <- function(x, ...) {
  UseMethod("adrf")
}

#' @exportS3Method adrf lm
#' @rdname adrf
adrf.lm <- function(x, treat, vcov = "unconditional", values = NULL, n = 51,
                    subset = NULL, by = NULL, wts = NULL, ...) {

  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x)

  grad_fun <- .get_grad_fun(x)

  pred_fun <- .get_pred_fun(x)

  .est_internal(
    .est_type = "adrf", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method adrf glm
adrf.glm <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                     subset = NULL, by = NULL, wts = NULL, ...) {
  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  chk::chk_string(type)

  subset_substitute <- substitute(subset)

  if (type == "link" ||
      (type == "response" && identical(x$family[["link"]], "identity"))) {
    estimator <- .get_estimator.lm(x)

    grad_fun <- .get_grad_fun.lm(x)

    pred_fun <- .get_pred_fun.lm(x)
  }
  else if (type == "response") {
    estimator <- .get_estimator(x)

    grad_fun <- .get_grad_fun(x)

    pred_fun <- .get_pred_fun(x, type = type)
  }
  else {
    estimator <- .get_estimator.default(x)

    grad_fun <- .get_grad_fun.default(x)

    pred_fun <- .get_pred_fun.default(x, type = type)
  }

  .est_internal(
    .est_type = "adrf", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method adrf default
#' @rdname adrf
adrf.default <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                         subset = NULL, by = NULL, wts = NULL, ...) {
  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  chk::chk_string(type)

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x)

  grad_fun <- .get_grad_fun(x)

  pred_fun <- .get_pred_fun(x, type = type)

  .est_internal(
    .est_type = "adrf", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method adrf mira
adrf.mira <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                      subset = NULL, by = NULL, wts = NULL, ...) {
  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  chk::chk_string(type)

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x[["analyses"]][[1L]])

  grad_fun <- .get_grad_fun(x[["analyses"]][[1L]])

  pred_fun <- .get_pred_fun(x[["analyses"]][[1L]], type = type)

  .est_internal_mi(
    .est_type = "adrf", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method print adrf_est
#' @rdname adrf
.print.adrf_est <- function(x, digits = max(3L, getOption("digits") - 3L), topn = 5, ...) {
  chk::chk_whole_number(digits)
  chk::chk_count(topn)

  .print_internal(x = x, digits = digits, topn = topn, ...)
}

#' @exportS3Method print adrf_est
print.adrf_est <- function(x, ...) {
  cl <- class(x)[1L]
  cat(sprintf("%s %s object\n\n",
              switch(cl, amef_est = , adrf_est = "An", "A"),
              .it(cl)))

  cat(sprintf(" - treatment: %s\n", x[["treat"]]))
  cat(sprintf("   + range: %s to %s\n",
              round(min(x$values), 4),
              round(max(x$values), 4)))
  cat(sprintf("   + evaluated at %s points\n",
              length(x$values)))

  if (is_not_null(x[["by_grid"]])) {
    cat(sprintf(" - by: %s\n", toString(names(x[["by_grid"]]))))
  }

  if (is_not_null(x[["contrast"]])) {
    cat(sprintf(" - %s: %s\n",
                ngettext(length(x[["contrast"]]), "contrast", "contrasts"),
                toString(x[["contrast"]])))
  }

  cat(sprintf(" - inference: %s\n",
              if (is_null(x[["vcov_type"]])) "none" else x[["vcov_type"]]))

  cat(.it("\nUse `plot()` to plot the curve.\n"))
}
