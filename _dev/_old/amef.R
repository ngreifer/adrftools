#' Average Marginal Effect Function Estimation
#'
#' @description
#' Estimates the average marginal effect function (AMEF) for a fitted model object,
#' evaluating average derivative effects over a range of treatment values. The AMEF is the derivative of the average dose-response function (ADRF).
#'
#' @inheritParams adrf
#' @param x a fitted model object (e.g., from `lm()`, `glm()`, etc.). For `print()`, an `amef_est` object.
#' @param values an optional numeric vector of treatment values at which to compute the AMEF.
#'   If `NULL`, a sequence of `n` values spanning the treatment range is used.
#' @param n integer specifying the number of equally spaced values for which to compute the AMEF when `values` is `NULL`.
#' @param wts optional numeric vector of weights to generalize the AMEF to a weighted population.
#' @param eps the constant used to perform numeric differentiation. Default is `2e-5`, which corresponds to .00002 times the standard deviation of the treatment variable.
#'
#' @returns
#' An object of class `"amef_est"` with elements:
#' \describe{
#'   \item{est}{Estimated average derivative effect values. These should be examined using [summary.amef_est()].}
#'   \item{values}{The treatment values used.}
#'   \item{treat}{The treatment variable name.}
#'   \item{vcov}{The estimated variance-covariance matrix of the estimated values.}
#'   \item{vcov_type}{Type of variance estimation used.}
#'   \item{by_grid}{Grid of combinations for `by` variables, if used.}
#'   \item{call}{The function call.}
#' }
#'
#' @details
#' `amef()` estimates the AMEF by numerically differentiating average
#' predicted outcomes in the sample for counterfactual treatment values, optionally stratified
#' by grouping variables and accounting for estimation uncertainty via
#' unconditional or conditional variance estimation or bootstrapping. Unconditional variance estimation and bootstrapping treat the sample as random. Unconditional variance estimation requires [sandwich::estfun()] and [sandwich::bread()] methods for the supplied object to be available.
#'
#' The AMEF is computed similarly to the ADRF: for each supplied value of the treatment, the treatment variable for all units is set to that value, and the average of the predictions from the model (weighted when `wts` is specified) is used as the estimated expected potential outcome under that treatment value. The difference between the average prediction adding `eps/2` to each treatment value and the average prediction subtracting `eps/2` from each treatment value divided by `eps` is used to numerically approximate the AMEF. The variance is computed using the formula in Hansen et al. (2024) when `vcov = "unconditional"`, which involves augmenting the influence function with a term to account for sampling from the superpopulation, or the delta method based on the original parameter estimates and their covariance otherwise.
#'
#' ## Specific methods
#'
#' * for `lm` objects, `type` is not used; predictions are automatically on the response scale.
#' * for `svyglm` for \pkg{survey}, `wts` is automatically set to `TRUE` and cannot be changed.
#' * for `mira` objects from \pkg{mice} and `mimira` objects from \pkg{MatchThem}, analyses are applied to each imputed dataset and pooled using Rubin's rules. Bootstrapping is not allowed with with objects.
#'
#' @seealso
#' * [adrf()] for the average dose-response function (ADRF)
#' * [summary.amef_est()] and [plot.amef_est()] to summarize and plot the AMEF
#' * [curve_contrast()] to compare multiple AMEFs when `by` is specified
#' * [curve_test()] for omnibus tests of the AMEF
#' * [point_contrast()] for comparing estimates along an AMEF
#' * [marginaleffects::avg_slopes()] for a more general implementation of computing average slopes
#'
#' @examples
#'

#' @export
amef <- function(x, ...) {
  UseMethod("amef")
}

#' @exportS3Method amef lm
#' @rdname amef
amef.lm <- function(x, treat, vcov = "unconditional", values = NULL, n = 51,
                     subset = NULL, by = NULL, wts = NULL, eps = 2e-5, ...) {

  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x)

  grad_fun <- .get_grad_fun(x)

  pred_fun <- .get_pred_fun(x)

  .est_internal(
    .est_type = "amef", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts, eps = eps,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method amef glm
amef.glm <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                      subset = NULL, by = NULL, wts = NULL, eps = 2e-5, ...) {
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
    .est_type = "amef", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts, eps = eps,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method amef default
#' @rdname amef
amef.default <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                          subset = NULL, by = NULL, wts = NULL, eps = 2e-5, ...) {
  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  chk::chk_string(type)

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x)

  grad_fun <- .get_grad_fun(x)

  pred_fun <- .get_pred_fun(x, type = type)

  .est_internal(
    .est_type = "amef", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts, eps = eps,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method amef mira
amef.mira <- function(x, treat, vcov = "unconditional", values = NULL, n = 51, type = "response",
                      subset = NULL, by = NULL, wts = NULL, eps = 2e-5, ...) {
  mcall <- process_call()

  chk::chk_not_missing(treat, "`treat`")

  chk::chk_string(type)

  subset_substitute <- substitute(subset)

  estimator <- .get_estimator(x[["analyses"]][[1L]])

  grad_fun <- .get_grad_fun(x[["analyses"]][[1L]])

  pred_fun <- .get_pred_fun(x[["analyses"]][[1L]], type = type)

  .est_internal_mi(
    .est_type = "amef", x = x, treat = treat, vcov = vcov,
    values = values, n = n,
    subset_substitute = subset_substitute, by = by, wts = wts, eps = eps,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    mcall = mcall,
    ...
  )
}

#' @exportS3Method print amef_est
#' @rdname amef
print.amef_est <- print.adrf_est
