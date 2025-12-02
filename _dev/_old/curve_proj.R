#' Effect Curve Projection

#' @description
#' `curve_proj()` produces a projection of an estimated ADRF or AMEF onto a specified linear model that is a function only of the treatment to act as a more interpretable summary of the original effect curve.
#'
#' @inheritParams summary.adrf_est
#' @param x for `curve_proj()`an object of class `adrf_est`, `amef_est`, or `curve_contrast`; the result of a call to [adrf()], [amef()], or [curve_contrast()]. For `plot()`, see `object` below.
#' @param formula a one-sided formula, the terms of which must only contain the treatment variable.
#' @param object a `curve_proj` object; the output of a call to `curve_proj()`.
#' @param \dots for `plot()`, further arguments passed to [plot.adrf_est()] or [plot.amef_est()]. Ignored otherwise.
#'
#' @returns
#' `curve_proj()` returns a `curve_proj` object, which contains components similar to that of a call to [lm()], including `coefficients`, `residuals`, `fitted.values`, and `terms`. These refer to the projection model. The `vcov` component contains the covariance matrix of the projection parameters. The coefficients and covariance matrix should be extracted with `coef()` and `vcov()`, respectively.
#'
#' @details
#' The projection model can be thought of as a linear regression of the effect curve estimates on the treatment. Whereas the original effect curve may be complicated and nonlinear, the projection model can be simple and easily interpretable, though it must be understood as a summary of the original effect curve. For example, the original ADRF might have been computed from an outcome model that involves treatment splines, covariates, and treatment-covariate interactions. Though the ADRF is a univariable function (i.e., of only the treatment), it isn't described by a single set of parameters. The linear projection of the ADRF, though, could be a simple linear model, described by an intercept and the slope on treatment. Though only a rough approximation to the ADRF, the linear projection may be more easily interpreted. This concept is described in Neugebauer and van der Laan (2007).
#'
#' `curve_proj()` fits this projection model and accounts for the uncertainty in the estimates of the effect curve in computing uncertainty estimates for the projection model parameters. Because the true effect curve is continuous, the model is fit minimizing \deqn{\int_{a_\text{lower}}^{a_\text{upper}}{\left(\hat{\theta}(a)-\hat{\mathbf{\beta}} B(a) \right)^2 da}} where \eqn{\hat{\theta}(a)} is the effect curve estimate at treatment value \eqn{a}, \eqn{B(a)} is the basis function representation of \eqn{a} (i.e., as specified in `formula`), and \eqn{\hat{\mathbf{\beta}}} is the vector of projection parameters to be estimated. This integral is approximated using a trapezoidal Riemann sum over the effect curve points evaluated in `x`. Because this projection is a linear operator over the effect curve, the covariance of the projection parameters can be computed using the delta method applied to the estimated covariance of the original effect curve estimates.
#'
#' `summary()` produces a summary table of the projection coefficient estimates and associated uncertainty estimates. `plot()` plots the projection curve as if it was an effect curve; see [plot.adrf_est()] for adding the projection curve to the plot of the original effect curve.
#'
#' @seealso
#' * [adrf()], [amef()], and [curve_contrast()] for estimating effect curves
#' * [curve_test()] for testing omnibus null hypotheses about effect curves, such as whether they are flat
#' * [ggplot2::geom_line()] and [ggplot2::geom_ribbon()] for the plotting functions
#'
#' @references
#' Neugebauer, R., & van der Laan, M. (2007). Nonparametric causal effects based on marginal structural models. *Journal of Statistical Planning and Inference*, 137(2), 419–434. \doi{/10.1016/j.jspi.2005.12.008}
#'
#' @examples
#'

#' @export
curve_proj <- function(x, formula, transform = TRUE) {
  mcall <- process_call()

  chk::chk_not_missing(x, "`x`")

  .chk_is(x, c("adrf_est", "amef_est"))

  chk::chk_not_missing(formula, "`formula`")

  if (!rlang::is_formula(formula, lhs = FALSE)) {
    .err("`formula` must be a one-sided formula with the projection model on the right-hand side")
  }

  empirical <- FALSE
  chk::chk_flag(empirical)

  n <- length(x[["values"]])

  by_list <- {
    if (is_null(x$by_grid) && is_null(x$contrast)) list(seq_len(n))
    else if (is_not_null(x$contrast)) split(seq_along(x$est), rep(seq_along(x$contrast), each = n))
    else split(seq_along(x$est), rep(seq_row(x$by_grid), each = n))
  }

  #Check that only treat is in formula
  treat <- x[["treat"]]
  vars_in_formula <- get_varnames(formula)

  if (!all(get_varnames(formula) %in% treat)) {
    .err(sprintf("only the treatment variable `%s` is allowed to appear in `formula`",
                 treat))
  }

  proj_data <- data.frame(x[["values"]]) |>
    setNames(treat) |>
    model.frame(formula = formula)

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  nm <- colnames(mm)

  if (!all(is.finite(mm))) {
    .err(sprintf("evaluation of the formula produced non-finite values of `%s`, which is not allowed",
                 treat))
  }

  if (length(by_list) > 1L) {
    mm <- block_diag(mm, length(by_list))
  }

  colnames(mm) <- NULL

  vcov_type <- x[["vcov_type"]]

  if (inherits(x, "adrf_est") && is_null(x$contrast)) {
    # transform
    transform_list <- process_transform(transform, .attr(x, "model")$family)
    transform <- transform_list$transform

    est <- transform(x$est)

    transformed <- !all(check_if_zero(est - x$est))
  }
  else {
    est <- x$est

    transformed <- FALSE
  }

  if (empirical) {
    treat_var <- .attr(x, "treat_var")
    by_id <- .attr(x, "by_id")

    ni <- length(treat_var)

    w <- rep_with(0, est)

    kw_mat <- matrix(0, nrow = ni, ncol = n)

    b_split <- gsplit(seq_len(ni), by_id)

    for (b in b_split) {
      kw_b <- get_kernel_w(treat_var[b], x[["values"]])

      kw_mat[b, ] <- kw_b
      w <- w + fsum(kw_b) / sum(kw_b)
    }

    w <- w * rep.int(get_trapezoidal_w(x[["values"]]),
                     length(by_list))
  }
  else {
    # Trapezoidal weights for intergal projection
    w <- rep.int(get_trapezoidal_w(x[["values"]]),
                 length(by_list))
  }

  fit <- lm.wfit(mm, y = est, w = w)

  out <- list(
    coefficients = fit$coefficients,
    residuals = fit$residuals,
    fitted.values = fit$fitted.values,
    x = mm,
    terms = mt
  )

  if (identical(vcov_type, "bootstrap")) {
    out$boot <- x$boot
    out$boot[["t0"]] <- fit$coefficients
    out$boot[["t"]] <- {
      if (transformed) lm.wfit(mm, y = t(transform(x$boot[["t"]])), w = w)$coefficients |> t()
      else lm.wfit(mm, y = t(x$boot[["t"]]), w = w)$coefficients |> t()
    }
    out$boot[["call"]] <- NULL

    V <- vcov(out$boot)
  }
  else {
    qr <- fit$qr
    ww <- qr.coef(qr, diag(sqrt(w)))

    if (transformed) {
      ww <- ww %*% diag(transform_list$d_transform(x$est))
    }

    V <- quad_mult(ww, vcov(x))

    # if (empirical && identical(vcov_type, "unconditional")) {
    #   S_proj <- tcrossprod(x$S, ww) + kw_mat %*% (fit$residuals * mm)
    #
    #   V <- crossprod(S_proj / nrow(S_proj))
    # }
  }

  dimnames(V) <- list(names(fit$coefficients),
                      names(fit$coefficients))

  out$vcov <- V

  out$coef_names <- nm
  out$vcov_type <- x[["vcov_type"]]
  out$call <- mcall

  out$by_grid <- x[["by_grid"]]
  out$contrast <- x[["contrast"]]

  attr(out, "by_id") <- .attr(x, "by_id")
  attr(out, "obj") <- x

  if (transformed) {
    attr(out, "inv_transform") <- transform_list$inv_transform
  }

  class(out) <- "curve_proj"

  out
}

curve_proj2 <- function(x, formula, transform = TRUE) {
  mcall <- process_call()

  chk::chk_not_missing(x, "`x`")

  .chk_is(x, c("adrf_est", "amef_est"))

  chk::chk_not_missing(formula, "`formula`")

  if (!rlang::is_formula(formula, lhs = FALSE)) {
    .err("`formula` must be a one-sided formula with the projection model on the right-hand side")
  }

  empirical <- FALSE
  chk::chk_flag(empirical)

  n <- length(x[["values"]])

  by_list <- {
    if (is_null(x$by_grid) && is_null(x$contrast)) list(seq_len(n))
    else if (is_not_null(x$contrast)) split(seq_along(x$est), rep(seq_along(x$contrast), each = n))
    else split(seq_along(x$est), rep(seq_row(x$by_grid), each = n))
  }

  #Check that only treat is in formula
  treat <- x[["treat"]]
  vars_in_formula <- get_varnames(formula)

  if (!all(get_varnames(formula) %in% treat)) {
    .err(sprintf("only the treatment variable `%s` is allowed to appear in `formula`",
                 treat))
  }

  proj_data <- data.frame(x[["values"]]) |>
    setNames(treat) |>
    model.frame(formula = formula)

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  nm <- colnames(mm)

  if (!all(is.finite(mm))) {
    .err(sprintf("evaluation of the formula produced non-finite values of `%s`, which is not allowed",
                 treat))
  }

  if (length(by_list) > 1L) {
    mm <- block_diag(mm, length(by_list))
  }

  colnames(mm) <- NULL

  vcov_type <- x[["vcov_type"]]

  family <- NULL

  if (inherits(x, "adrf_est") && is_null(x$contrast)) {
    family <- .attr(x, "model")$family
  }

  if (is_null(family)) {
    family <- gaussian()
  }

  est <- x$est

  if (empirical) {
    treat_var <- .attr(x, "treat_var")
    by_id <- .attr(x, "by_id")

    ni <- length(treat_var)

    w <- rep_with(0, est)

    kw_mat <- matrix(0, nrow = ni, ncol = n)

    b_split <- gsplit(seq_len(ni), by_id)

    for (b in b_split) {
      kw_b <- get_kernel_w(treat_var[b], x[["values"]])

      kw_mat[b, ] <- kw_b
      w <- w + fsum(kw_b) / sum(kw_b)
    }

    w <- w * rep.int(get_trapezoidal_w(x[["values"]]),
                     length(by_list))
  }
  else {
    # Trapezoidal weights for intergal projection
    w <- rep.int(get_trapezoidal_w(x[["values"]]),
                 length(by_list))
  }

  fit <- glm.fit(mm, y = est, weights = w,
                 family = family)

  out <- list(
    coefficients = fit$coefficients,
    # residuals = fit$residuals,
    fitted.values = fit$fitted.values,
    x = mm,
    terms = mt
  )

  if (identical(vcov_type, "bootstrap")) {
    out$boot <- x$boot
    out$boot[["t0"]] <- fit$coefficients
    out$boot[["t"]] <- apply(x$boot[["t"]], 2L, function(y) glm.fit(mm, y = y, weights = w, family = family, start = fit$coefficients)$coefficients)
    out$boot[["call"]] <- NULL

    V <- vcov(out$boot)
  }
  else {
    mu   <- fit$fitted.values
    eta  <- fit$linear.predictors

    # Derivative of mean wrt linear predictor: dmu/deta
    mu_eta <- family$mu.eta(eta)

    # Variance function diag(V)
    V_diag <- family$variance(mu)

    # Construct matrices D, V, and H
    D <- mu_eta * mm
    V_inv <- diag(1 / V_diag)
    H <- quad_mult(t(D), V_inv)

    # Extra covariance contribution from Sigma
    V <- quad_mult(solve(H, crossprod(D, V_inv)), vcov(x))

    # if (empirical && identical(vcov_type, "unconditional")) {
    #   S_proj <- tcrossprod(x$S, ww) + kw_mat %*% (fit$residuals * mm)
    #
    #   V <- crossprod(S_proj / nrow(S_proj))
    # }
  }

  dimnames(V) <- list(names(fit$coefficients),
                      names(fit$coefficients))

  out$vcov <- V

  out$coef_names <- nm
  out$vcov_type <- x[["vcov_type"]]
  out$call <- mcall

  out$by_grid <- x[["by_grid"]]
  out$contrast <- x[["contrast"]]

  attr(out, "by_id") <- .attr(x, "by_id")
  attr(out, "obj") <- x

  class(out) <- "curve_proj"

  out
}

#' @exportS3Method stats::coef curve_proj
coef.curve_proj <- function(x, ...) {
  x[["coefficients"]] |>
    setNames(.get_coef_labels(x))
}

#' @exportS3Method stats::vcov curve_proj
vcov.curve_proj  <- function(x, ...) {
  dimnames(x[["vcov"]]) <- rep.int(list(.get_coef_labels(x)), 2L)

  x[["vcov"]]
}

#' @exportS3Method print curve_proj
print.curve_proj <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  res <- .summary_internal(x, conf_level = 0, null = NULL,
                           simultaneous = FALSE)

  class(res) <- c("summary.curve_proj", class(res))

  print(res, digits = digits, ...)

  invisible(x)
}

#' @exportS3Method summary curve_proj
#' @rdname curve_proj
summary.curve_proj <- function(object, conf_level = 0.95, df = NULL, ci.type = "perc", ...) {

  vcov_type <- object[["vcov_type"]]

  if (is_null(vcov_type)) {
    conf_level <- 0
  }

  chk::chk_number(conf_level)
  chk::chk_range(conf_level, c(0, 1))

  null <- 0

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  if (do_inference) {
    if (vcov_type == "bootstrap") {
      chk::chk_string(ci.type)
    }
    else {
      ci.type <- "wald"
    }

    # Process df
    if (is_null(df)) {
      df <- .attr(object, "obj") |> .attr("model") |> get_df()
    }
    else {
      chk::chk_number(df)
      chk::chk_gt(df, 2)
    }
  }

  transform_list <- process_transform()

  res <- .summary_internal(object, conf_level, null = null,
                           simultaneous = FALSE,
                           transform_list = transform_list,
                           df = df,
                           ci.type = ci.type)

  class(res) <- c("summary.curve_proj", class(res))

  res
}

#' @exportS3Method print summary.curve_proj
print.summary.curve_proj <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  chk::chk_whole_number(digits)

  .print_summary_internal(x = x, digits = digits, topn = Inf,
                          ...)

  invisible(x)
}

#' @exportS3Method plot curve_proj
#' @rdname curve_proj
plot.curve_proj <- function(x, ...) {

  proj <- .attr(x, "obj")
  proj$est <- x$fitted.values
  proj$vcov <- quad_mult(x$x, vcov(x))
  proj$boot <- x$boot

  if (is_not_null(.attr(x, "inv_transform"))) {
    transform <- list(transform = identity,
                      inv_transform = .attr(x, "inv_transform"),
                      d_transform = identity)

    plot(x = proj, transform = transform, ...)
  }
  else {
    plot(x = proj, ...)
  }
}

.get_coef_labels <- function(x) {
  if (is_null(x$by_grid) && is_null(x$contrast)) {
    return(x[["coef_names"]])
  }

  if (is_not_null(x$contrast)) {
    labels <- x$contrast
  }
  else {
    labels <- do.call(function(...) paste(..., sep = ", "),
                      lapply(names(x$by_grid), function(i) {
                        sprintf("%s = %s",
                                i,
                                add_quotes(x$by_grid[[i]], chk::vld_character_or_factor(x$by_grid[[i]])))
                      }))
  }

  sprintf("%s | %s",
          rep(labels, each = length(x[["coef_names"]])),
          rep.int(x[["coef_names"]], length(labels)))
}
