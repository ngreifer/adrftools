#' Project an effect curve onto a simpler model
#'
#' @description
#' `curve_projection()` produces a projection of an estimated effect curve onto a specified linear model that is a function only of the treatment to act as a more interpretable summary of the original effect curve.
#'
#' @inheritParams summary.curve_est
#' @param x an [`effect_curve`] object; the output of a call to [adrf()].
#' @param model the projection model to be fit. Can be a one-sided formula corresponding to the projection model or one of the following strings: `"flat"`, `"linear"`, `"quadratic"`, `"cubic"`.
#' @param transform whether to compute the projection using a transformation of the linear predictor. Allowable options include `TRUE`, `FALSE`, or a function specifying a transformation (of which the inverse is used as the inverse link of the projection model). Ignored unless `object` is an ADRF. See Details.
#' @param subset an optional logical expression indicating the subset of subgroups for which to compute the projection. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param object,object2 a `curve_projection` object; the output of a call to `curve_projection()`.
#' @param null the null value for hypothesis test. Default is 0. Set to `NA` to omit tests.
#' @param df the "denominator" degrees of freedom to use for the test. Default is to use the residual degrees of freedom from the original model if it is a linear model (in which case an F-test is used) and `Inf` otherwise (in which case a \eqn{\chi^2} test is used).
#' @param \dots ignored.
#'
#' @returns
#' `curve_projection()` returns an `curve_projection` object, which inherits from [`effect_curve`]. This object is a function that produces estimates of the effect curve projection when called with values of the treatment as inputs. See [`effect_curve`] for details on calling this function.
#'
#' The coefficients and covariance matrix of the fitted projection should be extracted with `coef()` and `vcov()`, respectively. `summary()` produces the coefficients and quantities of interest for inference (test statistics, p-values, and confidence intervals). Using `plot()` on a `curve_projection` object plots the projection curve as if it was an effect curve; see [plot.effect_curve()] for adding the projection curve to the plot of the original effect curve.
#'
#' @details
#' The projection model can be thought of as a linear regression of the effect curve estimates on the treatment. Whereas the original effect curve may be complicated and nonlinear, the projection model can be simple and easily interpretable, though it must be understood as a summary of the original effect curve. For example, the original ADRF might have been computed from an outcome model that involves treatment splines, covariates, and treatment-covariate interactions. Though the ADRF is a univariable function (i.e., of only the treatment), it isn't described by a single set of parameters. The linear projection of the ADRF, though, could be a simple linear model, described by an intercept and the slope on treatment. Though only a rough approximation to the ADRF, the linear projection may be more easily interpreted. This concept is described in Neugebauer and van der Laan (2007).
#'
#' `curve_projection()` fits this projection model and accounts for the uncertainty in the estimates of the effect curve in computing uncertainty estimates for the projection model parameters. Because the true effect curve is continuous, the model is fit minimizing
#' \deqn{\int_{a_\text{lower}}^{a_\text{upper}}{\left(\hat{\theta}(a)-\hat{\mathbf{\beta}} B(a) \right)^2 da}}
#' where \eqn{\hat{\theta}(a)} is the effect curve estimate at treatment value \eqn{a}, \eqn{B(a)} is the basis function representation of \eqn{a} (i.e., as specified in `model`), and \eqn{\hat{\mathbf{\beta}}} is the vector of projection parameters to be estimated. This integral is approximated using a trapezoidal Riemann sum over the effect curve grid points.
#'
#' The covariance of the projection parameters can be computed using the delta method applied to the estimated covariance of the original effect curve estimates. When bootstrapping or posterior inference are used, the projection is applied to each bootstrap or posterior draw, respectively.
#'
#' ## Transform
#'
#' When `transform` is specified, the projection minimizes the distance between the original effect curve and the transformed linear predictor; that is, it minimizes
#'
#' \deqn{\int_{a_\text{lower}}^{a_\text{upper}}{\left(\hat{\theta}(a)-f^{-1}\left(\hat{\mathbf{\beta}} B(a) \right) \right)^2 da}}
#'
#' where \eqn{f^{-1}(y)} is the inverse of the transformation supplied to `transform` (i.e., corresponding to the inverse link function of a generalized linear model), essentially using nonlinear least squares (NLS) to estimate the effect curve projection. This make the coefficients in the projection model correspond to the coefficients on the linear predictor \eqn{\hat{\mathbf{\beta}} B(a)}. In this case, the projection is not simply a linear projection, but it may still be more interpretable than the original ADRF. For example, if the outcome model was originally fit using logistic regression and `transform = TRUE` in the call to `curve_projection()` with `model = "linear"`, the resulting projection would be a logistic curve governed by the intercept and slope of the linear predictor. See Examples for an example of this.
#'
#' By default, `transform` is `TRUE`, which means that when the original outcome model had a `family` component (e.g., a generalized linear model) and an ADRF is supplied to `curve_projection()`, the link is automatically supplied to `transform` and the projection model will be a nonlinear function of the linear predictor. Set `transform` to `FALSE` to require that the projection curve be simply the linear predictor with no transformation. Note this can lead to invalid estimates when the outcome is bounded.
#'
#' ## Comparing projection models
#'
#' `anova()` performs a Wald test comparing two nested projection models. The null hypothesis is that the simpler model is sufficient, i.e., that the coefficients on the terms in the larger model (supplied to `object`) that are not in the smaller model (supplied to `object2`) are all zero. Rejecting the null hypothesis implies that the larger model fits better.
#'
#' @seealso
#' * [plot.effect_curve()] for plotting the effect curve and its projection
#' * [`effect_curve`] for computing point estimates along the effect curve projection
#' * [summary.effect_curve()] for testing hypotheses about the effect curve, such as whether a given projection is sufficient
#' * [anova()] for comparing linear models
#'
#' @references
#' Neugebauer, R., & van der Laan, M. (2007). Nonparametric causal effects based on marginal structural models. *Journal of Statistical Planning and Inference*, 137(2), 419–434. \doi{/10.1016/j.jspi.2005.12.008}
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
#' # Linear projection is sufficient for
#' # characterizing the ADRF
#' summary(adrf1, hypothesis = "linear")
#'
#' # Compute the linear projection
#' proj1 <- curve_projection(adrf1, "linear")
#' # proj1 <- curve_projection(adrf1, ~logBLL) #same model
#'
#' proj1
#'
#' # Coefficients of the projection model
#' coef(proj1)
#' summary(proj1)
#'
#' # Plot the projection
#' plot(proj1)
#'
#' # Plot the projection over the ADRF
#' plot(adrf1, proj = proj1)
#'
#' # Compute a cubic projection
#' proj2 <- curve_projection(adrf1, "cubic")
#' # proj2 <- curve_projection(adrf1, ~poly(logBLL, 3)) #same model
#'
#' # Compare cubic to linear projection
#' anova(proj2, proj1)

#' @export
curve_projection <- function(x, model, transform = TRUE) {
  chk::chk_not_missing(x, "`x`")
  check_effect_curve(x, projection_ok = FALSE)

  chk::chk_not_missing(model, "`model`")

  .treat <- .attr(x, ".treat")

  if (chk::vld_string(model)) {
    model <- model |>
      match_arg(c("flat", "linear", "quadratic", "cubic")) |>
      switch(flat = ~1,
             linear = as.formula(sprintf("~ %1$s", .treat)),
             quadratic = as.formula(sprintf("~ %1$s + I(%1$s^2)", .treat)),
             cubic = as.formula(sprintf("~ %1$s + I(%1$s^2) + I(%1$s^3)", .treat)))
  }
  else if (!rlang::is_formula(model, lhs = FALSE)) {
    .err("{.arg model} must be a string or a one-sided formula with the projection model on the right-hand side")
  }
  else if (!all(get_varnames(model) %in% .treat)) {
    #Check that only treat is in model
    .err("only the treatment variable `{(.treat)}` is allowed to appear in {.arg model}")
  }

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")
  .vcov_type <- .attr(x, ".vcov_type")

  proj_data <- data.frame(.values) |>
    setNames(.treat) |>
    model.frame(formula = model)

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  nm <- colnames(mm)

  if (!all(is.finite(mm))) {
    .err("evaluation of the projection model produced non-finite values of `{(.treat)}`, which is not allowed")
  }

  n_by <- get_n_by(.contrast, .by_grid)

  mm <- block_diag(mm, n_by)

  colnames(mm) <- NULL

  transform_list <- process_transform(transform, x, .est)

  # Trapezoidal weights for integral projection
  w <- rep.int(get_trapezoidal_w(.values),
               n_by)

  fit <- .nls_fit(x = mm, y = .est, w = w,
                  transform = transform_list,
                  jac = .vcov_type != "none")

  proj_coefficients <- fit$coefficients

  if (.vcov_type == "none") {
    proj_fn <- .make_fun(x,
                         .est = fit$fitted.values)
  }
  else {
    if (.vcov_type == "bootstrap") {
      proj_boot <- .attr(x, ".boot")

      proj_boot[["t0"]] <- fit$coefficients
      proj_boot[["t"]] <- proj_boot[["t"]] |>
        dapply(MARGIN = 1L, return = "matrix",
               FUN = function(esti) {
                 .nls_fit(x = mm, y = esti, w = w,
                          transform = transform_list,
                          start = proj_coefficients)$coefficients
               }) |>
        unname()

      proj_boot[["call"]] <- NULL

      proj_vcov <- stats::vcov(proj_boot)

      .boot <- .attr(x, ".boot")

      .boot[["t0"]] <- fit$fitted.values
      .boot[["t"]] <- tcrossprod(proj_boot[["t"]], mm) |>
        transform_list$inv_transform()

      .boot[["call"]] <- NULL

      .vcov <- stats::vcov(.boot)

      proj_fn <- .make_fun(x,
                           .est = fit$fitted.values,
                           .vcov = .vcov,
                           .boot = .boot)

      attr(proj_fn, ".proj_boot") <- proj_boot
    }
    else if (.vcov_type == "posterior") {
      proj_draws <- .attr(x, ".draws") |>
        dapply(MARGIN = 1L, return = "matrix",
               FUN = function(esti) {
                 .nls_fit(x = mm, y = esti, w = w,
                          transform = transform_list,
                          start = proj_coefficients)$coefficients
               }) |>
        unname()

      proj_vcov <- stats::cov(proj_draws)

      .draws <- tcrossprod(proj_draws, mm) |>
        transform_list$inv_transform()

      .vcov <- stats::vcov(.draws)

      proj_fn <- .make_fun(x,
                           .est = fit$fitted.values,
                           .vcov = .vcov,
                           .draws = .draws)

      attr(proj_fn, ".proj_draws") <- proj_draws
    }
    else {
      J  <- fit$jacobian
      JW <- J * w

      ww <- solve(crossprod(J, JW), t(JW))

      proj_vcov <- quad_mult(ww, .attr(x, ".vcov"))

      .vcov <- quad_mult(J, proj_vcov)

      proj_fn <- .make_fun(x,
                           .est = fit$fitted.values,
                           .vcov = .vcov)
    }

    attr(proj_fn, ".proj_vcov") <- proj_vcov
  }

  attr(proj_fn, ".proj_coefficients") <- proj_coefficients
  attr(proj_fn, ".proj_coefficient_names") <- nm
  attr(proj_fn, ".proj_formula") <- deparse1(mt)

  class(proj_fn) <- unique(c("curve_projection", class(proj_fn)))

  proj_fn
}

#' @exportS3Method summary curve_projection
#' @rdname curve_projection
summary.curve_projection <- function(object, conf_level = 0.95, null = 0,
                                     df = NULL, ci.type = "perc", subset = NULL, ...) {
  .vcov_type <- .attr(object, ".vcov_type")
  .contrast <- .attr(object, ".contrast")
  .by_grid <- .attr(object, ".by_grid")

  .proj_coefficients <- stats::coef(object)
  .proj_vcov <- stats::vcov(object)
  .proj_boot <- .attr(object, ".proj_boot")
  .proj_draws <- .attr(object, ".proj_draws")
  .proj_coefficient_names <- .attr(object, ".proj_coefficient_names")

  if (.vcov_type == "none") {
    conf_level <- 0
    null <- NA
  }
  else {
    conf_level <- process_conf_level(conf_level)
    null <- process_null(null)
  }

  # subset
  if (is_null(.by_grid) || is_not_null(.contrast)) {
    .s <- seq_along(.proj_coefficients)
  }
  else {
    .subset <- process_subset_by_grid(substitute(subset),
                                      .by_grid = .by_grid,
                                      .contrast = .contrast)

    .s <- which(rep(.subset, each = length(.proj_coefficients) / nrow(.by_grid)))

    .by_grid <- ss(.by_grid, .subset)

    .proj_coefficients <- .proj_coefficients[.s]
  }

  # Only do inference if vcov_type is not empty and conf_level is
  # nonzero or null hypothesis is specified
  do_inference <- (.vcov_type != "none") &&
    is_not_null(.proj_vcov) &&
    (conf_level > 0 || !anyNA(null))

  if (do_inference) {
    ci.type <- process_ci.type(ci.type, .vcov_type)

    if (ci.type == "wald") {
      # Process df
      df <- df %or% .attr(object, ".df")

      chk::chk_number(df)
      chk::chk_gt(df, 0)

      stat <- if (is.finite(df)) "t" else "z"

      res_names <- c("term", "estimate",
                     "std.error",
                     stat[!anyNA(null)],
                     "p.value"[!anyNA(null)],
                     "conf.low"[conf_level > 0],
                     "conf.high"[conf_level > 0])

      res <- make_df(res_names, length(.s))

      se <- sqrt(diag(ss(.proj_vcov, .s, .s)))
      zeros <- se < 1e-10

      res$std.error <- se
    }
    else {
      res_names <- c("term", "estimate",
                     "p.value"[!anyNA(null)],
                     "conf.low"[conf_level > 0],
                     "conf.high"[conf_level > 0])

      res <- make_df(res_names, length(.s))
    }

    t_crit <- NULL

    # Compute p-values
    if (!anyNA(null)) {
      if (ci.type == "wald") {
        # Compute test statistic on transformed estimates if null is specified
        res[[stat]][!zeros] <- (.proj_coefficients[!zeros] - null) / se[!zeros]

        res$p.value[!zeros] <- switch(stat,
                                      "z" = 2 * pnorm(-abs(res[[stat]][!zeros])),
                                      "t" = 2 * pt(-abs(res[[stat]][!zeros]), df))
      }
      else if (.vcov_type == "bootstrap") {
        s <- summary(.proj_boot,
                     conf = 0,
                     parm = .s,
                     ci.type = ci.type,
                     p.value = TRUE,
                     null = null,
                     simultaneous = FALSE)

        res$p.value <- s[, ncol(s)]
      }
      else if (.vcov_type == "posterior") {
        res$p.value <- posterior_p_value(.proj_draws, parm = .s, null = null,
                                         simultaneous = FALSE)
      }
    }

    # Compute CIs
    if (conf_level > 0) {
      if (ci.type == "wald") {
        t_crit <- switch(stat,
                         "z" = abs(qnorm((1 - conf_level) / 2)),
                         "t" = abs(qt((1 - conf_level) / 2, df)))

        # Reverse transformation
        res$conf.low[] <- .proj_coefficients - t_crit * se
        res$conf.high[] <- .proj_coefficients + t_crit * se
      }
      else if (.vcov_type == "bootstrap") {
        boot_ci <- confint(.proj_boot,
                           level = conf_level,
                           parm = .s,
                           ci.type = ci.type,
                           simultaneous = FALSE)

        res$conf.low[] <- boot_ci[, 1L]
        res$conf.high[] <- boot_ci[, 2L]
      }
      else if (.vcov_type == "posterior") {
        draws_ci <- posterior_ci(.proj_draws,
                                 level = conf_level,
                                 parm = .s,
                                 simultaneous = FALSE)

        res$conf.low[] <- draws_ci[, 1L]
        res$conf.high[] <- draws_ci[, 2L]
      }
    }
  }
  else {
    res_names <- c("term", "estimate")

    res <- make_df(res_names, length(.s))
  }

  res$estimate <- .proj_coefficients

  res <- add_est_labels(res, .contrast, .by_grid,
                        .proj_coefficient_names, "term")

  attr(res, "conf_level") <- conf_level

  if (do_inference) {
    if (conf_level > 0 && is_not_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(res, "crit") <- t_crit
    }

    attr(res, "ci.type") <- ci.type

    if (!anyNA(null)) {
      attr(res, "null") <- null
    }
  }

  attr(res, ".curve_type") <- get_curve_type(object)
  attr(res, ".contrast") <- .attr(object, ".contrast")
  attr(res, ".by_grid") <- .attr(object, ".by_grid")
  attr(res, ".vcov_type") <- .attr(object, ".vcov_type")

  attr(res, ".proj_coefficients") <- .proj_coefficients
  attr(res, ".proj_vcov") <- .proj_vcov
  attr(res, ".proj_coefficient_names") <- .proj_coefficient_names

  class(res) <- unique(c("summary.curve_projection", class(res)))

  res
}

#' @exportS3Method print summary.curve_projection
print.summary.curve_projection <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  chk::chk_whole_number(digits)

  .reference <- .attr(x, ".reference")
  .contrast <- .attr(x, ".contrast")

  main <- paste(c(get_curve_type(x),
                  if (is_not_null(.reference)) "Reference",
                  if (is_not_null(.contrast)) "Contrast",
                  "Projection Coefficients"),
                collapse = " ")

  .print_estimate_table(x = x, digits = digits, topn = Inf,
                        main = main,
                        rownames = FALSE,
                        ...)

  .print_inference(x)

  invisible(x)
}

#' @exportS3Method stats::coef curve_projection
#' @rdname curve_projection
coef.curve_projection <- function(object, ...) {
  .attr(object, ".proj_coefficients") |>
    setNames(.get_curve_est_labels(object))
}

#' @exportS3Method stats::vcov curve_projection
#' @rdname curve_projection
vcov.curve_projection  <- function(object, ...) {
  v <- .attr(object, ".proj_vcov")

  dimnames(v) <- rep.int(list(.get_curve_est_labels(object)), 2L)

  v
}

#' @exportS3Method stats::coef summary.curve_projection
coef.summary.curve_projection <- coef.curve_projection

#' @exportS3Method stats::vcov summary.curve_projection
vcov.summary.curve_projection <- vcov.curve_projection

#' @exportS3Method stats::model.matrix curve_projection
model.matrix.curve_projection <- function(object, ...) {
  .treat <- .attr(object, ".treat")
  .values <- .attr(object, ".values")
  .by_grid <- .attr(object, ".by_grid")
  .contrast <- .attr(object, ".contrast")

  proj_data <- data.frame(.values) |>
    setNames(.treat) |>
    model.frame(formula = formula(object))

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  if (!all(is.finite(mm))) {
    .err("evaluation of the projection model produced non-finite values of `{(.treat)}`, which is not allowed")
  }

  n_by <- get_n_by(.contrast, .by_grid)

  mm <- block_diag(mm, n_by)

  colnames(mm) <- .get_curve_est_labels(object)

  mm
}

#' @exportS3Method stats::formula curve_projection
formula.curve_projection <- function(x, ...) {
  .attr(x, ".proj_formula") |> as.formula()
}

#' @exportS3Method stats::anova curve_projection
#' @rdname curve_projection
anova.curve_projection <- function(object, object2, df = NULL, ...) {
  chk::chk_not_missing(object, "`object`")
  .chk_is(object, "curve_projection")

  chk::chk_not_missing(object2, "`object2`")
  .chk_is(object2, "curve_projection")

  if (any_apply(c(".values", ".treat", ".vcov_type", ".curve_type", ".response",
                  ".by_grid", ".contrast", ".reference", ".family", ".df"),
                function(a) !identical(.attr(object, a), .attr(object2, a)))) {
    .err("{.arg object} and {.arg object2} must be the outputs of {.fun curve_projection} applied to the same `effect_curve` object")
  }

  # Process df
  df <- df %or% .attr(object, ".df")

  chk::chk_number(df)
  chk::chk_gt(df, 0)

  test <- if (is.finite(df)) "F" else "Chisq"

  tolerance <- ...get("tolerance", 1e-7)

  chk::chk_number(tolerance)
  chk::chk_gt(tolerance, 0)

  b1 <- coef(object, complete = FALSE)
  b2 <- coef(object2, complete = FALSE)

  if (length(b2) >= length(b1)) {
    .err("{.arg object2} does not appear to be nested within {.arg object1}")
  }

  Z1 <- .lm.fit(x = model.matrix(object2)[, names(b2), drop = FALSE],
                y = model.matrix(object)[, names(b1), drop = FALSE])[["residuals"]]

  Z1_svd <- svd(Z1)
  keep <- Z1_svd[["d"]] > tolerance
  .q <- sum(keep)

  if (.q > length(b1) - length(b2)) {
    .err("`object2` does not appear to be nested within `object`")
  }

  L <- t(Z1_svd[["v"]][, keep, drop = FALSE])

  V <- stats::vcov(object)

  value.hyp <- L %*% b1
  vcov.hyp <- quad_mult(L, V)

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  width <- min(55L, getOption("width", 55L))

  .title <- "Wald Test"

  .topnote <- sprintf("Model 1: %s\nModel 2: %s\n",
                      deparse1(formula(object)) |>
                        strwrap(width, exdent = 11L),
                      deparse1(formula(object2)) |>
                        strwrap(width, exdent = 11L))

  .nullhyp <- "H\u2080: Extra terms in Model 1 have coefficients equal to 0 (i.e., Model 2 is sufficient)"

  result <- make_df(c("df", test, "p.value"), "")

  result[[1L]] <- as.integer(.q)
  if (test == "Chisq") {
    result[[2L]] <- SSH
    result[[3L]] <- pchisq(result[[2L]], .q, lower.tail = FALSE)
  }
  else {
    result[[2L]] <- SSH / .q
    result[[3L]] <- pf(result[[2L]], .q, df, lower.tail = FALSE)
  }


  attr(result, "heading") <- c(.title, .topnote, .nullhyp)
  attr(result, "value") <- value.hyp
  attr(result, "vcov") <- vcov.hyp
  attr(result, ".vcov_type") <- .attr(object, ".vcov_type")
  attr(result, ".df") <- df

  class(result) <- c("anova.curve_projection", "anova", class(result))

  result
}

#' @exportS3Method print anova.curve_projection
print.anova.curve_projection <- function(x, digits = max(4L, getOption("digits") - 3L), ...) {
  width <- min(55L, getOption("width", 55L))

  out <- c(center_just(attr(x, "heading")[1L], wrt = space(width)),
           txtbar(width),
           attr(x, "heading")[2L],
           cli::ansi_strwrap(attr(x, "heading")[3L],
                             width, exdent = 4L),
           "")

  cat(out, sep = "\n")

  .print_estimate_table(x, digits = digits, topn = Inf, bar = FALSE)

  .df <- .attr(x, ".df")

  cat(c(txtbar(width),
        .it(sprintf("Inference: %s%s",
                    .get_vcov_type_name(attr(x, ".vcov_type")),
                    if (is.finite(.df)) sprintf(" (df = %s)", round(.df, 1L)) else "")) |>
          cli::ansi_strwrap(width, exdent = 2L)),
      sep = "\n")

  invisible(x)
}

.nls_fit <- function(x, y, w = NULL, offset = 0, transform, start = NULL,
                     maxit = 50, tol = 1e-8, jac = FALSE) {

  f <- transform$transform
  finv <- transform$inv_transform
  finv_prime <- transform$d_inv_transform

  n <- nrow(x)

  if (is_null(w)) {
    wsqt <- rep.int(1, n)
  }
  else {
    wsqt <- sqrt(w)
  }

  fy <- f(y)

  # If identity transform, just use WLS
  if (all(check_if_zero(fy - y))) {
    fit <- .lm.fit(x = wsqt * x, y = wsqt * (fy - offset))

    out <- list(
      coefficients = fit$coefficients,
      fitted.values = drop(x %*% fit$coefficients) + offset,
      jacobian = if (jac) x
    )

    return(out)
  }

  if (is_null(start)) {
    start <- .lm.fit(x = wsqt * x, y = wsqt * (fy - offset))$coefficients
  }

  .theta <- start

  for (iter in seq_len(maxit)) {
    eta <- drop(x %*% .theta)
    mu  <- finv(eta)

    J <- finv_prime(eta) * x

    .step <- .lm.fit(x = wsqt * J, y = wsqt * (y - mu))$coefficients

    .theta <- .theta + .step

    if (max(abs(.step)) < tol) {
      break
    }
  }

  .eta <- drop(x %*% .theta)

  list(
    coefficients = .theta,
    fitted.values = finv(.eta),
    jacobian = if (jac) finv_prime(.eta) * x
  )
}
