#' Project an effect curve onto a simpler model
#'
#' @description
#' `curve_projection()` produces a projection of an estimated effect curve onto a specified linear model that is a function only of the treatment to act as a more interpretable summary of the original effect curve.
#'
#' @inheritParams summary.curve_est
#' @param x an [`effect_curve`] object; the output of a call to [adrf()].
#' @param model the projection model to be fit. Can be a one-sided formula corresponding to the projection model or one of the following strings: `"flat"`, `"linear"`, `"quadratic"`, `"cubic"`.
#' @param transform whether to compute the projection using a transformation of the linear predictor. Allowable options include `TRUE`, `FALSE`, a family object (e.g., `binomial()`), a link object (e.g., the output of [make.link()]), or one of the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`; the inverse of the transformation is used as the inverse link of the projection model. Ignored unless `x` is an ADRF. See Details.
#' @param integrate string; how treatment values are weighted in fitting the projection. Allowable options include `"uniform"` (the default), which weights all treatment values in the range of the effect curve equally, and `"empirical"`, which weights them by an estimate of the density of the treatment. The weighting defines the projection, so the coefficients differ between the two options. Only `"uniform"` can be used with a shift curve (see [shift_curve()]). See Details.
#' @param subset an optional logical expression indicating the subset of subgroups for which to compute the projection. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param outcome an optional character vector naming the categories of the outcome for which to compute the projection. Can only be used when the model gives a separate prediction for each category of the outcome (e.g., multinomial and ordinal models); see the section *Categorical outcomes* in Details of [adrf()].
#' @param object,object2 a `curve_projection` object; the output of a call to `curve_projection()`.
#' @param null the null value for hypothesis test. Default is 0. Set to `NA` to omit tests.
#' @param df the "denominator" degrees of freedom to use for the test. Default is to use the residual degrees of freedom from the original model if it is a linear model (in which case an F-test is used) and `Inf` otherwise (in which case a \eqn{\chi^2} test is used), or, for models fit to multiply imputed data, degrees of freedom that account for the number of imputations (see [adrf()]).
#' @param \dots ignored.
#'
#' @returns
#' `curve_projection()` returns a `curve_projection` object, which inherits from [`effect_curve`]. This object is a function that produces estimates of the effect curve projection when called with values of the treatment as inputs. See [`effect_curve`] for details on calling this function.
#'
#' The coefficients and covariance matrix of the fitted projection should be extracted with `coef()` and `vcov()`, respectively. `summary()` produces the coefficients and quantities of interest for inference (test statistics, p-values, and confidence intervals). Using `plot()` on a `curve_projection` object plots the projection curve as if it was an effect curve; see [plot.effect_curve()] for adding the projection curve to the plot of the original effect curve.
#'
#' @details
#' The projection model can be thought of as a linear regression of the effect curve estimates on the treatment. Whereas the original effect curve may be complicated and nonlinear, the projection model can be simple and easily interpretable, though it must be understood as a summary of the original effect curve. For example, the original ADRF might have been computed from an outcome model that involves treatment splines, covariates, and treatment-covariate interactions. Though the ADRF is a univariable function (i.e., of only the treatment), it isn't described by a single set of parameters. The linear projection of the ADRF, though, could be a simple linear model, described by an intercept and the slope on treatment. Though only a rough approximation to the ADRF, the linear projection may be more easily interpreted. This concept is described in Neugebauer and van der Laan (2007).
#'
#' `curve_projection()` fits this projection model and accounts for the uncertainty in the estimates of the effect curve in computing uncertainty estimates for the projection model parameters. Because the true effect curve is continuous, the model is fit minimizing
#' \deqn{\int_{a_\text{lower}}^{a_\text{upper}}{\left(\hat{\theta}(a)- B(a)\hat{\mathbf{\beta}} \right)^2 da}}
#' where \eqn{\hat{\theta}(a)} is the effect curve estimate at treatment value \eqn{a}, \eqn{B(a)} is the basis function representation of \eqn{a} (i.e., as specified in `model`), and \eqn{\hat{\mathbf{\beta}}} is the vector of projection parameters to be estimated. This integral is approximated using a trapezoidal Riemann sum over the effect curve grid points.
#'
#' The covariance of the projection parameters can be computed using the delta method applied to the estimated covariance of the original effect curve estimates. When bootstrapping or posterior inference are used, the projection is applied to each bootstrap or posterior draw, respectively.
#'
#' ## Transform
#'
#' When `transform` is specified, the projection minimizes the distance between the original effect curve and the transformed linear predictor; that is, it minimizes
#'
#' \deqn{\int_{a_\text{lower}}^{a_\text{upper}}{\left(\hat{\theta}(a)-f^{-1}\left( B(a) \hat{\mathbf{\beta}} \right) \right)^2 da}}
#'
#' where \eqn{f^{-1}(y)} is the inverse of the transformation supplied to `transform` (i.e., corresponding to the inverse link function of a generalized linear model), essentially using nonlinear least squares (NLS) to estimate the effect curve projection. This makes the coefficients in the projection model correspond to the coefficients on the linear predictor \eqn{B(a)\hat{\mathbf{\beta}}}. In this case, the projection is not simply a linear projection, but it may still be more interpretable than the original ADRF. For example, if the outcome model was originally fit using logistic regression and `transform = TRUE` in the call to `curve_projection()` with `model = "linear"`, the resulting projection would be a logistic curve governed by the intercept and slope of the linear predictor. See Examples for an example of this.
#'
#' By default, `transform` is `TRUE`, which means that when the original outcome model had a `family` component (e.g., a generalized linear model) and an ADRF is supplied to `curve_projection()`, the link is automatically supplied to `transform` and the projection model will be a nonlinear function of the linear predictor. Set `transform` to `FALSE` to require that the projection curve be simply the linear predictor with no transformation. Note this can lead to invalid estimates when the outcome is bounded.
#'
#' ## Weighting treatment values
#'
#' With `integrate = "empirical"`, the integrals above are weighted by the density of the treatment, so the projection is the best approximation to the effect curve over the distribution of the treatment rather than over its range, and parts of the curve where the treatment is common count more than parts where it is sparse. The density is estimated as described in [summary.effect_curve()]. Because the weighting changes the projection itself, `anova()` requires both projections to use the same `integrate`.
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
#' Neugebauer, R., & van der Laan, M. (2007). Nonparametric causal effects based on marginal structural models. *Journal of Statistical Planning and Inference*, 137(2), 419–434. \doi{10.1016/j.jspi.2005.12.008}
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
curve_projection <- function(x, model, transform = TRUE, integrate = "uniform") {
  arg::arg_supplied(x)

  check_effect_curve(x, projection_ok = FALSE)

  arg::arg_supplied(model)

  .treat <- .attr(x, ".treat")

  if (rlang::is_string(model)) {
    model <- model |>
      arg::match_arg(c("flat", "linear", "quadratic", "cubic")) |>
      switch(flat = ~1,
             linear = as.formula(sprintf("~ %1$s", .treat)),
             quadratic = as.formula(sprintf("~ %1$s + I(%1$s^2)", .treat)),
             cubic = as.formula(sprintf("~ %1$s + I(%1$s^2) + I(%1$s^3)", .treat)))
  }
  else if (!rlang::is_formula(model, lhs = FALSE)) {
    arg::err("{.arg model} must be a string or a one-sided formula with the projection model on the right-hand side")
  }
  else if (!all(get_varnames(model) %in% .treat)) {
    #Check that only treat is in model
    arg::err("only {.var {(.treat)}} is allowed to appear in {.arg model}")
  }

  integrate <- process_integrate(integrate, x)

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")
  .vcov_type <- .attr(x, ".vcov_type")

  proj_data <- qDF(.values) |>
    setNames(.treat) |>
    model.frame(formula = model, na.action = stats::na.pass)

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  nm <- colnames(mm)

  if (!all(is.finite(mm))) {
    arg::err("evaluating the projection model at the values of {.var {(.treat)}} produced non-finite values, which is not allowed")
  }

  .check_proj_mm(mm)

  n_by <- get_n_by(.contrast, .by_grid)

  mm <- block_diag(mm, n_by)

  colnames(mm) <- NULL

  transform_list <- process_transform(transform, x, .est)

  # Weights of the grid points for the integral projection, one block per
  # subgroup or contrast; with the block-diagonal model matrix, fitting all
  # blocks together is the same as fitting each with its own weights
  w <- get_grid_w(.values, integrate, .attr(x, ".treat_density"), seq_len(n_by)) |>
    unlist()

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
      proj_boot[["t"]] <- .proj_coef_draws(mm, proj_boot[["t"]], w,
                                           transform_list, proj_coefficients)

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
      proj_draws <- .proj_coef_draws(mm, .attr(x, ".draws"), w,
                                     transform_list, proj_coefficients)

      proj_vcov <- stats::cov(proj_draws)

      .draws <- tcrossprod(proj_draws, mm) |>
        transform_list$inv_transform()

      .vcov <- stats::cov(.draws)

      proj_fn <- .make_fun(x,
                           .est = fit$fitted.values,
                           .vcov = .vcov,
                           .draws = .draws)

      attr(proj_fn, ".proj_draws") <- proj_draws
    }
    else {
      J <- fit$jacobian

      ww <- .proj_influence(mm, .est, fit, w, transform_list)

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

  # The transformed estimates identify the transformation, so anova() can check
  # that two projections used the same one
  attr(proj_fn, ".proj_transformed_est") <- transform_list$transform(.est)
  attr(proj_fn, ".proj_formula") <- deparse1(mt)
  attr(proj_fn, ".proj_integrate") <- integrate

  class(proj_fn) <- unique(c("curve_projection", class(proj_fn)))

  proj_fn
}

#' @exportS3Method summary curve_projection
#' @rdname curve_projection
summary.curve_projection <- function(object, conf_level = 0.95, null = 0,
                                     df = NULL, ci.type = "perc", subset = NULL,
                                     outcome = NULL, ...) {
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

  # subset and outcome; process_subset_by_grid() gives an error when subset is
  # supplied without by or for a contrast
  .subset <- process_subset_by_grid(substitute(subset),
                                    .by_grid = .by_grid,
                                    .contrast = .contrast,
                                    outcome = outcome,
                                    .outcome = .attr(object, ".outcome"))

  if (all(.subset)) {
    .s <- seq_along(.proj_coefficients)
  }
  else {
    .s <- which(rep(.subset, each = length(.proj_coefficients) / length(.subset)))

    .by_grid <- ss(.by_grid, .subset)

    if (!is_null(.contrast)) {
      .contrast <- .contrast[.subset]
    }

    .proj_coefficients <- .proj_coefficients[.s]

    if (!is_null(.proj_vcov)) {
      .proj_vcov <- ss(.proj_vcov, .s, .s)
    }
  }

  # Only do inference if vcov_type is not empty and conf_level is
  # nonzero or null hypothesis is specified
  do_inference <- (.vcov_type != "none") &&
    !is_null(.proj_vcov) &&
    (conf_level > 0 || !anyNA(null))

  if (do_inference) {
    ci.type <- process_ci.type(ci.type, .vcov_type)

    if (ci.type == "wald") {
      # Process df
      df <- df %or% .attr(object, ".df")

      arg::arg_number(df)
      arg::arg_gt(df, 0)

      stat <- if (is.finite(df)) "t" else "z"

      res_names <- c("term", "estimate",
                     "std.error",
                     stat[!anyNA(null)],
                     "p.value"[!anyNA(null)],
                     "conf.low"[conf_level > 0],
                     "conf.high"[conf_level > 0])

      res <- make_df(res_names, length(.s))

      se <- sqrt(diag(.proj_vcov))
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
        # Compute test statistics when a null value is specified
        res[[stat]][!zeros] <- (.proj_coefficients[!zeros] - null) / se[!zeros]

        res$p.value[!zeros] <- switch(stat,
                                      "z" = 2 * pnorm(-abs(res[[stat]][!zeros])),
                                      "t" = 2 * pt(-abs(res[[stat]][!zeros]), df))
      }
      else if (.vcov_type == "bootstrap") {
        s <- summary(.proj_boot,
                     conf = 0,
                     index = .s,
                     ci.type = ci.type,
                     p.value = TRUE,
                     null = null,
                     simultaneous = FALSE)

        res$p.value <- s[, fncol(s)]
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
    if (conf_level > 0 && !is_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(res, "crit") <- t_crit
    }

    attr(res, "ci.type") <- ci.type

    if (!anyNA(null)) {
      attr(res, "null") <- null
    }
  }

  attr(res, ".curve_type") <- get_curve_type(object)
  attr(res, ".shift") <- .attr(object, ".shift")
  attr(res, ".contrast") <- .contrast
  attr(res, ".by_grid") <- .by_grid
  attr(res, ".outcome") <- .attr(object, ".outcome")
  attr(res, ".vcov_type") <- .attr(object, ".vcov_type")

  attr(res, ".proj_coefficients") <- .proj_coefficients
  attr(res, ".proj_vcov") <- .proj_vcov
  attr(res, ".proj_coefficient_names") <- .proj_coefficient_names

  class(res) <- unique(c("summary.curve_projection", class(res)))

  res
}

#' @exportS3Method print summary.curve_projection
print.summary.curve_projection <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  arg::arg_whole_number(digits)

  .reference <- .attr(x, ".reference")
  .contrast <- .attr(x, ".contrast")

  main <- paste(c(.curve_type_label(x, title = TRUE),
                  if (!is_null(.reference)) "Reference",
                  if (!is_null(.contrast)) "Contrast",
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

  # No covariance matrix when vcov = "none" in the original call to adrf()
  if (is_null(v)) {
    return(NULL)
  }

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

  proj_data <- qDF(.values) |>
    setNames(.treat) |>
    model.frame(formula = formula(object), na.action = stats::na.pass)

  mm <- .attr(proj_data, "terms") |>
    model.matrix(data = proj_data)

  if (!all(is.finite(mm))) {
    arg::err("evaluating the projection model at the values of {.var {(.treat)}} produced non-finite values, which is not allowed")
  }

  n_by <- get_n_by(.contrast, .by_grid)

  block_diag(mm, n_by) |>
    setColnames(.get_curve_est_labels(object))
}

#' @exportS3Method stats::formula curve_projection
formula.curve_projection <- function(x, ...) {
  .attr(x, ".proj_formula") |> as.formula()
}

#' @exportS3Method stats::anova curve_projection
#' @rdname curve_projection
anova.curve_projection <- function(object, object2, df = NULL, ...) {
  arg::arg_supplied(object)
  arg::arg_is(object, "curve_projection")

  arg::arg_supplied(object2)
  arg::arg_is(object2, "curve_projection")

  if (any_apply(c(".values", ".treat", ".vcov_type", ".curve_type", ".response",
                  ".by_grid", ".contrast", ".outcome", ".reference", ".family", ".df",
                  ".shift"),
                function(a) !identical(.attr(object, a), .attr(object2, a)))) {
    arg::err("{.arg object} and {.arg object2} must be the outputs of {.fun curve_projection} applied to the same {.cls effect_curve} object")
  }

  if (!isTRUE(all.equal(.attr(object, ".proj_transformed_est"),
                        .attr(object2, ".proj_transformed_est")))) {
    arg::err("{.arg object} and {.arg object2} must use the same {.arg transform}")
  }

  if (!identical(.attr(object, ".proj_integrate"), .attr(object2, ".proj_integrate"))) {
    arg::err("{.arg object} and {.arg object2} must use the same {.arg integrate}")
  }

  if (is_null(stats::vcov(object))) {
    arg::err("{.fun anova} cannot be used when {.code vcov = \"none\"} in the original call to {.fun {(.curve_fun(object))}}")
  }

  # Process df
  df <- df %or% .attr(object, ".df")

  arg::arg_number(df)
  arg::arg_gt(df, 0)

  test <- if (is.finite(df)) "F" else "Chisq"

  tolerance <- ...get("tolerance", 1e-7)

  arg::arg_number(tolerance)
  arg::arg_gt(tolerance, 0)

  b1 <- coef(object, complete = FALSE)
  b2 <- coef(object2, complete = FALSE)

  if (length(b2) >= length(b1)) {
    arg::err("{.arg object2} does not appear to be nested within {.arg object}")
  }

  Z1 <- .lm.fit(x = model.matrix(object2)[, names(b2), drop = FALSE],
                y = model.matrix(object)[, names(b1), drop = FALSE])[["residuals"]]

  Z1_svd <- svd(Z1)
  keep <- Z1_svd[["d"]] > tolerance
  .q <- sum(keep)

  if (.q > length(b1) - length(b2)) {
    arg::err("{.arg object2} does not appear to be nested within {.arg object}")
  }

  L <- t(Z1_svd[["v"]][, keep, drop = FALSE])

  V <- stats::vcov(object)

  value.hyp <- L %*% b1
  vcov.hyp <- quad_mult(L, V)

  # The effect curve varies only in the span of the treatment terms of the
  # outcome model, so a projection with more terms than that has coefficients
  # whose covariance matrix is singular
  ev <- eigen(vcov.hyp, symmetric = TRUE, only.values = TRUE)[["values"]]

  if (min(ev) <= sqrt(.Machine$double.eps) * max(ev, 0)) {
    arg::err(c("The extra terms in {.arg object} cannot be tested because the covariance matrix of their coefficients is singular.",
               i = "This happens when the projection in {.arg object} has more terms than the effect curve can identify, e.g., more than the number of terms involving the treatment in the outcome model."))
  }

  SSH <- drop(crossprod(value.hyp, solve(vcov.hyp, value.hyp)))

  width <- min(55L, getOption("width", 55L))

  .title <- "Wald Test"

  # Each formula is wrapped as a whole and the lines joined, so a long formula
  # stays in a single heading element
  .topnote <- c(sprintf("Model 1: %s", deparse1(formula(object))) |>
                  strwrap(width, exdent = 9L),
                sprintf("Model 2: %s", deparse1(formula(object2))) |>
                  strwrap(width, exdent = 9L)) |>
    paste(collapse = "\n") |>
    paste0("\n")

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

  cli::cat_line(out)

  .print_estimate_table(x, digits = digits, topn = Inf, bar = FALSE)

  .df <- .attr(x, ".df")

  c(txtbar(width),
    .it(sprintf("Inference: %s%s",
                .get_vcov_type_name(attr(x, ".vcov_type")),
                if (is.finite(.df)) sprintf(" (df = %s)", round(.df, 1L)) else "")) |>
      cli::ansi_strwrap(width, exdent = 2L)) |>
    cli::cat_line()

  invisible(x)
}

# The derivative of the projection coefficients with respect to the curve
# estimates `y`. The projection solves J'W(y - mu) = 0 with J = diag(mu') X,
# so the derivative is H^-1 J'W with H = J'WJ - X' diag(w r mu'') X, where
# r = y - mu are the residuals and mu'' is the second derivative of the inverse
# transformation; the second term is zero for the identity, and is not zero in
# general because a projection is misspecified by design. Columns are
# equilibrated before solving, so raw powers of a treatment on a large scale
# stay well conditioned.
.proj_influence <- function(x, y, fit, w, transform) {
  J <- fit$jacobian

  eta <- drop(x %*% fit$coefficients)
  r <- y - fit$fitted.values

  h <- 1e-4 * pmax(1, abs(eta))
  mu2 <- (transform$d_inv_transform(eta + h) - transform$d_inv_transform(eta - h)) / (2 * h)

  d <- sqrt(colSums(J^2 * w))
  d[d == 0] <- 1

  Js <- J / rep(d, each = fnrow(J))
  xs <- x / rep(d, each = fnrow(x))

  H <- crossprod(Js, Js * w) - crossprod(xs, xs * (w * r * mu2))

  qr.solve(H, t(Js * w)) / d
}

# Coefficients of the projection fit to each row of `Y` (the bootstrap or
# posterior draws of the effect curve), one row per draw. With the identity
# transformation the projection is linear in the draws, so every draw is fit at
# once from a single QR decomposition rather than one fit per draw.
.proj_coef_draws <- function(x, Y, w, transform, start) {
  if (.is_identity_transform(transform)) {
    wsqt <- sqrt(w)

    return(t(qr.coef(qr(wsqt * x), wsqt * t(Y))))
  }

  # With one coefficient, dapply() returns a vector, so the result is coerced
  # to a matrix with one column per coefficient
  Y |>
    dapply(MARGIN = 1L, return = "matrix",
           FUN = function(y) {
             .nls_fit(x = x, y = y, w = w, transform = transform,
                      start = start)$coefficients
           }) |>
    matrix(ncol = length(start))
}

.is_identity_transform <- function(transform) {
  v <- c(-2, .3, 5)

  tv <- try(suppressWarnings(transform$transform(v)), silent = TRUE)

  !null_or_error(tv) && isTRUE(all.equal(tv, v))
}

.nls_fit <- function(x, y, w = NULL, offset = 0, transform, start = NULL,
                     maxit = 50, tol = 1e-8, jac = FALSE) {

  f <- transform$transform
  finv <- transform$inv_transform
  finv_prime <- transform$d_inv_transform

  n <- fnrow(x)

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

# Checks that the coefficients of the projection model can be estimated from
# the values of the treatment at which the effect curve was computed. Columns
# are scaled to unit length first so that raw powers of a treatment on a large
# scale are not mistaken for linearly dependent terms.
.check_proj_mm <- function(mm) {
  n <- fnrow(mm)

  if (fncol(mm) > n) {
    arg::err(c("The projection model has {fncol(mm)} coefficients, but the effect curve was computed at only {n} grid points.",
               i = "Use a simpler projection model or increase {.arg n} in the call to {.fun adrf} or {.fun shift_curve}."))
  }

  d <- sqrt(colSums(mm^2))
  d[d == 0] <- 1

  if (qr(mm / rep(d, each = n))$rank < fncol(mm)) {
    arg::err("the terms in the projection model are linearly dependent at the grid points at which the effect curve was computed")
  }

  invisible(mm)
}
