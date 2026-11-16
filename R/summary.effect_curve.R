#' Test omnibus hypotheses about an effect curve
#'
#' `summary()` tests an omnibus hypothesis about an effect curve. For example, it can be used to test that the ADRF is flat, that the contrast between two ADRFs is 0 everywhere, or that the AMEF is 0 everywhere.
#'
#' @param object an object of class [`effect_curve`]; the result of a call to [adrf()] or a function that modifies it.
#' @param hypothesis the hypothesis to be tested. Allowable options include `"flat"`, `"linear"`, `"quadratic"`, `"cubic"`, a one-sided formula corresponding to a projection model, or a single number (e.g., 0). See Details. The default is `"flat"` for ADRFs and 0 otherwise.
#' @param method string; the method used to compute the p-value of the test. Allowable options include `"sim"` for simulation-based, `"imhof"` for Imhof's approximation, `"davies"` for Davies's approximation, `"liu"` for Liu's approximation, `"satterthwaite"` for Satterthwaite's approximation, and `"saddlepoint"` for a saddlepoint approximation. Default is `"imhof"` when the \pkg{CompQuadForm} package is installed and `"saddlepoint"` otherwise. See Details.
#' @param subset an optional logical expression indicating the subset of the subgroups for which to perform tests. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param transform whether to perform the test on the transformed estimates. Allowable options include `TRUE`, `FALSE`, a family object (e.g., `binomial()`), a link object (e.g., the output of [make.link()]), or one of the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`. Ignored unless `object` is an ADRF object. See Details.
#' @param df the "denominator" degrees of freedom to use for the tests. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise, or, for models fit to multiply imputed data, degrees of freedom that account for the number of imputations (see [adrf()]).
#' @param nsim when `method` is `"sim"`, the number of iterations used to simulate the p-values. Higher numbers give more accurate p-values subject to less Monte Carlo error but are slower and require more memory. Default is 1,000,000.
#' @param \dots when `method` is `"imhof"` or `"davies"`, further arguments passed to \pkgfun{CompQuadForm}{imhof} or \pkgfun{CompQuadForm}{davies}, respectively.
#' @param x a `summary.effect_curve` object.
#' @param digits numeric; the number of digits to print.
#'
#' @returns
#' An object of class `"summary.effect_curve"`, which is a data.frame with a column for the p-value and, for stratified effect curves or contrasts thereof, additional columns identifying the subset to which the p-value refers.
#'
#' @details
#' `summary()` performs an omnibus test for an effect curve. The hypothesis tested is determined by the argument to `"hypothesis"`. When supplied as a single number, `summary()` tests whether all values on the effect curve are equal to that number. When supplied as a one-sided formula, `summary()` tests whether the projection of the effect curve model onto the model represented in the formula is sufficient to describe the effect curve. The test itself tests whether the residuals around the projection are all equal to 0, incorporating the uncertainty in estimating the effect curve. See [curve_projection()] for more information on how the curve projection and the uncertainty in the residuals are computed.
#'
#'  `"flat"` tests whether all values on the curve are equal to each other (i.e., whether the curve is flat), without specifying what value they are equal to. This is equivalent to testing whether the variance around the mean estimate is different from 0 or whether an intercept-only projection model is sufficient. `"linear"` tests whether the curve is linear, i.e., whether the residuals around linear projection of the curve are all 0. `"quadratic"` and `"cubic"` test whether the curve is quadratic or cubic, respectively, using the same method.
#'
#'  Rejecting the null hypothesis means that the curve is more complicated than the specified model. For example, rejecting the null hypothesis that the curve is linear implies that the curve is nonlinear (and, therefore, not flat either).
#'
#' The test involves computing a test statistic, specifying its distribution under the null hypothesis, and computing the p-value using the complement of the cumulative distribution function of the distribution evaluated at the test statistic value. The test statistic depends on `"hypothesis"`. For `hypothesis` equal to a constant, say, \eqn{c}, the test statistic is
#'  \deqn{T^* = \int_\mathcal{A} (\theta(a) - c)^2 \ da}
#' Otherwise, the test statistic is
#' \deqn{T^* = \int_\mathcal{A} (\theta(a) - \hat{\theta}_0(a))^2 \ da}
#' where \eqn{\hat{\theta}_0} is the projection of \eqn{\theta(a)} onto the null subspace specified by `hypothesis`.
#'
#' Each of these can be approximated as a quadratic form, \eqn{T=\mathbf{\theta}'\mathbf{W}\mathbf{\theta}} where \eqn{\mathbf{\theta}} is a vector of estimates at evaluation points along the curve and \eqn{\mathbf{W}} is a diagonal matrix of weights implementing a trapezoidal approximation to the integral. The null hypothesis is that \eqn{T=0}, which approximates the null hypothesis that \eqn{T^*=0}. Each of the allowable options to `method` corresponds to a method of approximating the distribution of \eqn{T} under the null hypothesis:
#' * `"sim"` simulates the null distribution by simulating from a multivariate normal distribution under the null hypothesis and computing the test statistic in each simulation. The p-value is the proportion of simulated estimates greater than the observed test statistic. When `df` is not `Inf`, the simulation is done from a multivariate t-distribution.
#' * `"imhof"`, `"davies"`, and `"liu"` assume the test statistic follows a generalized \eqn{\chi^2}-distribution and approximate its CDF numerically. `"imhof"` tends to be the most accurate and is recommended. Because Imhof's and Davies's methods control the absolute rather than the relative error of the computed probability, p-values below .0001, and any that these methods fail to compute, are computed using the saddlepoint approximation instead. `"liu"` can be inaccurate when `df` is not `Inf`. These methods require the \CRANpkg{CompQuadForm} package to be installed.
#' * `"satterthwaite"` also assumes the test statistic follows a generalized \eqn{\chi^2}-distribution, but this distribution is approximated using a scaled F-distribution with the same first two moments.
#' * `"saddlepoint"` also assumes the test statistic follows a generalized \eqn{\chi^2}-distribution, and this distribution is approximated using the saddlepoint method described by Kuonen (1999), which remains accurate far into the tails of the distribution. Near the center of the distribution, where the saddlepoint approximation is not defined, the Satterthwaite approximation is used instead.
#'
#' In general, we recommend using `method = "imhof"`, though this requires \pkg{CompQuadForm} to be installed (and is the default when it is). `method = "saddlepoint"` is also accurate, is very fast, and doesn't require any other packages. When using `"sim"`, increasing `nsim` improves the accuracy of the p-value by reducing Monte Carlo error. The default value of `1e6` ensures that the simulated p-value is within .0005 of the true p-value with at least 97.8% confidence for p-values less than .05. For greater precision at higher p-value thresholds, `nsim` should be increased. The other methods exist primarily for comparison purposes. `"satterthwaite"` is particularly fast and doesn't require any other packages, but it can be somewhat inaccurate.
#'
#' ## Transform
#'
#' When the effect curve is an ADRF and the outcomes are bounded (e.g., a probability between 0 and 1), the `transform` argument can be specified, which changes the details of the tests.
#'
#' The tests above assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. How `transform` is used depends on `hypothesis`. When `hypothesis` is a number, the estimates and that number are transformed to be on an unbounded scale, and the test is performed on the transformed estimates, whose distribution is computed using the delta method (or from the transformed bootstrap or posterior draws). When `hypothesis` is a string or formula, the estimates are not transformed; instead, the projection is fit with the inverse of the transformation as its inverse link, as described in [curve_projection()], and the residuals of the estimates around it are tested. By default, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simply a transformation of the estimated points along the curve already computed. A value of `hypothesis` outside the domain of the transformation (e.g., `hypothesis = 0` with a logit transformation) gives an error; set `transform = FALSE` to perform the test on the original scale. A transformation supplied to `transform` that cannot be applied to every estimate or draw also gives an error, and when `transform = TRUE`, the estimates are left untransformed with a warning in that case.
#'
#' Tests on the transformed and untransformed ADRFs correspond to different hypotheses; the difference is not simply a matter of the appropriate distribution of the statistic. For example, for a binary outcome model with a logistic transformation, testing that the transformed ADRF is linear corresponds to testing whether the ADRF is a sigmoid function, whereas testing that the untransformed ADRF is linear corresponds to testing whether the ADRF is a straight line. These choices correspond to how the projection of the ADRF is formed; see [curve_projection()] for details. See Examples for an example with a binary outcome. In the example, there is evidence to reject that the transformed ADRF is linear for one of the groups, indicating that a sigmoid function is not sufficient for describing the ADRF, but there is indeterminate evidence (at the .05 level) to reject that the untransformed ADRF is linear for either group, indicating that a linear function is sufficient for describing the ADRF (at least in the treatment range examined).
#'
#' @references
#' Kuonen, D. (1999). Saddlepoint approximations for distributions of quadratic forms in normal variables. *Biometrika*, 86(4), 929–935. \doi{10.1093/biomet/86.4.929}
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [summary.curve_est()] for performing inference on individual points on an effect curve
#' * [plot.effect_curve()] for plotting the effect curve
#' * [curve_projection()] for projecting a simpler model onto an effect curve
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- glm(Block >= 12 ~ poly(logBLL, 3) *
#'              Male * (Age + Race + PIR + NICU +
#'                        Smoke_Pregnant),
#'            data = nhanes3lead,
#'            family = binomial)
#'
#' # ADRFs of logBLL on P(Block >= 12) within
#' # groups defined by `Male`
#' adrf1 <- adrf(fit, treat = "logBLL",
#'               by = ~Male)
#'
#' adrf1
#'
#' # Test if ADRFs are flat
#' summary(adrf1)
#'
#' # Test if logit-transformed ADRFs are linear
#' # (i.e., if ADRFs have sigmoid shape)
#' summary(adrf1, hypothesis = "linear")
#' # summary(adrf1, hypothesis = ~logBLL) # same test
#'
#' proj1 <- curve_projection(adrf1, "linear")
#'
#' plot(adrf1, proj = proj1, conf_level = 0)
#'
#' # Test if un-transformed ADRFs are linear
#' summary(adrf1, hypothesis = "linear",
#'         transform = FALSE)
#'
#' proj2 <- curve_projection(adrf1, "linear",
#'                           transform = FALSE)
#'
#' plot(adrf1, proj = proj2, conf_level = 0)
#'
#' ## Test if ADRFs differ from each other by
#' ## testing if the ADRF contrast is 0
#' curve_contrast(adrf1) |>
#'   summary()

#' @exportS3Method summary effect_curve
summary.effect_curve <- function(object, hypothesis, method, subset = NULL,
                                 transform = TRUE, df = NULL, nsim = 1e6, ...) {

  .vcov_type <- .attr(object, ".vcov_type")

  if (.vcov_type == "none") {
    arg::err('{.fun summary} cannot be used when {.code vcov = "none"} in the original call to {.fun adrf}')
  }

  if (missing(hypothesis) || is_null(hypothesis)) {
    hypothesis <- if (.is_pure_adrf(object)) "flat" else 0
  }
  else if (rlang::is_string(hypothesis)) {
    hypothesis <- arg::match_arg(hypothesis, c("flat", "linear", "quadratic", "cubic"))
  }
  else if (rlang::is_formula(hypothesis)) {
    if (!rlang::is_formula(hypothesis, lhs = FALSE)) {
      arg::err("if {.arg hypothesis} is a formula, it must be a one-sided formula with the projection model on the right-hand side")
    }

    .treat <- .attr(object, ".treat")

    #Check that only treat is in formula
    vars_in_formula <- get_varnames(hypothesis)

    if (!all(get_varnames(hypothesis) %in% .treat)) {
      arg::err("only the treatment variable {.var {(.treat)}} is allowed to appear in {.arg hypothesis} when supplied as a formula")
    }
  }
  else if (!is_number(hypothesis)) {
    arg::err("{.arg hypothesis} must be a string, formula, or a number")
  }

  if (missing(method) || is_null(method)) {
    method <- {
      if (rlang::is_installed("CompQuadForm")) "imhof"
      else "saddlepoint"
    }
  }

  method <- arg::match_arg(method, c("sim", "imhof", "davies", "liu",
                                     "satterthwaite", "saddlepoint"))

  if (method %in% c("imhof", "davies", "liu")) {
    rlang::check_installed("CompQuadForm")
  }

  if (method == "sim") {
    arg::arg_count(nsim)
    arg::arg_gt(nsim, 100)
  }
  else {
    nsim <- NULL
  }

  df <- df %or% .attr(object, ".df")

  arg::arg_number(df)
  arg::arg_gt(df, 0)

  .est0 <- .attr(object, ".est")
  .vcov0 <- .attr(object, ".vcov")
  .by_grid <- .attr(object, ".by_grid")
  .contrast <- .attr(object, ".contrast")
  .draws <- .attr(object, ".draws")
  .boot <- .attr(object, ".boot")
  .values <- .attr(object, ".values")

  # transform
  transform_list <- process_transform(transform, object, .est0,
                                      null = if (is.numeric(hypothesis)) hypothesis,
                                      null_arg = "hypothesis")

  transformed <- !all(check_if_zero(transform_list$transform(.est0) - .est0))

  # subset; process_subset_by_grid() gives an error when subset is supplied
  # without by or for a contrast
  .subset <- process_subset_by_grid(substitute(subset),
                                    .by_grid = .by_grid,
                                    .contrast = .contrast)

  if (is_null(.by_grid) || !is_null(.contrast)) {
    .s <- seq_along(.est0)
  }
  else {
    .s <- which(rep(.subset, each = length(.est0) / fnrow(.by_grid)))

    .by_grid <- ss(.by_grid, .subset)

    .est0 <- .est0[.s]
    .vcov0 <- ss(.vcov0, .s, .s)

    if (!is_null(.boot)) {
      .boot[["t"]] <- ss(.boot[["t"]], j = .s)
      .boot[["t0"]] <- .boot[["t0"]][.s]
    }

    if (!is_null(.draws)) {
      .draws <- ss(.draws, j = .s)
    }
  }

  .values <- .attr(object, ".values")

  n <- length(.values)

  by_list <- {
    if (is_null(.by_grid) && is_null(.contrast)) list(seq_len(n))
    else if (!is_null(.contrast)) gsplit(seq_along(.est0), rep(seq_along(.contrast), each = n))
    else gsplit(seq_along(.est0), rep(seq_row(.by_grid), each = n))
  }

  if (rlang::is_string(hypothesis)) {
    n_coef <- switch(hypothesis, flat = 1L, linear = 2L, quadratic = 3L, cubic = 4L)

    if (n_coef > n) {
      .check_proj_mm(matrix(1, nrow = n, ncol = n_coef))
    }

    mm <- switch(hypothesis,
                 flat = matrix(1, nrow = n, ncol = 1L),
                 linear = cbind(1, .values),
                 quadratic = cbind(1, poly(.values, 2)),
                 cubic = cbind(1, poly(.values, 3)))
  }
  else if (rlang::is_formula(hypothesis)) {
    proj_data <- qDF(.values) |>
      setNames(.treat) |>
      model.frame(formula = hypothesis, na.action = stats::na.pass)

    mm <- .attr(proj_data, "terms") |>
      model.matrix(data = proj_data)

    if (!all(is.finite(mm))) {
      arg::err("evaluating the projection model at the values of {.var {(.treat)}} produced non-finite values, which is not allowed")
    }

    .check_proj_mm(mm)
  }

  #Weights for trapezoidal Riemann sum
  w <- get_trapezoidal_w(.values)

  wsq <- sqrt(w)

  p <- vapply(by_list, function(.b) {

    if (!is_null(.boot)) {
      resid_boot <- .boot

      if (is.numeric(hypothesis)) {
        resid_boot[["t0"]] <- est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(hypothesis)
        resid_boot[["t"]] <- transform_list$transform(ss(resid_boot[["t"]], j = .b)) - transform_list$transform(hypothesis)

        vcov.b <- stats::vcov(resid_boot)
      }
      else {
        fit <- .nls_fit(x = mm, y = .est0[.b], w = w,
                        transform = transform_list)

        resid_boot[["t0"]] <- est.b <- .est0[.b] - fit$fitted.values
        resid_boot[["t"]] <- .proj_resid_draws(mm, ss(resid_boot[["t"]], j = .b), w,
                                               transform_list, fit$coefficients)

        vcov.b <- stats::vcov(resid_boot)
      }
    }
    else if (!is_null(.draws)) {
      if (is.numeric(hypothesis)) {
        est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(hypothesis)

        resid_draws <- transform_list$transform(ss(.draws, j = .b)) - transform_list$transform(hypothesis)

        vcov.b <- stats::cov(resid_draws)
      }
      else {
        fit <- .nls_fit(x = mm, y = .est0[.b], w = w,
                        transform = transform_list)

        resid_draws <- .proj_resid_draws(mm, ss(.draws, j = .b), w,
                                         transform_list, fit$coefficients)

        est.b <- .est0[.b] - fit$fitted.values

        vcov.b <- stats::cov(resid_draws)
      }
    }
    else {
      if (is.numeric(hypothesis)) {
        est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(hypothesis)
        vcov.b <- quad_mult(diag(transform_list$d_transform(.est0[.b]), nrow = length(.b)),
                            ss(.vcov0, .b, .b))
      }
      else {
        fit <- .nls_fit(x = mm, y = .est0[.b], w = w,
                        transform = transform_list,
                        jac = TRUE, tol = 0, maxit = 20)

        H <- fit$jacobian %*% .proj_influence(mm, .est0[.b], fit, w, transform_list)

        est.b <- .est0[.b] - fit$fitted.values

        vcov.b <- quad_mult(diag(1, n) - H, ss(.vcov0, .b, .b))

        ## Uncomment below to perform test on difference between transformed
        ## estimates and transformed fitted values from projection
        # est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(fit$fitted.values)
        #
        # vcov.b <- quad_mult(diag(transform_list$d_transform(.est0[.b])) - transform_list$d_transform(fit$fitted.values) * H,
        #                     ss(.vcov0, .b, .b))
      }
    }

    # The covariance of the estimates themselves sets the scale against which
    # numerical tolerances are judged, so results do not depend on the units of
    # the outcome
    ref.b <- {
      if (is.numeric(hypothesis)) vcov.b
      else if (!is_null(.boot)) stats::cov(ss(.boot[["t"]], j = .b))
      else if (!is_null(.draws)) stats::cov(ss(.draws, j = .b))
      else ss(.vcov0, .b, .b)
    }

    stat <- sum((est.b * wsq)^2)

    Wsq <- diag(wsq, nrow = length(wsq))

    ev <- quad_mult(Wsq, vcov.b) |>
      eigen(symmetric = TRUE, only.values = TRUE) |>
      getElement("values")

    scale <- quad_mult(Wsq, ref.b) |>
      eigen(symmetric = TRUE, only.values = TRUE) |>
      getElement("values") |>
      max()

    .omnibus_p_value(stat, ev, scale, method = method,
                     df = df, nsim = nsim, ...)
  }, numeric(1L))

  res_names <- "p.value"

  out <- make_df(res_names, length(p)) |>
    ftransform(p.value = squish(p, lo = 1e-12, hi = 1)) |>
    add_est_labels(.contrast, .by_grid)

  attr(out, "method") <- method
  attr(out, "hypothesis") <- hypothesis

  attr(out, "transformed") <- transformed
  attr(out, "nsim") <- nsim
  attr(out, "df") <- df

  attr(out, ".treat") <- .attr(object, ".treat")
  attr(out, ".values") <- .values
  attr(out, ".curve_type") <- get_curve_type(object)
  attr(out, ".contrast") <- .contrast
  attr(out, ".by_grid") <- .by_grid
  attr(out, ".reference") <- .attr(object, ".reference")

  class(out) <- c("summary.effect_curve", class(out))

  out
}

#' @exportS3Method print summary.effect_curve
#' @rdname summary.effect_curve
print.summary.effect_curve <- function(x, digits = max(4L, getOption("digits") - 3L), ...) {
  .treat <- .attr(x, ".treat")
  .values <- .attr(x, ".values")
  .contrast <- .attr(x, ".contrast")
  .reference <- .attr(x, ".reference")
  nsim <- .attr(x, "nsim")
  transformed <- .attr(x, "transformed")
  method <- .attr(x, "method")
  hypothesis <- .attr(x, "hypothesis")

  width <- min(55L, getOption("width", 55L))

  qoi <- paste(c(if (isTRUE(transformed) && is.numeric(hypothesis)) "transformed",
                 get_curve_type(x),
                 if (!is_null(.contrast)) "contrast",
                 if (!is_null(.reference)) sprintf("difference from reference (%s = %s)",
                                                      .treat, round(.reference, 4L))),
               collapse = " ")

  out <- c(center_just("Omnibus Curve Test", wrt = space(width)),
           txtbar(width),
           sprintf("H\u2080: %s is %s%s for values of %s between %s and %s",
                   qoi,
                   if (rlang::is_formula(hypothesis)) deparse1(hypothesis)
                   else hypothesis,
                   if (transformed && !is.numeric(hypothesis)) " (transformed)"
                   else "",
                   .treat,
                   round(min(.values), digits),
                   round(max(.values), digits)) |>
             cli::ansi_strwrap(width, exdent = 4))

  cli::cat_line(c(out, ""))

  .print_estimate_table(x, digits = digits, topn = Inf, bar = FALSE)

  method <- switch(method,
                   "satterthwaite" = "the Satterthwaite approximation",
                   "imhof" = "the Imhof approximation (the saddlepoint approximation for p-values below .0001)",
                   "davies" = "the Davies approximation (the saddlepoint approximation for p-values below .0001)",
                   "liu" = "the Liu approximation",
                   "sim" = sprintf("a simulation approximation with %s replicates",
                                   prettyNum(nsim, big.mark = ",", scientific = FALSE)),
                   "saddlepoint" = "Kuonen's saddlepoint approximation")

  c(txtbar(width),
    cli::format_inline("{.emph Computed using {method}}") |>
      cli::ansi_strwrap(width, exdent = 2L)) |>
    cli::cat_line()

  invisible(x)
}

# Upper-tail probability of a quadratic form in normal variables: P(Q > stat),
# where Q = sum_k ev_k * chi^2_1 and ev are the (positive) eigenvalues of its
# covariance. With finite `df`, the covariance is estimated, and the
# probability is P(Q > stat * chi^2_df / df).
.quad_form_tail <- function(stat, ev, method = "imhof", df = Inf, nsim = 1e6, ...) {
  if (method == "sim") {
    return(.sim_p_value(stat, ev, n = nsim, df = df))
  }

  if (method == "satterthwaite") {
    # Moment matching
    mu <- sum(ev)
    s2 <- sum(ev^2)

    return(pf(stat / mu, df1 = mu^2 / s2, df2 = df, lower.tail = FALSE))
  }

  # Rescale to make stat = 1 (improves numerical performance)
  ev <- ev / stat

  if (is.finite(df)) {
    .q <- 0
    .lambda <- c(ev, -1 / df)
    .h <- c(alloc(1.0, length(ev)), df)
  }
  else {
    .q <- 1
    .lambda <- ev
    .h <- alloc(1.0, length(ev))
  }

  if (method == "saddlepoint") {
    p <- .saddlepoint_tail(.q, .lambda, .h)

    # The approximation is singular near the mean of the distribution, where
    # the Satterthwaite approximation is used instead, as in survey
    if (is.na(p)) {
      p <- pf(1 / sum(ev), df1 = sum(ev)^2 / sum(ev^2), df2 = df, lower.tail = FALSE)
    }

    return(p)
  }

  out <- suppressWarnings({
    switch(method,
           "imhof" = CompQuadForm::imhof(q = .q, lambda = .lambda,
                                         h = .h,
                                         epsabs = ...get("epsabs", 1e-6),
                                         epsrel = ...get("epsrel", 1e-6),
                                         limit = ...get("limit", 1e4)),
           "davies" = CompQuadForm::davies(q = .q, lambda = .lambda,
                                           h = .h,
                                           lim = ...get("lim", 1e4),
                                           acc = ...get("acc", 0.0001)),
           "liu" = list(Qq = CompQuadForm::liu(q = .q, lambda = .lambda,
                                               h = .h)))
  })

  p <- out[["Qq"]]

  # Imhof's and Davies's methods control the absolute error of the
  # probability, so small p-values can be inaccurate, negative, or flagged as
  # unreliable; the saddlepoint approximation, whose error is relative, is used
  # for them instead
  if (method %in% c("imhof", "davies") &&
      (!is.finite(p) || p < 1e-4 || p > 1 ||
       (!is_null(out[["ifault"]]) && out[["ifault"]] != 0))) {
    p_sp <- .saddlepoint_tail(.q, .lambda, .h)

    if (!is.na(p_sp)) {
      p <- p_sp
    }
  }

  p
}

# P-value of the omnibus test given the test statistic and the eigenvalues of
# the weighted covariance of the tested vector. Tolerances are relative to
# `scale`, the largest eigenvalue of the weighted covariance of the estimates,
# so the result does not depend on the units of the outcome. The result is
# restricted to [0, 1].
.omnibus_p_value <- function(stat, ev, scale, method, df = Inf, nsim = 1e6,
                             tol = 1e-10, ...) {
  scale <- max(scale, ev, 0)

  if (!is.finite(stat) || !is.finite(scale)) {
    return(NaN)
  }

  ev <- ev[ev > tol * scale]

  # When the tested quantity has no sampling variance, the hypothesis either
  # holds to numerical precision or does not hold at all (e.g., testing the
  # linearity of a curve that is linear by construction)
  if (is_null(ev)) {
    if (stat <= tol * max(scale, .Machine$double.xmin)) {
      return(1)
    }

    return(0)
  }

  if (stat <= 0) {
    return(1)
  }

  .quad_form_tail(stat, ev, method = method, df = df, nsim = nsim, ...) |>
    squish(lo = 0, hi = 1)
}

# Saddlepoint approximation to P(sum_j lambda_j chi^2_{h_j} > q) (Kuonen, 1999),
# using the Barndorff-Nielsen form, as in survey:::saddle(). Negative lambdas
# are allowed, which is how finite df are handled. Returns NA near the mean of
# the distribution, where the approximation is singular.
.saddlepoint_tail <- function(q, lambda, h = rep(1, length(lambda))) {
  d <- max(abs(lambda))
  lambda <- lambda / d
  q <- q / d

  K <- function(zeta) -sum(h * log1p(-2 * zeta * lambda)) / 2
  K1 <- function(zeta) sum(h * lambda / (1 - 2 * zeta * lambda))
  K2 <- function(zeta) 2 * sum(h * lambda^2 / (1 - 2 * zeta * lambda)^2)

  # The saddlepoint lies where 1 - 2 * zeta * lambda > 0 for every lambda
  upper <- {
    if (any(lambda > 0)) .99999 / (2 * max(lambda))
    else 1e3
  }

  lower <- {
    if (any(lambda < 0)) .99999 / (2 * min(lambda))
    else -sum(h) / (2 * q)
  }

  root <- try(uniroot(function(z) K1(z) - q, lower = lower, upper = upper,
                      tol = 1e-12)$root, silent = TRUE)

  if (null_or_error(root) || abs(root) < 1e-4) {
    return(NA_real_)
  }

  w <- sign(root) * sqrt(2 * (root * q - K(root)))
  v <- root * sqrt(K2(root))

  pnorm(w + log(v / w) / w, lower.tail = FALSE)
}

# Residuals of each draw (row of `Y`) around its own projection
.proj_resid_draws <- function(x, Y, w, transform, start) {
  coefs <- .proj_coef_draws(x, Y, w, transform, start)

  (Y - transform$inv_transform(tcrossprod(coefs, x))) |>
    unname()
}
