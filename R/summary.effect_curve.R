#' Test omnibus hypotheses about an effect curves
#'
#' `summary()` tests an omnibus hypothesis about an effect curve. For example, it can be used to test that the ADRF is flat, that the contrast between two ADRFs is 0 everywhere, or that the AMEF is 0 everywhere.
#'
#' @param object an object of class [`effect_curve`]; the result of a call to [adrf()] or a function that modifies it.
#' @param hypothesis the hypothesis to be tested. Allowable options include `"flat"` (the default), `"linear"`, `"quadratic"`, `"cubic"`, a one-sided formula corresponding to a projection model, or a single number (e.g., 0). See Details. The default is `"flat"` for ADRFs and 0 otherwise.
#' @param method string; the method used to compute the p-value of the test. Allowable options include `"sim"` for simulation-based, `"imhof"` for Imhof's approximation, `"davies"` for Davies's approximation, `"liu"` for Liu's approximation, `"satterthwaite"` for Satterthwaite's approximation, and `"saddlepoint"` for a saddlepoint approximation. Default is `"imhof"` when the \pkg{CompQuadForm} package is installed, otherwise `"saddlepoint"` when the \pkg{survey} package is installed, and `"sim"` otherwise. See Details.
#' @param subset an optional logical expression indicating the subset of the subgroups for which to perform tests. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param transform whether to perform the test on the transformed estimates. Allowable options include `TRUE`, `FALSE`, or a function specifying a transformation. Ignored unless `object` is an ADRF object. See Details.
#' @param df the "denominator" degrees of freedom to use for the tests. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise.
#' @param nsim when `method` is `"sim"`, the number of iterations used to simulate the p-values. Higher numbers give more accurate p-values subject to less Monte Carlo error but are slower and require more memory. Default is 1,000,000.
#' @param \dots when `method` is `"imhof"`, `"davies"`, or `"liu"`, further arguments passed to \pkgfun{CompQuadForm}{imhof}, \pkgfun{CompQuadForm}{davies}, or \pkgfun{CompQuadForm}{liu}, respectively.
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
#' The test involves computing a test statistic, specifying its distribution under the null hypothesis, and computing the p-value using the compliment of the cumulative density function of the distribution evaluated at the test statistic value. The test statistic depends on `"hypothesis"`. For `hypothesis` equal to a constant, say, \eqn{c}, the test statistic is
#'  \deqn{T^* = \int_\mathcal{A} (\theta(a) - c)^2 \ da}
#' Otherwise, the test statistic is
#' \deqn{T^* = \iint_\mathcal{A} (\theta(a) - \hat{\theta}_0(a))^2 \ da}
#' where \eqn{\hat{\theta}_0} is the projection of \eqn{\theta(a)} onto the null subspace specified by `hypothesis`.
#'
#' Each of these can be approximated as a quadratic form, \eqn{T=\mathbf{\Theta}'\mathbf{W}\mathbf{\Theta}} where \eqn{\mathbf{\Theta}} is a vector of estimates at evaluation points along the curve and \eqn{\mathbf{W}} is a diagonal matrix of weights implementing a trapezoidal approximation to the integral. The null hypothesis is that \eqn{T=0}, which approximates the null hypothesis that \eqn{T^*=0}. Each of the allowable options to `method` corresponds to a method of approximating the distribution of \eqn{T} under the null hypothesis:
#' * `"sim"` simulates the null distribution by simulating from a multivariate normal distribution under the null hypothesis and computing the test statistic in each simulation. The p-value is the proportion of simulated estimates greater than the observed test statistic. When `df` is not `Inf`, the simulation is done from a multivariate t-distribution.
#' * `"imhof"`, `"davies"`, and `"liu"` assume the test statistic follows a generalized \eqn{\chi^2}-distribution and approximate its CDF numerically. `"imhof"` tends to be the most accurate and is recommended. These methods require the \CRANpkg{CompQuadForm} package to be installed.
#' * `"satterthwaite"` also assumes the test statistic follows a generalized \eqn{\chi^2}-distribution, but this distribution is approximated using a scaled F-distribution with the same first two moments.
#' * `"saddlepoint"` also assumes the test statistic follows a generalized \eqn{\chi^2}-distribution, and this distribution is approximated using a saddlepoint method as implemented in \pkgfun{survey}{pFsum}. This method requires the \CRANpkg{survey} package to be installed.
#'
#' In general, we recommend using `method = "imhof"`, though this requires \pkg{CompQuadForm} to be installed (and is the default when it is). `method = "saddlepoint"` has also been shown to be quite accurate and very fast. When using `"sim"`, increasing `nsim` improves the accuracy of the p-value by reducing Monte Carlo error. The default value of `1e6` ensures that the simulated p-value is within .0005 of the true p-value with at least 98\% confidence for p-values less than .05. For greater precision at higher p-value thresholds, `nsim` should be increased. The other methods exist primarily for comparison purposes. `"satterthwaite"` is particularly fast and doesn't require any other packages, but it can be somewhat inaccurate.
#'
#' ## Transform
#'
#' When the effect curve is an ADRF and the outcomes are bounded (e.g., a probability between 0 and 1), the `transform` argument can be specified, which changes the details of the tests.
#'
#' The tests above assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. `transform` transforms the estimates to ones that are unbounded and computes the corresponding distribution of transformed estimates using the delta method. By default, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simple a transformation of the estimated points along the curve already computed. When `hypothesis` is a number, that number is also transformed.
#'
#' Tests on the transformed and untransformed ADRFs correspond to different hypotheses; the difference is not simply a matter of the appropriate distribution of the statistic. For example, for binary outcome model with a logistic transformation, testing that the transformed ADRF is linear corresponds to testing whether the ADRF is a sigmoid function, whereas testing that the untransformed ADRF is linear corresponds to testing whether the ADRF is a straight line. These choices correspond to how the projection of the ADRF is formed; see [curve_projection()] for details. See Examples for an example with a binary outcome. In the example, there is evidence to reject that the transformed ADRF is linear for one of the groups, indicating that a sigmoid function is not sufficient for describing the ADRF, but there is indeterminate evidence (at the .05) to reject that the untransformed ADRF is linear for either group, indicating that a linear function is sufficient for describing the ADRF (at least in the treatment range examined).
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
    .err('{.fun summary} cannot be used when `vcov = "none"` in the original call to {.fun adrf}')
  }

  if (missing(hypothesis) || is_null(hypothesis)) {
    hypothesis <- if (.is_pure_adrf(object)) "flat" else 0
  }
  else if (chk::vld_string(hypothesis)) {
    hypothesis <- match_arg(hypothesis, c("flat", "linear", "quadratic", "cubic"))
  }
  else if (rlang::is_formula(hypothesis)) {
    if (!rlang::is_formula(hypothesis, lhs = FALSE)) {
      .err("if {.arg hypothesis} is a formula, it must be a one-sided formula with the projection model on the right-hand side")
    }

    .treat <- .attr(object, ".treat")

    #Check that only treat is in formula
    vars_in_formula <- get_varnames(hypothesis)

    if (!all(get_varnames(hypothesis) %in% .treat)) {
      .err("only the treatment variable `{(.treat)}` is allowed to appear in {.arg hypothesis} when supplied as a formula")
    }
  }
  else if (!chk::vld_number(hypothesis)) {
    .err("{.arg hypothesis} must be a string, formula, or a number")
  }

  if (missing(method) || is_null(method)) {
    if (rlang::is_installed("CompQuadForm")) {
      method <- "imhof"
    }
    else if (rlang::is_installed("survey")) {
      method <- "saddlepoint"
    }
    else {
      method <- "sim"
    }
  }

  chk::chk_string(method)
  method <- tolower(method)
  method <- match_arg(method, c("sim", "imhof", "davies", "liu",
                                "satterthwaite", "saddlepoint"))

  if (method %in% c("imhof", "davies", "liu")) {
    rlang::check_installed("CompQuadForm")
  }
  else if (method == "saddlepoint") {
    rlang::check_installed("survey")
  }

  ev_tol <- 1e-10
  if (method %in% c("imhof", "davies", "liu", "sim", "satterthwaite",
                    "saddlepoint")) {
    chk::chk_number(ev_tol)
    chk::chk_gt(ev_tol, 0)
  }

  if (method == "sim") {
    chk::chk_count(nsim)
    chk::chk_gt(nsim, 100)
  }
  else {
    nsim <- NULL
  }

  df <- df %or% .attr(object, ".df")

  chk::chk_number(df)
  chk::chk_gt(df, 0)

  .est0 <- .attr(object, ".est")
  .vcov0 <- .attr(object, ".vcov")
  .by_grid <- .attr(object, ".by_grid")
  .contrast <- .attr(object, ".contrast")
  .draws <- .attr(object, ".draws")
  .boot <- .attr(object, ".boot")
  .values <- .attr(object, ".values")

  # transform
  transform_list <- process_transform(transform, object, .est0)

  transformed <- !all(check_if_zero(transform_list$transform(.est0) - .est0))

  # subset
  if (is_null(.by_grid) || is_not_null(.contrast)) {
    .s <- seq_along(.est0)
  }
  else {
    .subset <- process_subset_by_grid(substitute(subset),
                                      .by_grid = .by_grid,
                                      .contrast = .contrast)

    .s <- which(rep(.subset, each = length(.est0) / nrow(.by_grid)))

    .by_grid <- ss(.by_grid, .subset)

    .est0 <- .est0[.s]
    .vcov0 <- ss(.vcov0, .s, .s)
  }

  .values <- .attr(object, ".values")

  n <- length(.values)

  by_list <- {
    if (is_null(.by_grid) && is_null(.contrast)) list(seq_len(n))
    else if (is_not_null(.contrast)) gsplit(seq_along(.est0), rep(seq_along(.contrast), each = n))
    else gsplit(seq_along(.est0), rep(seq_row(.by_grid), each = n))
  }

  if (chk::vld_string(hypothesis)) {
    mm <- switch(hypothesis,
                 flat = matrix(1, nrow = n, ncol = 1L),
                 linear = cbind(1, .values),
                 quadratic = cbind(1, poly(.values, 2)),
                 cubic = cbind(1, poly(.values, 3)))
  }
  else if (rlang::is_formula(hypothesis)) {
    proj_data <- data.frame(.values) |>
      setNames(.treat) |>
      model.frame(formula = hypothesis)

    mt <- .attr(proj_data, "terms")

    mm <- model.matrix(mt, data = proj_data)
  }

  #Weights for trapezoidal Riemann sum
  w <- get_trapezoidal_w(.values)

  wsq <- sqrt(w)

  p <- vapply(by_list, function(.b) {

    if (is_not_null(.boot)) {
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
        resid_boot[["t"]] <- resid_boot[["t"]] |>
          dapply(MARGIN = 1L, return = "matrix",
                 FUN = function(.esti) {
                   fiti <- .nls_fit(x = mm, y = .esti[.b], w = w,
                                    transform = transform_list,
                                    start = fit$coefficients)

                   .esti[.b] - fiti$fitted.values
                 }) |>
          unname()

        vcov.b <- stats::vcov(resid_boot)
      }
    }
    else if (is_not_null(.draws)) {
      if (is.numeric(hypothesis)) {
        est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(hypothesis)

        resid_draws <- transform_list$transform(ss(.draws, j = .b)) - transform_list$transform(hypothesis)

        vcov.b <- stats::cov(resid_draws)
      }
      else {
        fit <- .nls_fit(x = mm, y = .est0[.b], w = w,
                        transform = transform_list)

        resid_draws <- .draws |>
          dapply(MARGIN = 1L, return = "matrix",
                 FUN = function(.esti) {
                   fiti <- .nls_fit(x = mm, y = .esti[.b], w = w,
                                    transform = transform_list,
                                    start = fit$coefficients)

                   .esti[.b] - fiti$fitted.values
                 }) |>
          unname()

        est.b <- .est0[.b] - fit$fitted.values

        vcov.b <- stats::cov(resid_draws)
      }
    }
    else {
      if (is.numeric(hypothesis)) {
        est.b <- transform_list$transform(.est0[.b]) - transform_list$transform(hypothesis)
        vcov.b <- quad_mult(diag(transform_list$d_transform(.est0[.b])),
                            ss(.vcov0, .b, .b))
      }
      else {
        fit <- .nls_fit(x = mm, y = .est0[.b], w = w,
                        transform = transform_list,
                        jac = TRUE, tol = 0, maxit = 20)

        J  <- fit$jacobian
        JW <- w * J

        H <- J %*% solve(crossprod(J, JW), t(JW))

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

    stat <- sum((est.b * wsq)^2)

    if (stat / n < 1e-10) {
      return(1)
    }

    vcov.b <- quad_mult(diag(wsq), vcov.b)

    eig_vcov <- eigen(vcov.b, symmetric = TRUE)

    pos <- eig_vcov$values > ev_tol  # Tolerance for nonnegative EVs

    if (!any(pos)) {
      return(NaN)
    }

    U_r <- ss(eig_vcov$vectors, j = pos)

    Lambda_r_sqrt <- diag(sqrt(eig_vcov$values[pos]), nrow = sum(pos))
    L <- U_r %*% Lambda_r_sqrt

    V <- crossprod(L)

    if (method == "sim") {
      return(.sim_p_value(stat, V, n = nsim, max_size = 1e8, df = df))
    }

    if (method == "satterthwaite") {
      # Moment matching
      mu <- sum(diag(V)) #sum(ev)
      s2 <- sum(V^2)     #sum(ev^2)

      ndf <- mu^2 / s2

      return(pf(stat / mu, df1 = ndf, df2 = df, lower.tail = FALSE))
    }

    #Rescale variance to make stat = 1 (improves numerical performance)
    V <- V / stat

    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values

    ev <- ev[abs(ev) > ev_tol * max(abs(ev))]

    if (method == "saddlepoint") {
      return(survey::pFsum(1, alloc(1.0, length(ev)), ev, df, lower.tail = FALSE,
                           method = "saddlepoint"))
    }

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

    suppressWarnings({
      switch(method,
             "imhof" = CompQuadForm::imhof(q = .q, lambda = .lambda,
                                           h = .h,
                                           epsabs = ...get("epsabs", 1e-6),
                                           epsrel = ...get("epsrel", 1e-6),
                                           limit = ...get("limit", 1e4))$Qq,
             "davies" = CompQuadForm::davies(q = .q, lambda = .lambda,
                                             h = .h,
                                             lim = ...get("lim", 1e4),
                                             acc = ...get("acc", 0.0001))$Qq,
             "liu" = CompQuadForm::liu(q = .q, lambda = .lambda,
                                       h = .h))
    })

  }, numeric(1L))

  res_names <- "p.value"

  out <- make_df(res_names, length(p)) |>
    ftransform(p.value = pmax(p, 1e-12)) |>
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
                 if (is_not_null(.contrast)) "contrast",
                 if (is_not_null(.reference)) sprintf("difference from reference (%s = %s)",
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

  cat(c(out, ""), sep = "\n")

  .print_estimate_table(x, digits = digits, topn = Inf, bar = FALSE)

  method <- switch(method,
                   "satterthwaite" = "the Satterthwaite approximation",
                   "imhof" = "the Imhof approximation",
                   "davies" = "the Davies approximation",
                   "liu" = "the Liu approximation",
                   "sim" = sprintf("a simulation approximation with %s replicates",
                                   prettyNum(nsim, big.mark = ",", scientific = FALSE)),
                   "saddlepoint" = "Kuonen's saddlepoint approximation")

  cat(c(txtbar(width),
        .it(sprintf("Computed using %s", method)) |>
          cli::ansi_strwrap(width, exdent = 2L)), sep = "\n")

  invisible(x)
}
