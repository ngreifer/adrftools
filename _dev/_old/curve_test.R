#' Test Omnibus Hypotheses About Effect Curves
#'
#' @description
#' `curve_test()` tests an omnibus hypothesis about an effect curve. For example, it can be used to test that the ADRF is flat, that the contrast between two ADRFs is 0 everywhere, or that the AMEF is 0 everywhere.
#'
#' @param x an object of class `adrf_est` or `amef_est`; the result of a call to [adrf()], [amef()], or [curve_contrast()].
#' @param hypothesis the hypothesis to be tested. Allowable options include `"flat"` (the default), `"linear"`, `"quadratic"`, `"cubic"`, or a single number (e.g., 0). See Details.
#' @param method string; the method used to compute the p-value of the test. Allowable options include `"sim"` for simulation-based, `"imhof"` for the generalized chi-square test using Imhof's approximation, `"davies"` for the generalized chi-square test using Davie's approximation, `"liu"` for the generalized chi-square test using Liu's approximation, and `"patnaik"` for the generalized chi-square test using Patnaik's approximation. Default is `"sim"`. See Details.
#' @param transform whether to perform the test on the transformed estimates. Allowable options include `TRUE`, `FALSE`, or a function specifying a transformation. Ignored unless `x` is an `adrf_est` object. See Details.
#' @param df the "denominator" degrees of freedom to use for the tests. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise. This only affects p-values when `method` is `"sim"` or `"patnaik"`.
#' @param nsim when `method` is `"sim"`, the number of iterations used to simulate the p-values. Higher numbers give more accurate p-values subject to less Monte Carlo error but are slower and require more memory. Default is 1,000,000.
#' @param \dots when `method` is `"imhof"`, `"davies"`, or `"liu"`, further arguments passed to \pkgfun{CompQuadForm}{imhof}, \pkgfun{CompQuadForm}{davies}, or \pkgfun{CompQuadForm}{liu}, respectively.
#'
#' @returns
#' An object of class `"curve_test"` with elements:
#' \describe{
#'   \item{p.value}{The p-value of the test.}
#'   \item{method}{The method used.}
#'   \item{hypothesis}{The hypothesis tested.}
#' }
#'
#' When `by` is supplied in the original call to [adrf()] or [amef()], a p-value is provided for each stratum of the `by` variables.
#'
#' The output has a `print()` method that displays information about the test.
#'
#' @details
#' `curve_test()` performs an omnibus test for an effect curve. The hypothesis tested is determined by the argument to `"hypothesis"`. The default, `"flat"`, tests whether all values on the curve are equal to each other (i.e., whether the curve is flat), without specifying what value they are equal to. This is equivalent to testing whether the variance around the mean estimate is different from 0. `"linear"` tests whether the curve is linear, i.e., whether the residuals around linear projection of the curve are flat. `"quadratic"` and `"cubic"` test whether the curve is quadratic or cubic, respectively, using the same method. Rejecting the null hypothesis means that the curve is more complicated than the specified function. For example, rejecting the null hypothesis that the curve is linear implies that the curve is nonlinear (and, therefore, not flat either). A number tests whether all values on the curve are equal to that number.
#'
#' The test involves computing a test statistic, specifying its distribution under the null hypothesis, and computing the p-value using the compliment of the cumulative density function of the distribution evaluated at the test statistic value. The test statistic depends on `"hypothesis"`. For `hypothesis` equal to a constant, say, \eqn{c}, the test statistic is
#'  \deqn{T^* = \int_\mathcal{A} (\theta(a) - c)^2 \ da}
#' Otherwise, the test statistic is
#' \deqn{T^* = \iint_\mathcal{A} (\theta(a) - \hat{\theta}_0(a))^2 \ da}
#' where \eqn{\hat{\theta}_0} is the projection of \eqn{\theta(a)} onto the null subspace specified by `hypothesis`.
#'
#' Each of these can be approximated as a quadratic form, \eqn{T=\mathbf{\Theta}'\mathbf{W}\mathbf{\Theta}} where \eqn{\mathbf{\Theta}} is a vector of estimates at evaluation points along the curve (corresponding to `values` or `n` supplied to [adrf()] or [amef()]) and \eqn{\mathbf{W}} is a symmetric matrix that approximates the integral using a trapezoidal approximation. The null hypothesis is that \eqn{T=0}, which approximates the null hypothesis that \eqn{T^*=0}. Each of the allowable options to `method` corresponds to a method of approximating the distribution of \eqn{T} under the null hypothesis:
#' * `"sim"` simulates the null distribution by simulating from a multivariate normal distribution under the null hypothesis and computing the test statistic in each simulation. The p-value is the proportion of simulated estimates greater than the observed test statistic. When `df` is not `Inf`, the simulation is done from a multivariate t-distribution.
#' * `"imhof"`, `"davies"`, and `"liu"` assume the test statistic follows a generalized \eqn{\chi^2}-distribution and approximate its CDF numerically. `"imhof"` tends to be the most accurate and is recommended. These methods require the \CRANpkg{CompQuadForm} package to be installed. This approximation can be off in small samples due to treating the parameter covariance matrix as if it were known rather than estimated.
#' * `"patnaik"` also assumes the test statistic follows a generalized \eqn{\chi^2}-distribution, but this distribution is approximated using a scaled \eqn{\chi^2}-distribution with the same first two moments. This approximation is fast to compute, does not require any other package, and is not subject to Monte Carlo error, but it can be inaccurate.
#'
#' In general, we recommend `"sim"` in small to moderate samples and `"imhof"` in large samples. When using `"sim"`, increasing `nsim` improves the accuracy of the p-value by reducing Monte Carlo error. The default value of `1e6` ensures that the simulated p-value is within .0005 of the true p-value with at least 98\% confidence for p-values less than .05. For greater precision at higher p-value thresholds, `nsim` should be increased.
#'
#' ## Transform
#'
#' The tests above assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. `transform` transforms the estimates to ones that are unbounded and computes the corresponding distribution of transformed estimates using the delta method. By default with `adrf_est` objects, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simple a transformation of the estimated points along the curve already computed. When `hypothesis` is a number, that number is also transformed.
#'
#' @seealso
#' * [summary.adrf_est()] to summarize and plot the effect curve and its confidence bands.
#'
#' @examples
#'

#' @export
curve_test <- function(x, hypothesis = "flat", method = "sim", transform = TRUE, df = NULL,
                       nsim = 1e6, ...) {

  chk::chk_not_missing(x, "`x`")

  .chk_is(x, c("adrf_est", "amef_est"))

  if (chk::vld_string(hypothesis)) {
    hypothesis <- match_arg(hypothesis, c("flat", "linear", "quadratic", "cubic"))
  }
  else if (!chk::vld_number(hypothesis)) {
    .err("`hypothesis` must be a string or a number")
  }

  vcov_type <- x[["vcov_type"]]

  if (is_null(vcov_type)) {
    .err(sprintf('`curve_test()` cannot be used when `vcov = "none"` in the original call to `%s()`',
                 if (inherits(x, "adrf_est")) "adrf" else "amef"))
  }

  chk::chk_string(method)
  method <- tolower(method)
  method <- match_arg(method, c("sim", "patnaik", "imhof", "davies", "liu"))

  if (is_null(df)) {
    df <- get_df(.attr(x, "model"))
  }
  else {
    chk::chk_number(df)
    chk::chk_gt(2)
  }

  if (method %in% c("imhof", "davies", "liu")) {
    rlang::check_installed("CompQuadForm")
  }

  ev_tol <- 1e-10
  if (method %in% c("imhof", "davies", "liu", "sim", "patnaik")) {
    chk::chk_number(ev_tol)
    chk::chk_gt(ev_tol, 0)
  }
  else if (method == "sim") {
    chk::chk_count(nsim)
    chk::chk_gt(nsim, 100)
  }

  if (inherits(x, "adrf_est") && is_null(x$contrast)) {
    # transform
    transform_list <- process_transform(transform, .attr(x, "model")$family)
    transform <- transform_list$transform

    est <- transform(x$est)

    transformed <- !all(check_if_zero(est - x$est))
  }
  else {
    transformed <- FALSE
  }

  if (transformed) {
    if (vcov_type == "bootstrap") {
      vcov <- cov(transform(x$boot[["t"]][, seq_along(est), drop = FALSE]))
    }
    else {
      d_transformed_est <- transform_list$d_transform(x$est)
      vcov <- quad_mult(diag(d_transformed_est), x$vcov)
    }

    if (chk::vld_number(hypothesis)) {
      hypothesis_ <- transform(hypothesis)
    }
  }
  else {
    est <- x$est
    vcov <- x$vcov

    if (chk::vld_number(hypothesis)) {
      hypothesis_ <- hypothesis
    }
  }

  n <- length(x$values)

  by_list <- {
    if (is_null(x$by_grid) && is_null(x$contrast)) list(seq_len(n))
    else if (is_not_null(x$contrast)) split(seq_along(x$est), rep(seq_along(x$contrast), each = n))
    else split(seq_along(x$est), rep(seq_row(x$by_grid), each = n))
  }

  #Weights for trapezoidal Riemann sum
  w <- get_trapezoidal_w(x$values)

  wsq <- sqrt(w)

  if (chk::vld_string(hypothesis)) {
    Xwsq <- wsq * switch(hypothesis,
                         flat = matrix(1, nrow = n, ncol = 1L),
                         linear = cbind(1, x$values),
                         quadratic = cbind(1, poly(x$values, 2)),
                         cubic = cbind(1, poly(x$values, 3)))
  }

  p <- vapply(by_list, function(.b) {
    if (chk::vld_string(hypothesis)) {
      fit <- .lm.fit(x = Xwsq, y = est[.b] * wsq)

      qr <- fit[c("qr", "qraux", "pivot", "tol", "rank")] |>
        structure(class = "qr")
      qq <- qr.qy(qr, diag(1, nrow = nrow(qr$qr), ncol = qr$rank))
      M <- diag(n) - (tcrossprod(qq) / wsq) %*% diag(wsq) #residual maker

      est.b <- fit$residuals / wsq
      vcov.b <- quad_mult(M, ss(vcov, .b, .b))
    }
    else {
      est.b <- est[.b] - hypothesis_
      vcov.b <- ss(vcov, .b, .b)
    }

    stat <- sum(est.b^2 * w)

    eig_vcov <- eigen(vcov.b, symmetric = TRUE)

    pos <- eig_vcov$values > ev_tol  # Tolerance for nonnegative EVs
    U_r <- eig_vcov$vectors[, pos, drop = FALSE]

    Lambda_r_sqrt <- diag(sqrt(eig_vcov$values[pos]))
    L <- U_r %*% Lambda_r_sqrt

    #Rescale variance to make stat = 1 (improves numerical performance)
    V <- crossprod(L * wsq) / stat

    if (startsWith(method, "sim")) {
      return(.sim_p_value(1, V, nsim, max_size = 1e8, df = df))
    }

    if (method == "patnaik") {
      # Moment matching
      mu <- sum(diag(V)) #sum(ev)
      s2 <- 2 * sum(V^2) #2*sum(ev^2)

      if (is.finite(df)) {
        s2 <- s2 * (1 + 1 / df)
      }

      uu <- 2 * mu / s2

      return(pchisq(uu, df = uu * mu, lower.tail = FALSE))
    }

    ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values

    ev <- ev[abs(ev) > 1e-10]

    suppressWarnings({
      switch(method,
             "imhof" = CompQuadForm::imhof(q = 1, lambda = ev,
                                           h = alloc(1.0, length(ev)),
                                           delta = alloc(0.0, length(ev)),
                                           epsabs = ...get("epsabs", 1e-6),
                                           epsrel = ...get("epsrel", 1e-6),
                                           limit = ...get("limit", 1e4))$Qq,
             "davies" = CompQuadForm::davies(q = 1, lambda = ev,
                                             h = alloc(1.0, length(ev)),
                                             delta = alloc(0.0, length(ev)),
                                             lim = ...get("lim", 1e4),
                                             acc = ...get("acc", 0.0001))$Qq,
             "liu" = CompQuadForm::liu(q = 1, lambda = ev,
                                       h = alloc(1.0, length(ev)),
                                       delta = alloc(0.0, length(ev))))
    })

  }, numeric(1L))

  out <- list(p.value = pmax(p, 1e-12),
              method = method,
              hypothesis = hypothesis,
              df = df)

  attr(out, "obj") <- x
  attr(out, "transformed") <- transformed

  class(out) <- "curve_test"
  out
}

#' @exportS3Method print curve_test
print.curve_test <- function(x, digits = max(4L, getOption("digits") - 3L), ...) {
  out <- c(center_just("Omnibus Curve Test", wrt = space(50L)),
           txtbar(50L))

  is_contrast <- inherits(.attr(x, "obj"), "curve_contrast")

  qoi <- {
    if (is_contrast) {
      if (inherits(.attr(x, "obj"), "adrf_est")) "ADRF difference"
      else "AMEF difference"
    }
    else if (inherits(.attr(x, "obj"), "adrf_est")) "ADRF"
    else if (inherits(.attr(x, "obj"), "amef_est")) "AMEF"
  }

  out <- c(out,
           sprintf("H\u2080: %s is %s for values of %s between %s and %s",
                   qoi,
                   x$hypothesis,
                   .attr(x, "obj")$treat,
                   round(min(.attr(x, "obj")$values), digits),
                   round(max(.attr(x, "obj")$values), digits)) |>
             strwrap(50, exdent = 4))

  if (is_contrast) {
    res <- data.frame(.attr(x, "obj")$contrast,
                      x$p.value) |>
      setNames(c("Contrast", "P-value"))
  }
  else if (is_not_null(.attr(x, "obj")$by_grid)) {
    res <- cbind(.attr(x, "obj")$by_grid,
                 x$p.value)
    names(res)[ncol(res)] <- "P-value"
  }
  else {
    res <- data.frame(x$p.value) |>
      setNames("P-value")
  }

  res[["P-value"]] <- format(res[["P-value"]], digits = digits)

  if (is_not_null(.attr(x, "nsim"))) {
    min_p <- 1 / .attr(x, "nsim")
    res[["P-value"]][x$p.value < min_p] <- sprintf("< %s", min_p)
  }

  tmp <- utils::capture.output({
    print.data.frame(res, row.names = FALSE, digits = digits)
  })

  out <- c(out, "", tmp,
           txtbar(50L),
           .it(sprintf("Computed using %s",
                       .method2phrase(x$method))))

  cat(out, sep = "\n")

  invisible(x)
}

.method2phrase <- function(x) {
  switch(x,
         "patnaik" = "the Patnaik approximation",
         "imhof" = "the Imhof approximation",
         "davies" = "the Davies approximation",
         "liu" = "the Liu approximation",
         "sim" = "a simulation approximation",
         "supt" = "an inversion of the sup-t confidence band")
}
