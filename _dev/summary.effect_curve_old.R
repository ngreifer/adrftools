.summary.effect_curve_old <- function(object, hypothesis, method, subset = NULL, transform = TRUE, df = NULL,
                                      nsim = 1e6, ...) {

  .vcov_type <- .attr(object, ".vcov_type")

  if (is_null(.vcov_type)) {
    .err("`summary()` cannot be used when `vcov = \"none\"` in the original call to `adrf()`")
  }

  if (missing(hypothesis) || is_null(hypothesis)) {
    hypothesis <- if (.is_pure_adrf(object)) "flat" else 0
  }
  else if (chk::vld_string(hypothesis)) {
    hypothesis <- match_arg(hypothesis, c("flat", "linear", "quadratic", "cubic"))
  }
  else if (rlang::is_formula(hypothesis)) {
    if (!rlang::is_formula(hypothesis, lhs = FALSE)) {
      .err("if `hypothesis` is a formula, it must be a one-sided formula with the projection model on the right-hand side")
    }

    .treat <- .attr(object, ".treat")

    #Check that only treat is in formula
    vars_in_formula <- get_varnames(hypothesis)

    if (!all(get_varnames(hypothesis) %in% .treat)) {
      .err(sprintf("only the treatment variable `%s` is allowed to appear in `hypothesis` when supplied as a formula",
                   .treat))
    }
  }
  else if (!chk::vld_number(hypothesis)) {
    .err("`hypothesis` must be a string, formula, or a number")
  }

  if (missing(method) || is_null(method)) {
    if (rlang::is_installed("CompQuadForm")) {
      method <- "imhof"
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

  # transform
  transform_list <- process_transform(transform, object, .est0)

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

  transform <- transform_list$transform
  .est <- transform(.est0)

  transformed <- !all(check_if_zero(.est - .est0))

  if (.vcov_type == "bootstrap") {
    .posterior <- .attr(object, ".boot")[["t"]] |>
      ss(j = .s) |>
      transform()

    if (method != "posterior") {
      .vcov <- stats::cov(.posterior)
    }
  }
  else if (.vcov_type == "posterior") {
    .posterior <- .attr(object, ".draws") |>
      ss(j = .s) |>
      transform()

    if (method != "posterior") {
      .vcov <- stats::cov(.posterior)
    }
  }
  else if (transformed) {
    d_transformed_est <- transform_list$d_transform(.est0)
    .vcov <- quad_mult(diag(d_transformed_est), .vcov0)
  }
  else {
    .vcov <- .vcov0
  }

  .values <- .attr(object, ".values")

  n <- length(.values)

  by_list <- {
    if (is_null(.by_grid) && is_null(.contrast)) list(seq_len(n))
    else if (is_not_null(.contrast)) gsplit(seq_along(.est), rep(seq_along(.contrast), each = n))
    else gsplit(seq_along(.est), rep(seq_row(.by_grid), each = n))
  }

  #Weights for trapezoidal Riemann sum
  w <- get_trapezoidal_w(.values)

  wsq <- sqrt(w)

  if (rlang::is_formula(hypothesis)) {
    proj_data <- data.frame(.values) |>
      setNames(.treat) |>
      model.frame(formula = hypothesis)

    if (is_not_null(model.offset(proj_data))) {
      .err("offsets are not allowed in `hypothesis` when supplied as a formula")
    }

    mt <- .attr(proj_data, "terms")

    if (is_null(.attr(mt, "term.labels")) && .attr(mt, "intercept") == 1) {
      hypothesis <- "flat"
    }

    Xwsq <- wsq * model.matrix(mt, data = proj_data)
  }
  else if (chk::vld_string(hypothesis)) {
    Xwsq <- wsq * switch(hypothesis,
                         flat = matrix(1, nrow = n, ncol = 1L),
                         linear = cbind(1, .values),
                         quadratic = cbind(1, poly(.values, 2)),
                         cubic = cbind(1, poly(.values, 3)))
  }

  p <- vapply(by_list, function(.b) {
    if (is.numeric(hypothesis)) {
      est.b <- .est[.b] - transform(hypothesis)
      vcov.b <- ss(.vcov, .b, .b)
    }
    else {
      fit <- .lm.fit(x = Xwsq, y = .est[.b] * wsq)

      qr <- fit[c("qr", "qraux", "pivot", "tol", "rank")] |>
        structure(class = "qr")
      qq <- qr.qy(qr, diag(1, nrow = nrow(qr$qr), ncol = qr$rank))
      M <- diag(n) - (tcrossprod(qq) / wsq) %*% diag(wsq) #residual maker

      est.b <- fit$residuals / wsq
      vcov.b <- quad_mult(M, ss(.vcov, .b, .b))
    }

    stat <- sum(est.b^2 * w)

    eig_vcov <- eigen(vcov.b, symmetric = TRUE)

    pos <- eig_vcov$values > ev_tol  # Tolerance for nonnegative EVs
    U_r <- ss(eig_vcov$vectors, j = pos)

    Lambda_r_sqrt <- diag(sqrt(eig_vcov$values[pos]))
    L <- U_r %*% Lambda_r_sqrt

    #Rescale variance to make stat = 1 (improves numerical performance)
    V <- crossprod(L * wsq) / stat

    if (method == "sim") {
      return(.sim_p_value(1, V, nsim, max_size = 1e8, df = df))
    }

    if (method == "satterthwaite") {
      # Moment matching
      mu <- sum(diag(V)) #sum(ev)
      s2 <- sum(V^2)     #sum(ev^2)

      ndf <- mu^2 / s2

      return(pf(1 / mu, df1 = ndf, df2 = df, lower.tail = FALSE))
    }

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

  res_names <- c("p.value")

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
