# Multivariate (starting with 2)

.effect_curve_internal_mv <- function(x, treat, vcov, ranges = .95, n,
                                      data, subset, by, wts, cluster = NULL, fwb.args = list(),
                                      estimator = NULL, grad_fun = NULL,
                                      pred_fun = NULL, mcall,
                                      .adrf_env = parent.frame(2L), ...) {

  # get data
  model_data <- process_model_data(x, data)

  # treat
  treat_vars <- process_treat_mv(treat, model_data)

  # wts
  wts <- process_wts(wts, x, data)

  # values
  values_list <- process_range_mv(ranges, n, treat_vars, w = wts)
  ns <- lengths(values_list)

  values <- do.call("expand.grid", c(values_list, list(KEEP.OUT.ATTRS = FALSE)))
  n < nrow(values)

  # bayes
  is_bayes <- .is_bayesian(x, model_data)

  # vcov + cluster
  cl_list <- process_vcov_and_cluster(vcov, x, cluster,
                                      is_bayes, model_data)

  vcov <- cl_list[["vcov"]]
  # vcov <- process_vcov(vcov, x, is_bayes)

  # cluster
  cl_list[["vcov"]] <- NULL
  # cl_list <- process_cluster(cluster, x, vcov, model_data)

  # subset
  s <- process_subset(subset, model_data, env = .adrf_env)
  model_data <- ss(model_data, s)

  # by
  processed_by_list <- process_by(by, model_data)
  by_grid <- processed_by_list$by_grid
  by_id <- processed_by_list$by_id

  # family
  .family <- process_family(x)

  # df
  .df <- process_df(x)

  # response
  .response <- process_response(x)

  .fmean <- collapse::fmean

  ni <- nrow(model_data)
  ii <- seq_len(ni)

  vi <- seq_len(n)
  val_id <- rep(vi, each = ni)

  cross_id <- GRP(list(by_id = rep.int(by_id, n), val_id = val_id),
                  sort = TRUE)

  wi <- rep.int(ss(wts, s), n)

  data_grid <- ss(model_data, rep.int(ii, n))

  data_grid[treat] <- values[val_id, ]

  if (!is_bayes) {
    beta0 <- marginaleffects::get_coef(x)
    aliased <- is.na(beta0)
  }

  if (!is_bayes && is_null(pred_fun)) {
    mm <- get_tt(x, model_data) |>
      model.matrix(data = data_grid) |>
      ss(j = !aliased)

    est_fun <- function(model = NULL, b = NULL, data0 = NULL) {
      estimator(mm = data0 %or% mm, beta0 = b, model = model)
    }
  }
  else {
    est_fun <- function(model = NULL, b = NULL, data0 = NULL) {
      pred_fun(data_grid = data0 %or% data_grid, beta0 = b, model = model, ...)
    }
  }

  vcov_type <- {
    if (identical(vcov, "none")) "none"
    else if (is_bayes) "posterior"
    else if (identical(vcov, "boot") || identical(vcov, "fwb")) "bootstrap"
    else if (identical(vcov, "unconditional")) "unconditional"
    else "conditional"
  }

  .boot <- .draws <- NULL

  if (is_bayes) {
    .draws <- rsplit(seq_row(data_grid), rep.int(val_id, nrow(data_grid) / (ni * n))) |>
      lapply(function(.j) {
        est_fun(x, data0 = ss(data_grid, .j)) |>
          .attr("posterior_draws") |>
          .fmean(g = by_id, w = ss(wi, .j))
      }) |>
      do.call(what = "rbind") |>
      t()

    .est <- colMeans(.draws)

    if (identical(vcov_type, "none")) {
      .vcov <- .draws <- NULL
    }
    else {
      .vcov <- stats::cov(.draws)
    }
  }
  else if (identical(vcov_type, "none")) {
    .est <- est_fun(x) |>
      .fmean(g = cross_id, w = wi)

    .vcov <- NULL
  }
  else if (identical(vcov_type, "bootstrap")) {

    bootfun <- function(data, w = alloc(1.0, ni), ...) {

      wts_boot <- wts * w

      x_boot <- do.call("update", list(x, weights = wts_boot))

      wi_boot <- {
        if (ni == 1L) alloc(1.0, n)
        else rep.int(ss(wts_boot, s), n)
      }

      est_fun(x_boot) |>
        .fmean(g = cross_id, w = wi_boot)
    }

    if (is_null(fwb.args)) {
      fwb.args <- list()
    }

    chk::chk_list(fwb.args)

    if (identical(vcov, "boot")) {
      fwb.args[["wtype"]] <- "multinom"
    }

    if (is_not_null(cluster)) {
      fwb.args[["cluster"]] <- cl_list[["cluster"]]
    }

    .boot <- do.call(fwb::fwb, c(list(data = process_model_data(x),
                                      statistic = bootfun,
                                      drop0 = FALSE),
                                 fwb.args))

    .est <- stats::coef(.boot)
    .vcov <- stats::vcov(.boot)
  }
  else if (identical(vcov_type, "unconditional")) {
    beta0 <- beta0[!aliased]

    phat <- est_fun(x)

    .est <- .fmean(phat, g = cross_id, w = wi)

    gr <- {
      if (is_null(grad_fun) || is_not_null(pred_fun)) {
        .gradient(function(b) {
          est_fun(x, b) |>
            .fmean(g = cross_id, w = wi) |>
            drop()
        }, .x = beta0)
      }
      else {
        grad_fun(mm, eta = drop(mm %*% beta0)) |>
          fmean(g = cross_id, w = wi)
      }
    }

    if (inherits(x, "svyglm")) {
      infl <- ni * (.attr(x, "influence") %or% tcrossprod(model.matrix(x) * resid(x, "working") * x$weights,
                                                          x$naive.cov))
    }
    else {
      infl <- sandwich::estfun(x) |>
        tcrossprod(sandwich::bread(x))
    }

    S <- tcrossprod(infl, gr)

    if (ni > 1L) {
      wts_sum <- sum(wts)

      phat_split <- gsplit(phat, cross_id)

      b_split <- gsplit(ii, by_id)

      for (k in seq_along(phat_split)) {
        b <- s[b_split[[cross_id$groups$by_id[k]]]]
        S[b, k] <- S[b, k] + (phat_split[[k]] - .est[k]) * wts_sum / sum(wts[b])
      }
    }

    if (is_null(cl_list)) {
      .vcov <- crossprod(S / nrow(S))
    }
    else {
      .vcov <- Reduce("+", lapply(seq_along(cl_list[["cluster"]]), function(j) {
        crossprod(fsum(S, cl_list[["cluster"]][[j]]) / nrow(S)) * cl_list[["adj"]][j]
      }))
    }
  }
  else {
    .est <- est_fun(x) |>
      .fmean(g = cross_id, w = wi)

    V <- marginaleffects::get_vcov(x, vcov = vcov)

    gr <- {
      if (is_null(grad_fun) || is_not_null(pred_fun)) {
        .gradient(function(b) {
          est_fun(x, b) |>
            .fmean(g = cross_id, w = wi) |>
            drop()
        }, .x = beta0)
      }
      else {
        grad_fun(mm, eta = drop(mm %*% beta0)) |>
          .fmean(g = cross_id, w = wi)
      }
    }

    .vcov <- quad_mult(gr, V)
  }

  attr(vcov_type, "clustered") <- is_not_null(cl_list[["cluster"]])

  .make_fun_mv(.est = drop(unname(.est)),
               .vcov = .vcov,
               .values = values,
               .treat = treat,
               .vcov_type = vcov_type,
               .boot = .boot,
               .draws = .draws,
               .curve_type = "ADRF",
               .family = .family,
               .df = .df,
               .response = .response,
               .by_grid = by_grid)
}

process_treat_mv <- function(treat, data) {
  chk::chk_character(treat)
  chk::chk_subset(treat, names(data))

  treat_vars <- data[treat]

  if (!all_apply(treat_vars, chk::vld_numeric)) {
    .err("{.arg treat} must be the names of one or more numeric variables in the original dataset")
  }

  treat_vars
}

process_range_mv <- function(ranges = .95, n, treat_vars, w = NULL, strict = FALSE) {
  z <- ncol(treat_vars)

  if (!is.list(ranges)) {
    chk::chk_numeric(ranges)

    if (z == 1L) {
      if (length(ranges) > 2L) {
        .err("when supplied as a vector, {.arg range} must have length 1 or 2")
      }

      ranges <- list(ranges) |>
        setNames(names(treat_vars))
    }
    else if (length(ranges) == 1L) {
      ranges <- rep.int(list(ranges), z) |>
        setNames(names(treat_vars))
    }
    else {
      .err("{.arg ranges} must be a list with an entry for each treatment")
    }
  }

  if (length(ranges) == 1L && z > 1L) {
    ranges <- rep.int(ranges, z) |>
      setNames(names(treat_vars))
  }

  if (length(ranges) != z) {
    .err("{.arg ranges} must be a list with an entry for each treatment")
  }

  if (is_null(names(ranges))) {
    names(ranges) <- names(treat_vars)
  }
  else if (!setequal(unique(names(ranges)), unique(names(treat_vars)))) {
    .err("the names of each entry in {.arg ranges} must correspond to one of the supplied treatments")
  }

  ranges <- ranges[names(treat_vars)]

  values <- make_list(names(ranges))

  for (i in names(treat_vars)) {
    if (!is.numeric(ranges[[i]])) {
      .err("each entry in {.arg ranges} must be numeric")
    }

    if (length(ranges[[i]]) > 1L) {
      ranges[[i]] <- sort(ranges[[i]])

      range_t <- .range(treat_vars[[i]])

      if (strict) {
        if (any(ranges[[i]] < range_t[1L]) || any(ranges[[i]] > range_t[2L])) {
          r <- format(range_t, digits = 4L, drop0trailing = TRUE)
          .err("no values in the {.arg ranges} entry for {i} can be outside the range of its observed observed treatment values ({r[1L]} to {r[2L]})")
        }
      }
      else {
        range_diff <- range_t[2L] - range_t[1L]

        if (any(ranges[[i]] < range_t[1L] - .1 * range_diff) ||
            any(ranges[[i]] > range_t[2L] + .1 * range_diff)) {
          .err("some values in the {.arg ranges} entry for {i} are outside the range of observed treatment values")
        }
      }
    }
    else if (ranges[[i]] == 1) {
      ranges[[i]] <- .range(treat_vars[[i]])
    }
    else if (chk::vld_range(ranges[[i]], c(0, 1), inclusive = FALSE)) {
      ranges[[i]] <- .quantile(treat_vars[[i]], probs = c(1 - ranges[[i]], 1 + ranges[[i]]) / 2,
                               w = w)
    }
    else {
      .err("if an entry in {.arg ranges} is a single number, it must be a quantile (i.e., between 0 and 1)")
    }

    values[[i]] <- seq(ranges[[i]][1L], ranges[[i]][2L], length.out = n)
  }

  values
}

.make_fun_mv <- function(fn = NULL, .contrast_mat = NULL,
                         .est, .vcov, .values, .treat, .vcov_type, .boot, .draws,
                         .curve_type, .family, .df, .response, .call,
                         .by_grid, .contrast, .reference) {

  if (missing(.curve_type) && is_not_null(fn)) {
    if (inherits(fn, "adrf_curve")) {
      .curve_type <- "ADRF"
    }
    else if (inherits(fn, "amef_curve")) {
      .curve_type <- "AMEF"
    }
  }

  for (i in rlang::fn_fmls_names()[-(1:2)]) {
    if (rlang::is_missing(environment()[[i]])) {
      if (is_null(fn)) {
        assign(i, NULL)
      }
      else {
        assign(i, .attr(fn, i))
      }
    }
  }

  if (is_not_null(.contrast_mat)) {
    # est is the contrast estimates
    .est <- drop(tcrossprod(.est, .contrast_mat))

    if (.vcov_type == "bootstrap") {
      .boot[["t"]] <- .boot[["t"]] |>
        ss(j = seq_along(.est)) |>
        tcrossprod(.contrast_mat)

      .boot[["t0"]] <- .boot[["t0"]][seq_along(.est)] |>
        tcrossprod(.contrast_mat) |>
        drop()

      .vcov <- stats::vcov(.boot)
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = seq_along(.est)) |>
        tcrossprod(.contrast_mat)

      .vcov <- stats::cov(.draws)
    }
    else if (.vcov_type != "none") {
      .vcov <- quad_mult(.contrast_mat, .vcov)
    }
  }

  fn <- function(..., subset = NULL) {
    values <- ...mget(.treat)

    for (.t in .treat) {
      if (is_null(values[[.t]])) {
        .err("values must be supplied for {(.t)}")
      }

      chk::chk_numeric(values[[.t]], add_quotes(.t, "`"))
    }

    values <- do.call("expand.grid", c(values[.treat], list(KEEP.OUT.ATTRS = FALSE)))

    # subset
    if (any(rlang::fn_fmls_names() == "subset") && is_null(.contrast)) {
      .subset <- process_subset_by_grid(substitute(subset),
                                        .by_grid = .by_grid,
                                        .contrast = .contrast)

      .s <- which(rep(.subset, each = length(.est) / nrow(.by_grid)))

      .by_grid <- ss(.by_grid, .subset)

      .est <- .est[.s]
    }
    else {
      .s <- seq_along(.est)
    }

    n <- nrow(values)

    val_mat <- matrix(0, nrow = n, ncol = nrow(.values))

    ## NEED TO DO MULTIVARIATE INTERPOLATION
    matches <- max.col(-abs(outer(values, .values, "-")), "first")

    to_interp <- !check_if_zero(values - .values[matches])

    if (any(to_interp)) {
      val_mat[to_interp, ] <- get_locpoly_w(x = values[to_interp,],
                                            v = .values)
    }

    if (!all(to_interp)) {
      val_mat[cbind(which(!to_interp), matches[!to_interp])] <- 1
    }

    n_by <- get_n_by(.contrast, .by_grid)

    val_mat <- block_diag(val_mat, n_by)

    # est0 is the estimates at values, including interp'd
    .est0 <- drop(val_mat %*% .est)

    res_names <- c(.treat, "estimate")
    res <- make_df(res_names, length(.est0))

    res$estimate <- .est0

    res <- add_est_labels_mv(res, .contrast, .by_grid, values, .treat)

    if (.vcov_type == "bootstrap") {
      .boot[["t"]] <- .boot[["t"]] |>
        ss(j = .s) |>
        tcrossprod(val_mat)

      .boot[["t0"]] <- .boot[["t0"]][.s] |>
        tcrossprod(val_mat) |>
        drop()

      attr(res, ".vcov") <- stats::vcov(.boot)
      attr(res, ".boot") <- .boot
    }
    else if (.vcov_type == "posterior") {
      .draws <- .draws |>
        ss(j = .s) |>
        tcrossprod(val_mat)

      attr(res, ".vcov") <- stats::cov(.draws)
      attr(res, ".draws") <- .draws
    }
    else if (.vcov_type != "none") {
      attr(res, ".vcov") <- quad_mult(val_mat, ss(.vcov, .s, .s))
    }

    attr(res, "values") <- values

    attr(res, ".treat") <- .treat
    attr(res, ".vcov_type") <- .vcov_type
    attr(res, ".curve_type") <- .curve_type
    attr(res, ".by_grid") <- .by_grid
    attr(res, ".contrast") <- .contrast
    attr(res, ".family") <- .family
    attr(res, ".df") <- .df
    attr(res, ".reference") <- .reference

    class(res) <- unique(c("curve_est", class(res)))

    res
  }

  # rlang::fn_fmls_names(fn)[1L] <- .treat

  if (is_null(.by_grid) || is_not_null(.contrast)) {
    rlang::fn_fmls(fn)[length(rlang::fn_fmls(fn))] <- NULL
  }

  attr(fn, ".est") <- .est
  attr(fn, ".vcov") <- .vcov
  attr(fn, ".values") <- .values
  attr(fn, ".treat") <- .treat
  attr(fn, ".vcov_type") <- .vcov_type
  attr(fn, ".boot") <- .boot
  attr(fn, ".draws") <- .draws
  attr(fn, ".reference") <- .reference
  attr(fn, ".by_grid") <- .by_grid
  attr(fn, ".contrast") <- .contrast
  attr(fn, ".family") <- .family
  attr(fn, ".df") <- .df
  attr(fn, ".family") <- .family
  attr(fn, ".response") <- .response
  attr(fn, ".call") <- .call

  classes <- c("effect_curve", class(fn))

  if (.curve_type == "ADRF") {
    classes <- c("adrf_curve", classes)
  }

  if (.curve_type == "AMEF") {
    classes <- c("amef_curve", classes)
  }

  if (is_not_null(.reference)) {
    classes <- c("reference_curve", classes)
  }

  if (is_not_null(.contrast)) {
    classes <- c("contrast_curve", classes)
  }

  class(fn) <- unique(classes)

  fn
}

add_est_labels_mv <- function(res, .contrast, .by_grid, .values = NULL, est_names = NULL) {

  add_values <- is_not_null(.values) && is_not_null(est_names)

  if (!add_values) {
    .values <- data.frame(NA)
  }

  # Apply by_grid and contrast labels if necessary
  if (is_not_null(.contrast)) {
    contrast <- factor(rep(.contrast, each = nrow(.values)),
                       levels = .contrast)

    if (add_values) {
      res[est_names] <- ss(.values, rep.int(seq_row(.values), length(.contrast)))
    }

    res <- cbind(contrast = contrast, res)
  }
  else if (is_not_null(.by_grid)) {
    res_names <- names(res)
    val_id <- seq_row(.values)
    by_id <- seq_row(.by_grid)

    res <- res |>
      ftransform(..vi = rep.int(val_id, length(by_id)),
                 ..by_id = rep(by_id, each = length(val_id))) |>
      merge(cbind(.by_grid, ..by_id = by_id),
            by = "..by_id",
            all.x = TRUE, all.y = FALSE,
            sort = FALSE) |>
      roworderv(c("..by_id", "..vi"))

    if (add_values) {
      res[est_names] <- ss(.values, res[["..vi"]])
    }

    res <- ss(res, j = c(names(.by_grid), res_names))
  }
  else if (add_values) {
    res[est_names] <- .values
  }

  res
}

curve_contrast_mv <- function(x) {
  chk::chk_not_missing(x, "`x`")
  check_effect_curve(x, contrast_ok = FALSE, projection_ok = FALSE)

  .by_grid <- .attr(x, ".by_grid")

  if (is_null(.by_grid)) {
    fn_name <- rlang::current_call() |> rlang::call_name()
    .err("{.fun adrf} must have been called with {.arg by} specified to use {.fun {fn_name}}")
  }

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")

  # All combinations of 'by' levels
  combos <- utils::combn(seq_row(.by_grid), 2L, simplify = FALSE)

  comp_grid <- expand.grid(val_id = seq_row(.values),
                           combo = seq_along(combos),
                           KEEP.OUT.ATTRS = FALSE)

  m <- matrix(0, nrow = nrow(comp_grid), ncol = length(.est))

  by_id <- rep(seq_row(.by_grid), each = nrow(.values))

  for (co in seq_along(combos)) {
    diag(m[comp_grid$combo == co, by_id == combos[[co]][2L]]) <- 1
    diag(m[comp_grid$combo == co, by_id == combos[[co]][1L]]) <- -1
  }

  .contrast <- vapply(combos, function(combo) {
    sprintf("[%s] - [%s]",
            get_by_grid_labels(ss(.by_grid, combo[2L])),
            get_by_grid_labels(ss(.by_grid, combo[1L]))
    )
  }, character(1L))

  .make_fun_mv(x,
               .contrast_mat = m,
               .contrast = .contrast)
}

curve_projection_mv <- function(x, model, transform = TRUE) {
  chk::chk_not_missing(x, "`x`")
  check_effect_curve(x, projection_ok = FALSE)

  chk::chk_not_missing(model, "`model`")

  .treat <- .attr(x, ".treat")

  if (chk::vld_string(model)) {
    model <- model |>
      match_arg(c("flat", "linear", "quadratic", "cubic")) |>
      switch(flat = ~1,
             linear = as.formula(sprintf("~ %s", paste(.treat, collapse = " + "))),
             quadratic = as.formula(sprintf("~ (%s)^2 + %s",
                                            paste(.treat, collapse = " + "),
                                            paste(sprintf("I(%s^2)", .treat), collapse = " + "))))
  }
  else if (!rlang::is_formula(model, lhs = FALSE)) {
    .err("{.arg model} must be a string or a one-sided formula with the projection model on the right-hand side")
  }
  else if (!all(get_varnames(model) %in% .treat)) {
    #Check that only treat is in model
    .err("only the treatment variables ({.and {.var {(.treat)}}}) are allowed to appear in {.arg model}")
  }

  .est <- .attr(x, ".est")
  .values <- .attr(x, ".values")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")
  .vcov_type <- .attr(x, ".vcov_type")

  proj_data <- model.frame(model, data = .values)

  mt <- .attr(proj_data, "terms")

  mm <- model.matrix(mt, data = proj_data)

  nm <- colnames(mm)

  if (!all(is.finite(mm))) {
    .err("evaluation of the projection model produced non-finite values of {.or {.var {(.treat)}}}, which is not allowed")
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

plot.effect_curve <- function(x, conf_level = 0.95, simultaneous = TRUE, null = NULL,
                              subset = NULL, proj = NULL, transform = TRUE,
                              ci.type = "perc", df = NULL, ...) {
  check_proj(x, proj)

  .est0 <- .attr(x, ".est")
  .vcov <- .attr(x, ".vcov")
  .vcov_type <- .attr(x, ".vcov_type")
  .treat <- .attr(x, ".treat")
  .values <- .attr(x, ".values")
  .boot <- .attr(x, ".boot")
  .draws <- .attr(x, ".draws")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")

  if (.vcov_type == "none") {
    conf_level <- 0
    transform <- FALSE
    simultaneous <- FALSE
  }
  else {
    conf_level <- process_conf_level(conf_level)
  }

  null <- process_null(null, x)

  transform_list <- process_transform(transform, x, .est0)

  # subset
  .subset <- process_subset_by_grid(substitute(subset),
                                    .by_grid = .by_grid,
                                    .contrast = .contrast)

  if (is_null(.by_grid) || is_not_null(.contrast)) {
    .s <- seq_along(.est0)
  }
  else {
    .s <- which(rep(.subset, each = length(.est0) / nrow(.by_grid)))
    .est0 <- .est0[.s]
  }

  # First, apply transforms if necessary
  transform <- transform_list$transform
  inv_transform <- transform_list$inv_transform

  .est <- transform(.est0)

  do_inference <- (.vcov_type != "none") &&
    is_not_null(.vcov) && (conf_level > 0)

  if (do_inference) {
    simultaneous <- process_simultaneous(simultaneous, .est)

    ci.type <- process_ci.type(ci.type, .vcov_type, simultaneous)

    res_names <- c(.treat, "estimate",
                   "conf.low", "conf.high")

    res <- make_df(res_names, length(.est))

    do_transform <- ci.type != "perc" &&
      !all(check_if_zero(.est - .est0)) &&
      all(check_if_zero(inv_transform(.est) - .est0))

    t_crit <- vp <- NULL

    # Compute CI
    if (conf_level > 0) {
      if (ci.type == "wald") {
        # Process df
        df <- df %or% .attr(x, ".df")

        chk::chk_number(df)
        chk::chk_gt(df, 0)

        stat <- if (is.finite(df)) "t" else "z"

        if (.vcov_type == "bootstrap") {
          .vcov <- .boot[["t"]] |>
            ss(j = .s) |>
            transform() |>
            cov()
        }
        else if (.vcov_type == "posterior") {
          .vcov <- .draws |>
            ss(j = .s) |>
            transform() |>
            cov()
        }
        else if (do_transform) {
          d_transform <- transform_list$d_transform
          .vcov <- quad_mult(diag(d_transform(.est0)), ss(.vcov, .s, .s))
        }
        else {
          .vcov <- ss(.vcov, .s, .s)
        }

        se <- sqrt(diag(.vcov))
        zeros <- se < 1e-10

        if (simultaneous) {
          vp <- ss(.vcov, !zeros, !zeros) |>
            cov2cor() |>
            .to_psd()

          t_crit <- try(mvtnorm::qmvt(conf_level,
                                      tail = "both.tails",
                                      df = df, corr = vp,
                                      keepAttr = FALSE,
                                      maxiter = 1e5,
                                      abseps = 1e-5)$quantile,
                        silent = TRUE)

          if (null_or_error(t_crit)) {
            .err("there was an error computing simultaneous confidence intervals")
          }
        }
        else {
          t_crit <- switch(stat,
                           "z" = abs(qnorm((1 - conf_level) / 2)),
                           "t" = abs(qt((1 - conf_level) / 2, df)))
        }
        # Reverse transformation
        res$conf.low[] <- inv_transform(.est - t_crit * se)
        res$conf.high[] <- inv_transform(.est + t_crit * se)
      }
      else if (.vcov_type == "bootstrap") {
        boot_ci <- confint(.boot, level = conf_level,
                           parm = .s,
                           ci.type = ci.type,
                           simultaneous = simultaneous)

        res$conf.low[] <- boot_ci[, 1L]
        res$conf.high[] <- boot_ci[, 2L]
      }
      else if (.vcov_type == "posterior") {
        draw_ci <- posterior_ci(.draws, level = conf_level,
                                parm = .s,
                                simultaneous = simultaneous)

        res$conf.low[] <- draw_ci[, 1L]
        res$conf.high[] <- draw_ci[, 2L]
      }
    }

    if (is_not_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(conf_level, "crit") <- t_crit
    }
  }
  else {
    res_names <- c(.treat, "estimate")
    res <- make_df(res_names, length(.est0))
  }

  # Use original estimate in output
  res$estimate <- inv_transform(.est)

  if (is_not_null(proj)) {
    res$proj_estimate <- .attr(proj, ".est")[.s]

    proj_lt <- "55"
  }

  # Apply by_grid and contrast labels if necessary
  if (is_not_null(.contrast)) {
    contrast <- factor(rep(.contrast, each = nrow(.values)),
                       levels = .contrast)

    res[.treat] <- ss(.values, rep.int(seq_row(.values), length(.contrast)), .treat)

    res <- cbind(contrast = contrast, res)

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$contrast)) +
      labs(color = "Contrast")

    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$contrast),
                    alpha = .2) +
        labs(fill = "Contrast")
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$contrast),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }

    p <- p +
      scale_color_brewer(palette = "Set1", aesthetics = c("color", "fill"),
                         breaks = .contrast,
                         limits = .contrast)
  }
  else if (is_not_null(.by_grid)) {
    res_names <- names(res)
    val_id <- seq_row(.values)
    by_id <- seq_row(.by_grid)

    res$.vi <- rep.int(val_id, sum(.subset))
    res$.by_id <- rep(which(.subset), each = length(val_id))

    res <- merge(res,
                 cbind(.by_grid, .by_id = by_id),
                 by = ".by_id",
                 all.x = TRUE, all.y = FALSE,
                 sort = FALSE) |>
      roworderv(c(".by_id", ".vi"))

    res[.treat] <- .values[res$.vi, .treat]

    .labels <- get_by_grid_labels(.by_grid, sep = ", ")

    res$.by_id <- factor(res$.by_id, levels = seq_row(.by_grid),
                         labels = .labels)

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$.by_id)) +
      labs(color = "Subgroup")


    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$.by_id),
                    alpha = .2) +
        labs(fill = "Subgroup")
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$.by_id),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }

    p <- p + scale_color_brewer(palette = "Set1", aesthetics = c("color", "fill"),
                                breaks = .labels[.subset],
                                limits = .labels)
  }
  else {
    res[.treat] <- .values[, .treat]

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate))

    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high),
                    alpha = .2)
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }
  }

  if (!anyNA(null)) {
    p <- p + geom_hline(yintercept = null, color = "black")
  }

  p + labs(x = .treat, y = .get_y_axis(x),
           title = .get_title(x),
           subtitle = .get_subtitle(x, conf_level, simultaneous)) +
    coord_cartesian(expand = c(left = FALSE, right = FALSE),
                    default = TRUE) +
    theme_bw() +
    theme(legend.position = "right",
          legend.direction = "vertical",
          geom = element_geom(ink = "#E41A1C"),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
}
