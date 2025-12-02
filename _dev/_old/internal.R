.est_internal <- function(.est_type = "adrf", x, treat, vcov, values, n, subset_substitute, by, wts, eps = NULL,
                          estimator = NULL, grad_fun = NULL,
                          pred_fun = NULL, mcall, env = parent.frame(2L), ...) {

  model_data <- process_data(x)

  # treat
  treat_var <- process_treat(treat, model_data)

  # values
  values <- process_values(values, n, treat_var)
  n <- length(values)

  # subset
  s <- process_subset(subset_substitute, model_data, env = env)
  model_data <- ss(model_data, s)

  # wts
  wts <- process_wts(wts, x)

  # by
  processed_by_list <- process_by(by, model_data)
  by_grid <- processed_by_list$by_grid
  by_id <- processed_by_list$by_id

  # vcov
  vcov <- process_vcov(vcov, x)

  # est
  beta0 <- marginaleffects::get_coef(x)

  # Speed up if treat is the only variable by treating dataset as a single observation
  if (allv(names(model_data), treat) && allv(by_id, by_id[1L]) && !isFALSE(...get(".quick"))) {
    model_data <- ss(model_data, 1L)
    wts <- 1.0
    s <- 1L
    by_id <- by_id[1L]

    .fmean <- function(x, ...) {x}
  }
  else {
    .fmean <- collapse::fmean
  }

  ni <- nrow(model_data)
  ii <- seq_len(ni)

  vi <- seq_len(n)
  val_id <- rep(vi, each = ni)

  cross_id <- GRP(list(by_id = rep.int(by_id, n), val_id = val_id),
                  sort = TRUE)

  wi <- rep.int(ss(wts, s), n)

  if (.est_type == "adrf") {
    eps <- NULL

    compute_E_amef <- FALSE

    data_grid <- ss(model_data, rep.int(ii, n))

    data_grid[[treat]] <- values[val_id]

    .f <- function(z) {z}
  }
  else if (.est_type == "amef") {
    eps <- process_eps(eps, treat_var)

    compute_E_amef <- min(values) <= min(treat_var) && max(values) >= max(treat_var)

    data_grid <- ss(model_data, c(rep.int(ii, n), rep.int(ii, n)))

    ## iii: units in `a + eps/2` group
    iiip <- seq_len(n * length(ii))
    iiim <- iiip + length(iiip)
    data_grid[[treat]][iiip] <- values[val_id] + eps / 2
    data_grid[[treat]][iiim] <- values[val_id] - eps / 2

    E_est_list <- NULL

    .f <- function(z) {
      nr <- {
        if (is_null(dim(z))) length(z) / 2L
        else nrow(z) / 2L
      }

      (ss(z, seq_len(nr)) - ss(z, nr + seq_len(nr))) / eps
    }
  }
  else {
    .err("unrecognized `.est_type`")
  }

  if (is_null(pred_fun)) {
    mm <- delete.response(x$terms) |>
      model.matrix(data = data_grid)

    eta <- drop(mm %*% beta0)

    est_fun <- function(b, model) {
      estimator(mm, b, model)
    }
  }
  else {
    est_fun <- function(b, model) {
      pred_fun(beta = b, mod = model, data_grid = data_grid, ...)
    }
  }

  .vcov <- vcov_type <- S <- boot <- NULL

  if (identical(vcov, "none")) {
    est <- est_fun(beta0, x) |>
      .f() |>
      .fmean(g = cross_id, w = wi)
  }
  else if (identical(vcov, "boot") || identical(vcov, "fwb")) {
    bootfun <- function(data, w = alloc(1.0, ni), ...) {
      wts_boot <- wts * w
      x_boot <- do.call("update", list(x, weights = wts_boot))

      beta0_boot <- marginaleffects::get_coef(x_boot)

      wi_boot <- {
        if (ni == 1L) alloc(1.0, n)
        else rep.int(ss(wts_boot, s), n)
      }

      est_fun(beta0_boot, x_boot) |>
        .f() |>
        .fmean(g = cross_id, w = wi_boot)
    }

    if (identical(vcov, "boot")) {
      rlang::local_options(fwb_wtype = "multinom")
    }

    boot <- fwb::fwb(data = process_data(x),
                     statistic = bootfun, drop0 = FALSE, ...)

    est <- coef(boot)
    .vcov <- vcov(boot)
    vcov_type <- "bootstrap"
  }
  else if (identical(vcov, "unconditional")) {
    phat <- est_fun(beta0, x) |> .f()

    est <- .fmean(phat, g = cross_id, w = wi)

    gr <- {
      if (is_null(grad_fun) || is_not_null(pred_fun)) {
        .gradient(function(b, mod) {
          est_fun(b, mod) |>
            .f() |>
            .fmean(g = cross_id, w = wi) |>
            drop()
        }, .x = beta0, mod = x)
      }
      else {
        grad_fun(mm, eta, x) |>
          .f() |>
          fmean(g = cross_id, w = wi)
      }
    }

    if (inherits(x, "svyglm")) {
      infl <- {
        if (is_null(.attr(x, "influence")))
          ni * tcrossprod((model.matrix(x) * resid(x, "working") * x$weights),
                          x$naive.cov)
        else
          ni * .attr(x, "influence")
      }
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
        S[b, k] <- S[b, k] + (phat_split[[k]] - est[k]) * wts_sum / sum(wts[b])
      }
    }

    .vcov <- crossprod(S / nrow(S))

    vcov_type <- "unconditional"
  }
  else {
    phat <- est_fun(beta0, x) |> .f()

    est <- .fmean(phat, g = cross_id, w = wi)

    V <- marginaleffects::get_vcov(x, vcov = vcov)

    gr <- {
      if (is_null(grad_fun) || is_not_null(pred_fun)) {
        .gradient(function(b, mod) {
          est_fun(b, mod) |>
            .f() |>
            .fmean(g = cross_id, w = wi) |>
            drop()
        }, .x = beta0, mod = x)
      }
      else {
        grad_fun(mm, eta, x) |>
          .f() |>
          .fmean(g = cross_id, w = wi)
      }
    }

    .vcov <- quad_mult(gr, V)

    vcov_type <- "conditional"
  }

  out <- list(est = unname(drop(est)),
              values = values,
              treat = treat,
              vcov = unname(.vcov),
              vcov_type = vcov_type,
              boot = boot,
              S = S)

  # #Average AMEF
  # if (compute_E_amef) {
  #   E_est_list <- .compute_E_amef(treat_var, values, est, .vcov, S, by_id,
  #                                 wts, vcov_type)
  #
  #   out$E_est <- E_est_list$E_est
  #   out$V_est <- E_est_list$V_est
  # }

  out$by_grid <- by_grid
  out$call <- mcall

  attr(out, "model") <- x
  attr(out, "treat_var") <- treat_var
  attr(out, "by_id") <- processed_by_list$by_id

  class(out) <- switch(.est_type,
                       adrf = "adrf_est",
                       amef = "amef_est")

  out
}

.est_internal_mi <- function(.est_type = "adrf", x, treat, vcov, values, n, subset_substitute, by, wts, eps = NULL,
                             estimator = NULL, grad_fun = NULL,
                             pred_fun = NULL, mcall, env = parent.frame(2L), ...) {

  model_data.complete <- lapply(x[["analyses"]], process_data) |>
    rowbind(idcol = ".imp")

  treat_var.complete <- process_treat(treat, model_data.complete)

  imp_split <- gsplit(seq_row(model_data.complete), model_data.complete[[".imp"]])

  #values
  values <- process_values(values, n, treat_var.complete)

  n <- length(values)
  vi <- seq_len(n)

  # by
  processed_by_list <- process_by(by, model_data.complete)
  by_grid <- processed_by_list$by_grid
  by_id.complete <- processed_by_list$by_id

  vcov.input.list <- process_vcov_mi(vcov, x[["analyses"]])

  if (.est_type == "adrf") {
    eps <- NULL

    compute_E_amef <- FALSE

    .f <- function(z) {z}
  }
  else if (.est_type == "amef") {
    eps <- process_eps(eps, treat_var.complete)

    compute_E_amef <- min(values) <= min(treat_var.complete) && max(values) >= max(treat_var.complete)

    .f <- function(z) {
      nr <- {
        if (is_null(dim(z))) length(z) / 2L
        else nrow(z) / 2L
      }

      (ss(z, seq_len(nr)) - ss(z, nr + seq_len(nr))) / eps
    }
  }
  else {
    .err("unrecognized `.est_type`")
  }

  m <- length(x[["analyses"]])

  est.list <- vcov.list <- make_list(m)

  if (compute_E_amef) {
    E_amef.list <- V_amef.list <- make_list(m)
  }

  for (.i in seq_len(m)) {
    .imp <- imp_split[[.i]]

    model_data <- ss(model_data.complete, .imp)

    # treat
    treat_var <- treat_var.complete[.imp]

    # subset
    s <- process_subset(subset_substitute, model_data, env = env)
    model_data <- ss(model_data, s)

    # wts
    wts <- process_wts(wts, x[["analyses"]][[.i]])

    # by
    by_id <- by_id.complete[.imp]

    # vcov
    vcov <- vcov.input.list[[.i]]

    # est
    beta0 <- marginaleffects::get_coef(x[["analyses"]][[.i]])

    # Speed up if treat is the only variable by treating dataset as a single observation
    if (allv(names(model_data), treat) && allv(by_id, by_id[1L]) && !isFALSE(...get(".quick"))) {
      model_data <- ss(model_data, 1L)
      wts <- 1.0
      s <- 1L
      by_id <- by_id[1L]

      .fmean <- function(x, ...) {x}
    }
    else {
      .fmean <- collapse::fmean
    }

    ni <- nrow(model_data)
    ii <- seq_len(ni)

    val_id <- rep(vi, each = ni)

    cross_id <- GRP(list(by_id = rep.int(by_id, n), val_id = val_id), sort = TRUE)

    wi <- rep.int(ss(wts, s), n)

    if (.est_type == "adrf") {
      data_grid <- ss(model_data, rep.int(ii, n))

      data_grid[[treat]] <- values[val_id]
    }
    else if (.est_type == "amef") {
      data_grid <- ss(model_data, c(rep.int(ii, n), rep.int(ii, n)))

      ## iii: units in `a + eps/2` group
      iiip <- seq_len(n * length(ii))
      iiim <- iiip + length(iiip)
      data_grid[[treat]][iiip] <- values[val_id] + eps / 2
      data_grid[[treat]][iiim] <- values[val_id] - eps / 2

      E_est_list <- NULL
    }

    if (is_null(pred_fun)) {
      mm <- delete.response(x[["analyses"]][[.i]]$terms) |>
        model.matrix(data = data_grid)

      eta <- drop(mm %*% beta0)

      est_fun <- function(b, model) {
        estimator(mm, b, model)
      }
    }
    else {
      est_fun <- function(b, model) {
        pred_fun(b, model, data_grid)
      }
    }

    vcov_type <- S <- boot <- NULL

    if (identical(vcov, "none")) {
      est <- est_fun(beta0, x[["analyses"]][[.i]]) |>
        .f() |>
        .fmean(g = cross_id, w = wi)
    }
    else if (identical(vcov, "unconditional")) {
      phat <- est_fun(beta0, x[["analyses"]][[.i]]) |> .f()

      est <- .fmean(phat, g = cross_id, w = wi)

      gr <- {
        if (is_null(grad_fun) || is_not_null(pred_fun)) {
          .gradient(function(b, mod) {
            est_fun(b, mod) |>
              .f() |>
              .fmean(g = cross_id, w = wi) |>
              drop()
          }, .x = beta0, mod = x[["analyses"]][[.i]])
        }
        else {
          grad_fun(mm, eta, x[["analyses"]][[.i]]) |>
            .f() |>
            .fmean(g = cross_id, w = wi)
        }
      }

      S <- sandwich::estfun(x[["analyses"]][[.i]]) |>
        tcrossprod(sandwich::bread(x[["analyses"]][[.i]])) |>
        tcrossprod(gr)

      if (ni > 1L) {
        wts_sum <- sum(wts)

        phat_split <- gsplit(phat, cross_id)

        b_split <- gsplit(ii, by_id)

        for (k in seq_along(phat_split)) {
          b <- s[b_split[[cross_id$groups$by_id[k]]]]
          S[b, k] <- S[b, k] + (phat_split[[k]] - est[k]) * wts_sum / sum(wts[b])
        }
      }

      vcov.list[[.i]] <- crossprod(S / nrow(S))

      vcov_type <- "unconditional"
    }
    else {
      phat <- est_fun(beta0, x[["analyses"]][[.i]]) |> .f()

      est <- .fmean(phat, g = cross_id, w = wi)

      V <- marginaleffects::get_vcov(x[["analyses"]][[.i]], vcov = vcov)

      gr <- {
        if (is_null(grad_fun) || is_not_null(pred_fun)) {
          .gradient(function(b, mod) {
            est_fun(b, mod) |>
              .f() |>
              .fmean(g = cross_id, w = wi) |>
              drop()
          }, .x = beta0, mod = x[["analyses"]][[.i]])
        }
        else {
          grad_fun(mm, eta, x[["analyses"]][[.i]]) |>
            .f() |>
            .fmean(g = cross_id, w = wi)
        }
      }

      vcov.list[[.i]] <- quad_mult(gr, V)

      vcov_type <- "conditional"
    }

    est.list[[.i]] <- est

    #Average AMEF
    # if (compute_E_amef) {
    #   E_est_list <- .compute_E_amef(treat_var, values, est, vcov.list[[.i]], S, by_id,
    #                                 wts, vcov_type)
    #
    #   E_amef.list[[.i]] <- E_est_list$E_est
    #   V_amef.list[[.i]] <- E_est_list$V_est
    # }
  }

  # Pool estimates
  pooled <- pool_est(est.list, vcov.list)

  out <- list(est = unname(pooled$est),
              values = values,
              treat = treat,
              vcov = unname(pooled$vcov),
              vcov_type = vcov_type)

  # if (compute_E_amef) {
  #   E_amef_pooled <- pool_est(E_amef.list, V_amef.list)
  #
  #   out$E_est <- E_amef_pooled$est
  #   out$V_est <- E_amef_pooled$vcov
  # }

  out$by_grid <- by_grid
  out$call <- mcall

  attr(out, "model") <- x
  attr(out, "treat_var") <- treat_var.complete
  attr(out, "by_id") <- by_id.complete

  class(out) <- switch(.est_type,
                       adrf = "adrf_est",
                       amef = "amef_est")

  out
}

.summary_internal <- function(object, conf_level = 0.95, null = NULL, simultaneous = TRUE,
                              transform_list = NULL, df = Inf, ci.type = "perc", values = NULL) {

  vcov_type <- object[["vcov_type"]]

  # Only do inference if vcov_type is not empty and conf_level is
  # nonzero or null hypothesis is specified
  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  is_curve <- inherits(object, c("adrf_est", "amef_est"))

  term_col <- if (is_curve) object$treat else "term"

  if (is_null(values)) {
    est0 <- coef(object)

    term_val <- if (is_curve) object$values else object$coef_names

    if (do_inference) {
      if (vcov_type != "bootstrap") {
        vcov0 <- stats::vcov(object)
      }
    }
  }
  else {
    n <- length(values)

    val_mat <- matrix(0, nrow = n, ncol = length(object[["values"]]))

    matches <- max.col(-abs(outer(values, object[["values"]], "-")), "first")

    to_interp <- !check_if_zero(values - object[["values"]][matches])

    if (any(to_interp)) {
      val_mat[to_interp,] <- get_locpoly_w(x = values[to_interp],
                                           v = object[["values"]],
                                           constant = 0.5)
    }

    if (!all(to_interp)) {
      val_mat[cbind(which(!to_interp), matches[!to_interp])] <- 1
    }

    # est0 is the estimates at values, including interp'd
    est0 <- drop(val_mat %*% coef(object))

    term_val <- values

    if (do_inference) {
      if (vcov_type == "bootstrap") {
        object$boot[["t"]] <- tcrossprod(object$boot[["t"]][, seq_along(coef(object)), drop = FALSE],
                                         val_mat)
        object$boot[["t0"]] <- tcrossprod(object$boot[["t0"]][seq_along(coef(object))], val_mat) |>
          drop()
      }
      else {
        vcov0 <- quad_mult(val_mat, stats::vcov(object))
      }
    }
  }

  if (do_inference) {
    stat <- if (is.finite(df)) "t" else "z"

    if (is_null(null)) {
      res_names <- c(term_col,
                     "estimate", "std.error",
                     "conf.low", "conf.high")
    }
    else {
      res_names <- c(term_col,
                     "estimate", "std.error",
                     stat, "p.value",
                     "conf.low", "conf.high")
    }

    res <- make_df(res_names, length(est0))

    # First, apply transforms if necessary
    if (is_null(transform_list)) {
      transform_list <- process_transform()
    }

    transform <- transform_list$transform
    inv_transform <- transform_list$inv_transform

    est <- transform(est0)

    do_transform <- (vcov_type != "bootstrap" || ci.type != "perc") &&
      !all(check_if_zero(est - est0)) && all(check_if_zero(inv_transform(est) - est0))

    if (do_transform) {
      d_transform <- transform_list$d_transform

      if (is_null(inv_transform)) {
        inv_transform <- make_inverse(transform, est0)
      }

      if (is_not_null(null)) {
        null <- transform(null)
      }

      if (vcov_type == "bootstrap") {
        vcov <- cov(transform(object$boot[["t"]][, seq_along(est), drop = FALSE]))
      }
      else {
        vcov <- quad_mult(diag(d_transform(est0)), vcov0)
      }
    }
    else {
      if (vcov_type == "bootstrap") {
        vcov <- cov(object$boot[["t"]][, seq_along(est), drop = FALSE])
      }
      else {
        vcov <- vcov0
      }
    }

    se <- sqrt(diag(vcov))
    zeros <- se < 1e-10

    # Compute test statistic on transformed estimates if null is specified
    if (is_not_null(null)) {
      res[[stat]][!zeros] <- (est[!zeros] - null) / se[!zeros]
    }

    t_crit <- t_dist <- vp <- NULL

    # Compute CI
    if (conf_level > 0) {
      if (vcov_type == "bootstrap" && ci.type != "wald") {
        boot_ci <- confint(object[["boot"]], level = conf_level,
                           parm = seq_along(est),
                           ci.type = ci.type,
                           simultaneous = simultaneous)

        res$conf.low[] <- boot_ci[, 1L]
        res$conf.high[] <- boot_ci[, 2L]
      }
      else if (simultaneous) {
        vp <- ss(vcov, !zeros, !zeros) |>
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

        # Reverse transformation
        res$conf.low[] <- inv_transform(est - t_crit * se)
        res$conf.high[] <- inv_transform(est + t_crit * se)

      }
      else {
        t_crit <- switch(stat,
                         "z" = abs(qnorm((1 - conf_level) / 2)),
                         "t" = abs(qt((1 - conf_level) / 2, df)))

        # Reverse transformation
        res$conf.low[] <- inv_transform(est - t_crit * se)
        res$conf.high[] <- inv_transform(est + t_crit * se)
      }
    }

    # Compute p-values
    if (is_not_null(null)) {
      if (vcov_type == "bootstrap" && ci.type != "wald") {
        s <- summary(object[["boot"]],
                     conf = 0,
                     parm = seq_along(est),
                     ci.type = ci.type,
                     p.value = TRUE,
                     null = null,
                     simultaneous = simultaneous)

        res$p.value <- s[, ncol(s)]
      }
      else if (simultaneous) {
        if (is_null(vp)) {
          vp <- ss(vcov, !zeros, !zeros) |>
            cov2cor() |>
            .to_psd()
        }

        pvals <- try(vapply(abs(res[[stat]][!zeros]),
                            function(zi) {
                              1 - mvtnorm::pmvt(lower = rep.int(-zi, sum(!zeros)),
                                                upper = rep.int(zi, sum(!zeros)),
                                                df = df, corr = vp,
                                                keepAttr = FALSE,
                                                abseps = 1e-5)
                            }, numeric(1L)),
                     silent = TRUE)

        if (null_or_error(pvals)) {
          .err("there was an error computing simultaneous p-values")
        }

        res$p.value[!zeros] <- pvals
      }
      else {
        res$p.value[!zeros] <- switch(stat,
                                      "z" = 2 * pnorm(-abs(res[[stat]][!zeros])),
                                      "t" = 2 * pt(-abs(res[[stat]][!zeros]), df))
      }
    }

    # Remove SE if transformation is not identity
    res$std.error <- {
      if (do_transform) NULL
      else se
    }

    # Use original estimate in output
    res$estimate <- {
      if (do_transform) inv_transform(est)
      else est0
    }
  }
  else {
    res_names <- c(term_col, "estimate")
    res <- make_df(res_names, length(est0))

    # Use original estimate in output
    res$estimate <- est0
  }

  # Apply by_grid and contrast labels if necessary
  if (is_null(object$by_grid) && is_null(object$contrast)) {
    res[[term_col]] <- term_val
    attr(res, "treat") <- object$treat
  }
  else if (is_not_null(object$by_grid)) {
    res_names <- names(res)
    val_id <- seq_along(term_val)
    by_id <- seq_row(object$by_grid)

    res$.vi <- rep.int(val_id, length(by_id))
    res$.by_id <- rep(by_id, each = length(val_id))

    res <- merge(res,
                 cbind(object$by_grid, .by_id = by_id),
                 by = ".by_id",
                 all.x = TRUE, all.y = FALSE,
                 sort = FALSE)

    res <- ss(res, order(res$.by_id, res$.vi))

    res[[term_col]] <- term_val[res$.vi]

    by_names <- names(object$by_grid)
    res <- res[c(by_names, res_names)]

    attr(res, "treat") <- object$treat
    attr(res, "by") <- by_names
  }
  else {
    contrast <- factor(rep(object$contrast, each = length(term_val)),
                       levels = object$contrast)

    res[[term_col]] <- rep.int(term_val, length(object$contrast))

    res <- cbind(contrast = contrast, res)

    attr(res, "treat") <- object$treat
    attr(res, "contrast") <- object$contrast
  }

  attr(res, "obj") <- object
  attr(res, "values") <- values

  if (do_inference) {
    attr(res, "simultaneous") <- simultaneous
    attr(res, "conf_level") <- conf_level

    if (conf_level > 0 && is_not_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(res, "crit") <- t_crit
    }

    if (vcov_type == "bootstrap") {
      attr(res, "ci.type") <- ci.type
    }

    attr(res, "null") <- null
  }

  res
}

.print_summary_internal <- function(x, digits, topn, ...) {
  is_contrast <- inherits(.attr(x, "obj"), "curve_contrast")

  summary_type <- {
    if (inherits(x, "point_contrast")) "point"
    else if (inherits(x, "basic_curve_summary")) "raw"
    else if (inherits(x, "summary.curve_proj")) "proj"
    else "estimate"
  }

  curve_obj <- {
    if (inherits(x, "summary.curve_proj")) .attr(x, "obj") |> .attr("obj")
    else .attr(x, "obj")
  }

  qoi <- {
    if (is_contrast) {
      if (inherits(curve_obj, "adrf_est")) "ADRF Difference"
      else "AMEF Difference"
    }
    else if (inherits(curve_obj, "adrf_est")) "ADRF"
    else if (inherits(curve_obj, "amef_est")) "AMEF"
  }

  main <- sprintf("%s %s",
                  qoi,
                  switch(summary_type,
                         point = "Pointwise Contrasts",
                         proj = "Projection Coefficients",
                         "Estimates"))

  if (utils::hasName(x, "comp.vals")) {
    x[["comp.vals"]] <- NULL
  }

  tmp <- .print_estimate_table(x = x, digits = digits, topn = topn,
                               help.fn = switch(summary_type,
                                                point = "point_contrast",
                                                raw = "print.adrf_est",
                                                proj = "summary.curve_proj",
                                                estimate = "summary.adrf_est"),
                               main = main,
                               rownames = FALSE,
                               ...)

  cat(tmp, sep = "\n")

  vcov_type <- .attr(x, "obj")[["vcov_type"]]

  if (is_not_null(vcov_type)) {
    #Additional info
    if (is_not_null(.attr(x, "crit"))) {
      df <- .attr(x, "crit") |> .attr("df")

      stat_msg <- {
        if (is.finite(df)) {
          sprintf(" (t* = %s, df = %s)",
                  round(.attr(x, "crit"), 3),
                  round(df, 1))
        }
        else {
          sprintf(" (z* = %s)",
                  round(.attr(x, "crit"), 3))
        }
      }
    }
    else {
      stat_msg <- ""
    }

    simultaneous <- .attr(x, "simultaneous")

    sprintf("Inference: %s%s\nConfidence level: %s%%%s\n",
            vcov_type,
            if (isTRUE(simultaneous)) ", simultaneous"
            else if (isFALSE(simultaneous)) ", pointwise"
            else "",
            round(100 * .attr(x, "conf_level"), 2),
            stat_msg) |>
      .it() |>
      cat()
  }

  computed_values <- .attr(x, "values")

  if (is_not_null(computed_values)) {
    original_values <- .attr(x, "obj")[["values"]]

    matches <- max.col(-abs(outer(computed_values, original_values, "-")), "first")

    to_interp <- !check_if_zero(computed_values - original_values[matches])

    if (any(to_interp)) {
      sprintf("Interpolated values: %s \u2208 {%s}\n",
              .attr(x, "treat"),
              toString(computed_values[to_interp])) |>
        .it() |>
        cat()
    }
  }

  invisible(x)
}

.print_internal <- function(x, digits, topn, ...) {
  res <- summary(x, conf_level = 0, transform = FALSE, null = NULL)

  class(res) <- c(class(res), "basic_curve_summary")

  print(res, digits = digits, topn = topn, ...)

  invisible(x)
}

.plot_internal <- function(x, ylab = "E[Y(a)]", null = NULL, proj = NULL) {
  treat <- .attr(x, "treat")
  by <- .attr(x, "by")
  contrast <- .attr(x, "contrast")

  ci <- is_not_null(x$conf.low) && is_not_null(x$conf.high)

  if (is_not_null(proj)) {
    x$proj_estimate <- proj$fitted.values

    inv_transform <- .attr(proj, "inv_transform")
    if (is_not_null(inv_transform) && is.function(inv_transform)) {
      x$proj_estimate <- inv_transform(x$proj_estimate)
    }

    proj_lt <- "55"
  }

  if (is_null(by) && is_null(contrast)) {
    p <- ggplot(x, aes(x = .data[[treat]])) +
      geom_line(aes(y = .data$estimate))

    if (ci) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high),
                    alpha = .2)
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate),
                         linetype = proj_lt)
    }
  }
  else if (is_not_null(by)) {
    by_data <- x[by]
    by_grid <- unique(by_data)
    by_grid$.by_id <- seq_row(by_grid)

    x <- merge(x, by_grid, by = by,
               all.x = TRUE, all.y = FALSE,
               sort = FALSE)

    labels <- do.call(function(...) paste(..., sep = ", "),
                      lapply(by, function(i) {
                        sprintf("%s = %s", i, by_grid[[i]])
                      }))

    x$.by_id <- factor(x$.by_id, levels = seq_row(by_grid),
                       labels = labels)

    p <- ggplot(x, aes(x = .data[[treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$.by_id))
    if (ci) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$.by_id),
                    alpha = .2)
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$.by_id),
                         linetype = proj_lt)
    }
  }
  else {
    p <- ggplot(x, aes(x = .data[[treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$contrast))
    if (ci) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$contrast),
                    alpha = .2)
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$contrast),
                         linetype = proj_lt)
    }
  }

  if (is_not_null(null)) {
    p <- p + geom_hline(yintercept = null, color = "black")
  }

  p + labs(x = treat, y = ylab, fill = NULL, color = NULL) +
    coord_cartesian(expand = c(left = FALSE, right = FALSE),
                    default = TRUE) +
    scale_color_brewer(palette = "Set1", aesthetics = c("color", "fill")) +
    theme_bw() +
    theme(legend.position = "bottom",
          geom = element_geom(ink = "#E41A1C"))
}
