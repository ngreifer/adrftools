#' Estimate an average dose-response function (ADRF)
#'
#' @description
#' Estimates the average dose-response function (ADRF) for a fitted model object.
#'
#' @param x a fitted model object (e.g., from [lm()] or [glm()]).
#' @param treat a string specifying the name of the treatment variable.
#' @param vcov how the covariance matrix of the estimates should be computed. If `"unconditional"` (default), use the sandwich estimator including sampling uncertainty. If `"boot"` or `"fwb"`, use the traditional or fractional weighted bootstrap, respectively (both of which require the \CRANpkg{fwb} package to be installed). Otherwise, may be a covariance matrix or other allowed input to the `vcov` argument of [marginaleffects::get_vcov()]. Can also be `"none"` to avoid computing the uncertainty.
#' @param cluster an optional data frame or one-sided formula with the clustering terms for cluster-robust inference.
#' @param range numeric; a numeric vector corresponding either to the lower and upper bounds of the treatment values for which to compute the affect curve or a single number corresponding to the middle quantile of the treatment. Default is .95 to use the .025 and .975 quantiles of the treatment. See Details.
#' @param n integer specifying the number of equally spaced grid points on which to compute the effect curve anchor points. Default is 70; higher numbers increase computation time and size of the resulting object but improve accuracy.
#' @param type character string indicating the type of prediction. Passed to [marginaleffects::get_predict()]. Default is `"response"` for predictions on the scale of the outcome variable. Other options might include `"link"` for the linear predictor. This argument is ignored for `lm` objects.
#' @param data an optional data frame containing the observations originally used to fit the outcome model supplied to `x`. This should only be used if the supplied model is not supported by \pkg{insight}. In most cases, this should not need to be supplied.
#' @param subset an optional logical expression indicating the subset of data to use for estimation. Will be evaluated in the environment of the original dataset supplied to the model fitting function.
#' @param by optional variable(s) over which to group the estimation. Can be a character vector or one-sided formula.
#' @param wts optional numeric vector of weights to generalize the effect curve to a weighted population.
#' @param fwb.args an option list of arguments to be passed to \pkgfun{fwb}{fwb} when `vcov` is `"boot"` or `"fwb"`.
#' @param ... further arguments passed to [marginaleffects::get_predict()].
#'
#' @returns
#' An object of class `effect_curve`. This object is a function with attributes. See [`effect_curve-class`] for details on this function and its outputs.
#'
#' @details
#' `adrf()` estimates the ADRF by computing average predicted outcomes in the sample for counterfactual treatment values, optionally stratified by grouping variables and accounting for estimation uncertainty via unconditional or conditional variance estimation or bootstrapping. Unconditional variance estimation and bootstrapping treat the sample as random. When `vcov = "unconditional"`, the variance is computed using the formula in Hansen et al. (2024), which involves augmenting the influence function with a term to account for sampling from the superpopulation. Unconditional variance estimation requires [sandwich::estfun()] and [sandwich::bread()] methods for the supplied object to be available.
#'
#' When a `mira` object from \pkg{mice} or a `mimira` object from \pkg{MatchThem} is supplied, analyses are applied to each imputed dataset and pooled using Rubin's rules. Bootstrapping is not allowed with such objects.
#'
#' When a `svyglm` object from \pkg{survey} is supplied, `adrf()` automatically incorporates the survey weights extracted from the object. The same is true for `glm_weightit` objects, etc., from \pkg{WeightIt} when `s.weights` are supplied in the original call to `weightit()`. See `vignette("adrftools")` for more details on using the `wts` argument.
#'
#' ## `range`
#'
#' The `range` argument controls for which range of the treatment the effect curve is to be evaluated. It can be supplied either as two numbers corresponding to the lower and upper bounds for the treatment (e.g., `range = c(0, 10)`) or as a single number corresponding to the middle quantile of the treatment (e.g., `range = .9`, which uses the .05 and .95 quantiles of the treatment as the bounds). The default is .95 to use the .025 and .975 quantiles of the treatment. When supplied as a quantile, the quantiles are evaluated incorporating the weights supplied to `wts`.
#'
#' A reason not to use the full treatment range (e.g., by setting `range = 1`) is that there is likely very little certainty about the effect curve at the treatment extremes. This uncertainty can muddy tests of the effect curve. However, limiting the treatment range means inferences about the effect curve are less generalizable to more extreme values of the treatment. Note that this does not change the data used to fit the effect curve, just the points along the effect curve for which inference and estimation are to take place.
#'
#' @seealso
#' * [plot.effect_curve()] for plotting the ADRF
#' * [summary.effect_curve()] for testing hypotheses about the ADRF
#' * [`effect_curve`] for computing point estimates along the ADRF
#' * [curve_projection()] for projecting a simpler model onto the ADRF
#' * [reference_curve()] for computing the difference between each point on the ADRF and a specific reference point
#' * [curve_contrast()] for contrasting ADRFs computed within subgroups
#' * [amef()] for computing the average marginal effect function (AMEF), the derivative of the ADRF
#' * [marginaleffects::avg_predictions()] for computing average adjusted predictions for fitted models (similar to the ADRF)
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             Male * (Age + Race + PIR +
#'                       Enough_Food),
#'           data = nhanes3lead)
#'
#' # ADRF of logBLL on Math, unconditional
#' # inference
#' adrf1 <- adrf(fit, treat = "logBLL")
#'
#' adrf1
#'
#' ## Plot the ADRF
#' plot(adrf1)
#'
#' ## ADRF estimates at given points
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary()
#'
#' # ADRF of logBLL on Math, unconditional
#' # inference; manual range
#' adrf2 <- adrf(fit, treat = "logBLL",
#'               range = c(0, 2))
#'
#' adrf2
#'
#' plot(adrf2)
#'
#' # ADRF of logBLL on Math, bootstrap
#' # inference
#' \dontrun{
#' adrf_b <- adrf(fit, treat = "logBLL",
#'               vcov = "fwb")
#'
#' adrf_b
#'
#' plot(adrf_b)
#' }
#'
#' # ADRF in subset
#' adrf_m <- adrf(fit, treat = "logBLL",
#'                subset = Male == 1)
#'
#' adrf_m
#'
#' # ADRFs in subgroups
#' adrf_by <- adrf(fit, treat = "logBLL",
#'                 by = ~Male)
#'
#' adrf_by

#' @export
adrf <- function(x, ...) {
  UseMethod("adrf")
}

#' @exportS3Method adrf default
#' @rdname adrf
adrf.default <- function(x, treat, vcov = "unconditional", cluster = NULL, type = "response",
                         data = NULL, subset = NULL, by = NULL, wts = NULL,
                         range = .95, n = 71, fwb.args = list(), ...) {

  arg_not_missing(treat)

  arg_string(type)

  subset_substitute <- substitute(subset)

  switch(.get_model_type(x, type = type),
         "lm" = {
           estimator <- .get_estimator.lm(x)
           grad_fun <- .get_grad_fun.lm(x)
           pred_fun <- .get_pred_fun.lm(x)
         },
         "glm" = {
           estimator <- .get_estimator.glm(x)
           grad_fun <- .get_grad_fun.glm(x)
           pred_fun <- .get_pred_fun.glm(x, type = type)
         },
         "default" = {
           estimator <- .get_estimator.default(x)
           grad_fun <- .get_grad_fun.default(x)
           pred_fun <- .get_pred_fun.default(x, type = type)
         })

  .effect_curve_internal(
    x = x, treat = treat, vcov = vcov,
    range = range, n = n,
    data = data, subset = subset_substitute, by = by, wts = wts,
    cluster = cluster, fwb.args = fwb.args,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    ...
  )
}

#' @exportS3Method adrf mira
adrf.mira <- function(x, treat, vcov = "unconditional", cluster = NULL, type = "response",
                      data = NULL, subset = NULL, by = NULL, wts = NULL,
                      range = .95, n = 71, ...) {

  arg_not_missing(treat)

  arg_string(type)

  subset_substitute <- substitute(subset)

  x1 <- x[["analyses"]][[1L]]

  switch(.get_model_type(x1, type = type),
         "lm" = {
           estimator <- .get_estimator.lm(x1)
           grad_fun <- .get_grad_fun.lm(x1)
           pred_fun <- .get_pred_fun.lm(x1)
         },
         "glm" = {
           estimator <- .get_estimator.glm(x1)
           grad_fun <- .get_grad_fun.glm(x1)
           pred_fun <- .get_pred_fun.glm(x1, type = type)
         },
         "default" = {
           estimator <- .get_estimator.default(x1)
           grad_fun <- .get_grad_fun.default(x1)
           pred_fun <- .get_pred_fun.default(x1, type = type)
         })

  .effect_curve_internal_mi(
    x = x, treat = treat, vcov = vcov,
    range = range, n = n,
    data = data, subset = subset_substitute, by = by, wts = wts,
    cluster = cluster,
    estimator = estimator,
    grad_fun = grad_fun,
    pred_fun = pred_fun,
    ...
  )
}

.effect_curve_internal <- function(x, treat, vcov, range = .95, n,
                                   data, subset, by, wts, cluster = NULL, fwb.args = list(),
                                   estimator = NULL, grad_fun = NULL,
                                   pred_fun = NULL,
                                   .adrf_env = parent.frame(2L), ...) {

  # get data
  model_data <- process_model_data(x, data)

  # treat
  treat_var <- process_treat(treat, model_data)

  # wts
  wts <- process_wts(wts, x, data)

  # values
  values <- process_range(range, n, treat_var, w = wts)
  n <- length(values)

  # bayes
  is_bayes <- .is_bayesian(x, model_data)

  # vcov + cluster
  cl_list <- process_vcov_and_cluster(vcov, x, cluster,
                                      is_bayes, model_data)

  vcov <- cl_list[["vcov"]]
  cl_list[["vcov"]] <- NULL

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

  # Speed up if treat is the only variable by treating dataset as a single observation
  if (allv(names(model_data), treat) && all_the_same(by_id) && !isFALSE(...get(".quick"))) {
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

  data_grid <- ss(model_data, rep.int(ii, n))

  data_grid[[treat]] <- values[val_id]

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

    arg_list(fwb.args)

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

  .make_fun(.est = drop(unname(.est)),
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

.effect_curve_internal_mi <- function(x, treat, vcov, range = .95, n,
                                      data, subset, by, wts, cluster = NULL,
                                      estimator = NULL, grad_fun = NULL,
                                      pred_fun = NULL,
                                      .adrf_env = parent.frame(2L), ...) {
  # get data
  model_data.complete <- process_model_data_mi(x, data)

  treat_var.complete <- process_treat(treat, model_data.complete)

  imp_split <- gsplit(seq_row(model_data.complete), model_data.complete[[".imp"]])

  # wts
  wts.complete <- process_wts_mi(wts, x, data)

  #values
  values <- process_range(range, n, treat_var.complete,
                          w = wts.complete)

  n <- length(values)
  vi <- seq_len(n)

  # by
  processed_by_list <- process_by(by, model_data.complete)
  by_grid <- processed_by_list$by_grid
  by_id.complete <- processed_by_list$by_id

  # bayes
  are_bayes <- vapply(seq_along(imp_split), function(.i) {
    .is_bayesian(x[["analyses"]][[.i]], ss(model_data.complete, imp_split[[.i]]))
  }, logical(1L))

  if (any(are_bayes) && !all(are_bayes)) {
    .err("either all models fit to the imputed datasets must be Bayesian or no models must be")
  }

  is_bayes <- any(are_bayes)

  # family
  .family <- process_family(x[["analyses"]][[1L]])

  if (any_apply(x[["analyses"]][-1L], function(mod) {
    !identical(.family, process_family(mod),
               ignore.environment = TRUE)
  })) {
    .err("all models must have the same family")
  }

  # df
  .df <- process_df(x)

  # response
  .response <- process_response(x[["analyses"]][[1L]])

  # vcov + cluster
  cl_list.list <- process_vcov_and_cluster_mi(vcov, x[["analyses"]], cluster, is_bayes,
                                              model_data.complete, imp_split)

  m <- length(x[["analyses"]])

  est.list <- vcov.list <- make_list(m)

  .draws <- NULL

  vcov_type <- {
    if (identical(cl_list.list[[1L]][["vcov"]], "none")) "none"
    else if (is_bayes) "posterior"
    else if (identical(cl_list.list[[1L]][["vcov"]], "unconditional")) "unconditional"
    else "conditional"
  }

  for (.i in seq_len(m)) {
    .imp <- imp_split[[.i]]

    model_data <- ss(model_data.complete, .imp)

    # treat
    treat_var <- treat_var.complete[.imp]

    # wts
    wts <- wts.complete[.imp]

    # vcov + cluster
    cl_list <- cl_list.list[[.i]]

    vcov <- cl_list[["vcov"]]

    cl_list[["vcov"]] <- NULL

    # subset
    s <- process_subset(subset, model_data, env = .adrf_env)
    model_data <- ss(model_data, s)

    # by
    by_id <- by_id.complete[.imp]

    # model
    xi <- x[["analyses"]][[.i]]

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

    data_grid <- ss(model_data, rep.int(ii, n))

    data_grid[[treat]] <- values[val_id]

    if (!is_bayes && is_null(pred_fun)) {
      beta0 <- marginaleffects::get_coef(xi)

      aliased <- is.na(beta0)

      mm <- get_tt(xi, model_data) |>
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

    if (is_bayes) {
      .draws <- .draws |>
        rbind(rsplit(seq_row(data_grid), rep.int(val_id, nrow(data_grid) / (ni * n))) |>
                lapply(function(.j) {
                  est_fun(xi, data0 = ss(data_grid, .j)) |>
                    .attr("posterior_draws") |>
                    .fmean(g = by_id, w = ss(wi, .j))
                }) |>
                do.call(what = "rbind") |>
                t())
    }
    else if (identical(vcov_type, "none")) {
      est <- est_fun(xi) |>
        .fmean(g = cross_id, w = wi)
    }
    else if (identical(vcov_type, "unconditional")) {
      beta0 <- beta0[!aliased]

      phat <- est_fun(xi)

      est <- .fmean(phat, g = cross_id, w = wi)

      gr <- {
        if (is_null(grad_fun) || is_not_null(pred_fun)) {
          .gradient(function(b) {
            est_fun(xi, b) |>
              .fmean(g = cross_id, w = wi) |>
              drop()
          }, .x = beta0)
        }
        else {
          grad_fun(mm, eta = drop(mm %*% beta0)) |>
            .fmean(g = cross_id, w = wi)
        }
      }

      if (inherits(xi, "svyglm")) {
        infl <- ni * (.attr(xi, "influence") %or% tcrossprod(model.matrix(xi) * resid(xi, "working") * xi$weights,
                                                             xi$naive.cov))
      }
      else {
        infl <- sandwich::estfun(xi) |>
          tcrossprod(sandwich::bread(xi))
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

      if (is_null(cl_list)) {
        vcov.list[[.i]] <- crossprod(S / nrow(S))
      }
      else {
        vcov.list[[.i]] <- Reduce("+", lapply(seq_along(cl_list[["cluster"]]), function(j) {
          crossprod(fsum(S, cl_list[["cluster"]][[j]]) / nrow(S)) * cl_list[["adj"]][j]
        }))
      }
    }
    else {
      est <- est_fun(xi) |>
        .fmean(g = cross_id, w = wi)

      V <- marginaleffects::get_vcov(xi, vcov = vcov)

      gr <- {
        if (is_null(grad_fun) || is_not_null(pred_fun)) {
          .gradient(function(b) {
            est_fun(xi, b) |>
              .fmean(g = cross_id, w = wi) |>
              drop()
          }, .x = beta0)
        }
        else {
          grad_fun(mm, eta = drop(mm %*% beta0)) |>
            .fmean(g = cross_id, w = wi)
        }
      }

      vcov.list[[.i]] <- quad_mult(gr, V)
    }

    est.list[[.i]] <- est
  }

  # Pool estimates
  if (is_bayes) {
    pooled <- list(est = colMeans(.draws))

    if (identical(vcov_type, "none")) {
      .draws <- NULL
    }
    else {
      pooled$vcov <- stats::cov(.draws)
    }
  }
  else {
    pooled <- pool_est(est.list, vcov.list)
  }

  attr(vcov_type, "clustered") <- is_not_null(cl_list.list[[1L]][["cluster"]])

  .make_fun(.est = unname(pooled$est),
            .vcov = unname(pooled$vcov),
            .values = values,
            .treat = treat,
            .vcov_type = vcov_type,
            .draws = .draws,
            .curve_type = "ADRF",
            .family = .family,
            .df = .df,
            .response = .response,
            .by_grid = by_grid)
}
