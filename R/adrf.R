#' Estimate an average dose-response function (ADRF)
#'
#' @description
#' Estimates the average dose-response function (ADRF) for a fitted model object.
#'
#' @param x a fitted model object (e.g., from [lm()] or [glm()]).
#' @param treat a string specifying the name of the treatment variable.
#' @param vcov how the covariance matrix of the estimates should be computed. If `"unconditional"` (the default for frequentist models), use the sandwich estimator including sampling uncertainty, with the HC2 small-sample adjustment when it is available (see Details). If `"boot"` or `"fwb"`, use the traditional or fractional weighted bootstrap, respectively (both of which require the \CRANpkg{fwb} package to be installed). Otherwise, may be a covariance matrix or other allowed input to the `vcov` argument of [marginaleffects::get_vcov()]. Can also be `"none"` to avoid computing the uncertainty. For Bayesian models, only `"posterior"`, which uses the posterior of the estimates, and `"none"` are allowed. For models fit to multiply imputed data, `"boot"` and `"fwb"` are not allowed.
#' @param cluster an optional data frame or one-sided formula with the clustering terms for cluster-robust inference. When `vcov` is `"unconditional"`, the CR2 small-sample adjustment is used when it is available (see Details). When `vcov` is `"conditional"`, the covariance matrix of the model coefficients is computed using [sandwich::vcovCL()] with its default type; `vcov` can also be `"HC0"`, `"HC1"`, `"HC2"`, or `"HC3"` to request that type from [sandwich::vcovCL()] (`"HC2"` and `"HC3"` are available only for `lm()` and `glm()` models). Other strings and functions supplied to `vcov` cannot be used with `cluster`.
#' @param range numeric; a numeric vector corresponding either to the lower and upper bounds of the treatment values for which to compute the effect curve or a single number corresponding to the middle quantile of the treatment. Default is .95 to use the .025 and .975 quantiles of the treatment. See Details.
#' @param n integer specifying the number of equally spaced grid points on which to compute the effect curve anchor points. Must be between 2 and 1000. Default is 51; higher numbers increase computation time and size of the resulting object but improve accuracy.
#' @param type character string indicating the type of prediction. Passed to [marginaleffects::get_predict()]. Default is `"response"` for predictions on the scale of the outcome variable. Other options might include `"link"` for the linear predictor. This argument is ignored for `lm` objects.
#' @param data an optional data frame containing the observations originally used to fit the outcome model supplied to `x`. This is needed when the dataset cannot be recovered from the model (e.g., when the model was fit inside a function or is not supported by \pkg{insight}), and it is used to find variables named in `subset`, `by`, `wts`, and `cluster` that are not in the dataset used to fit the model. In most cases, this does not need to be supplied.
#' @param subset an optional logical expression indicating the subset of data to use for estimation. Will be evaluated in the environment of the original dataset supplied to the model fitting function. Can also be a vector of row numbers or row names; logical vectors and row numbers refer to the units used to fit the model (i.e., excluding any dropped because of missing values), and row names refer to the original dataset.
#' @param by optional variable(s) over which to group the estimation. Can be a character vector or one-sided formula.
#' @param wts optional weights used to generalize the effect curve to a weighted population. Can be a numeric vector, the name of a variable in the original dataset, `TRUE` to use the weights used to fit the model, or `FALSE` for no weights. The default depends on the model; see Details.
#' @param fwb.args an optional list of arguments to be passed to \pkgfun{fwb}{fwb} when `vcov` is `"boot"` or `"fwb"`.
#' @param ... further arguments passed to [marginaleffects::get_predict()].
#'
#' @returns
#' An object of class `<effect_curve>`. This object is a function with attributes. See [`effect_curve-class`] for details on this function and its outputs.
#'
#' @details
#' `adrf()` estimates the ADRF by computing average predicted outcomes in the sample for counterfactual treatment values, optionally stratified by grouping variables and accounting for estimation uncertainty via unconditional or conditional variance estimation or bootstrapping. Unconditional variance estimation and bootstrapping treat the sample as random. When `vcov = "unconditional"`, the variance is computed using the formulas in Hansen et al. (2024), which involves augmenting the influence function with a term to account for sampling from the superpopulation. Unconditional variance estimation requires [sandwich::estfun()] and [sandwich::bread()] methods for the supplied object to be available.
#'
#' In the unconditional variance, each unit's contribution to the influence function of the model coefficients is divided by \eqn{\sqrt{1 - h_i}}, where \eqn{h_i} is the unit's leverage (i.e., hat value), as in the HC2 estimator of MacKinnon and White (1985). Without this adjustment, the variance would share the downward bias of the HC0 sandwich estimator, which is largest for estimates that depend on units with high leverage, such as the ADRF at the extremes of the treatment, and which can make the omnibus tests in [summary.effect_curve()] reject a true null hypothesis too often in small samples. The leverages are computed from the fitted model for `lm()` and `glm()` models and extracted using [hatvalues()] for other models; when they cannot be found, the HC1 adjustment is used instead with a message, which multiplies the variance contributed by the model coefficients by \eqn{n/(n-p)}, where \eqn{n} is the number of units and \eqn{p} the number of estimated parameters, and corrects the bias only for the average leverage. The adjustment is not used with `svyglm` models, whose variance is computed from the survey design.
#'
#' When `cluster` is supplied, the unconditional variance uses the CR2 adjustment described by Pustejovsky and Tipton (2018), the cluster-level analog of HC2: within each cluster, the vector of residuals is multiplied by \eqn{(I - H_{gg})^{-1/2}}, where \eqn{H_{gg}} is the block of the hat matrix for the units in the cluster, in place of the usual \eqn{G/(G-1)} adjustment, where \eqn{G} is the number of clusters. With clusters of one unit, it is the HC2 adjustment. CR2 is available only for `lm()` and `glm()` models (including those whose influence function is that of an ordinary `glm()` fit, such as `glm_weightit()` models with fixed weights); for other models, the CR1 adjustment, \eqn{\frac{G}{G-1}\frac{n-1}{n-p}} (the default of `sandwich::vcovCL()`), is used with a message. For generalized linear models, the residuals and hat matrix are those of the final iteratively reweighted least squares fit, so the adjusted variance can differ from that of `sandwich::vcovCL(type = "HC2")` when the working weights vary within clusters.
#'
#' When a `<mira>` object from \pkg{mice} or a `<mimira>` object from \pkg{MatchThem} is supplied, analyses are applied to each imputed dataset and pooled using Rubin's rules. The degrees of freedom used for tests and confidence intervals are those of Reiter (2007), or of Li, Raghunathan, and Rubin (1991) when the complete-data degrees of freedom are infinite (e.g., for generalized linear models). These account for the number of imputations and the fraction of missing information, which is assumed to be the same for all points on the curve, so a single value is used for every point and for joint inference. Bootstrapping is not allowed with such objects.
#'
#' ## `range`
#'
#' The `range` argument controls for which range of the treatment the effect curve is to be evaluated. It can be supplied either as two numbers corresponding to the lower and upper bounds for the treatment (e.g., `range = c(0, 10)`) or as a single number corresponding to the middle quantile of the treatment (e.g., `range = .9`, which uses the .05 and .95 quantiles of the treatment as the bounds). The default is .95 to use the .025 and .975 quantiles of the treatment. When supplied as a quantile, the quantiles are evaluated incorporating the weights supplied to `wts`.
#'
#' A reason not to use the full treatment range (e.g., by setting `range = 1`) is that there is likely very little certainty about the effect curve at the treatment extremes. This uncertainty can muddy tests of the effect curve. However, limiting the treatment range means inferences about the effect curve are less generalizable to more extreme values of the treatment. Note that this does not change the data used to fit the effect curve, just the points along the effect curve for which inference and estimation are to take place.
#'
#' ## `wts`
#'
#' The `wts` argument supplies weights that define the population over which the predicted outcomes are averaged, such as sampling weights. When supplied as a numeric vector, it must have one value for each unit used to fit the model (or, with multiply imputed data, either one value for each unit in each imputed dataset or one value for each row of the imputed datasets stacked together). When supplied as a string, it is the name of a variable in the original dataset. `wts = TRUE` uses the weights used to fit the model, as extracted by [insight::get_weights()], and `wts = FALSE` weights all units equally.
#'
#' When `wts` is not supplied, the default depends on the model. When a `<svyglm>` object from \pkg{survey} is supplied, the survey weights are used. When a `<glm_weightit>` object, etc., from \pkg{WeightIt} is supplied, the sampling weights are used if `s.weights` were supplied in the original call to `weightit()`; the balancing weights are not, because they are used to fit the outcome model rather than to define the target population. For other models, all units are weighted equally, even when the model was fit with weights. See `vignette("adrftools")` for more details on using the `wts` argument.
#'
#' @references
#' Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an *F* reference distribution. *Journal of the American Statistical Association*, 86(416), 1065–1073. \doi{10.1080/01621459.1991.10475152}
#'
#' MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305–325. \doi{10.1016/0304-4076(85)90158-7}
#'
#' Pustejovsky, J. E., & Tipton, E. (2018). Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *Journal of Business & Economic Statistics*, 36(4), 672–683. \doi{10.1080/07350015.2016.1247004}
#'
#' Reiter, J. P. (2007). Small-sample degrees of freedom for multi-component significance tests with multiple imputation for missing data. *Biometrika*, 94(2), 502–508. \doi{10.1093/biomet/asm028}
#'
#' @seealso
#' * [plot.effect_curve()] for plotting the ADRF
#' * [summary.effect_curve()] for testing hypotheses about the ADRF
#' * [`effect_curve`] for computing point estimates along the ADRF
#' * [curve_projection()] for projecting a simpler model onto the ADRF
#' * [reference_curve()] for computing the difference between each point on the ADRF and a specific reference point
#' * [curve_contrast()] for contrasting ADRFs computed within subgroups
#' * [amef()] for computing the average marginal effect function (AMEF), the derivative of the ADRF
#' * [shift_curve()] for estimating the mean outcome when every unit's treatment is shifted from its observed value
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
#' @examplesIf rlang::is_installed("fwb")
#' # ADRF of logBLL on Math, bootstrap
#' # inference
#' \donttest{adrf_b <- adrf(fit, treat = "logBLL",
#'                vcov = "fwb")
#'
#' adrf_b
#'
#' plot(adrf_b)}
#' @examples
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
                         range = .95, n = 51, fwb.args = list(), ...) {

  arg::arg_supplied(treat)

  .check_misplaced_args(...names(), "adrf")

  subset_substitute <- substitute(subset)

  .effect_curve_internal(
    x = x, treat = treat, vcov = vcov,
    range = range, n = n,
    data = data, subset = subset_substitute,
    by = by, wts = wts, cluster = cluster, fwb.args = fwb.args,
    type = type, ...
  )
}

#' @exportS3Method adrf mira
adrf.mira <- function(x, treat, vcov = "unconditional", cluster = NULL, type = "response",
                      data = NULL, subset = NULL, by = NULL, wts = NULL,
                      range = .95, n = 51, ...) {

  arg::arg_supplied(treat)

  .check_misplaced_args(...names(), "adrf")

  subset_substitute <- substitute(subset)

  .effect_curve_internal(
    x = x, treat = treat, vcov = vcov,
    range = range, n = n,
    data = data, subset = subset_substitute,
    by = by, wts = wts, cluster = cluster,
    type = type, ...
  )
}

# Estimates the ADRF for a single fitted model or for the models fit to each
# imputed dataset of a `mira` object. A single model is treated as a list of one
# model. Data, weights, the covariance specification, and `subset` are
# processed for each model, the subgroups defined by `by` are formed across all
# models, and .adrf_fit_one() computes the estimates for each model. Estimates
# from multiply imputed data are pooled using Rubin's rules, or, for Bayesian
# models, by stacking the posterior draws.
#
# With `shift_fun` (from process_shift_fun()), the curve is a shift curve: the
# counterfactual treatment of each unit is `shift_fun(A, delta)` for its
# observed treatment `A` and each grid value `delta` of `shift`, rather than
# the grid value itself.
.effect_curve_internal <- function(x, treat, vcov, range, n,
                                   data, subset, by, wts, cluster = NULL, fwb.args = list(),
                                   type, .adrf_env = parent.frame(2L),
                                   shift_fun = NULL, shift = NULL, ...) {
  # type
  arg::arg_string(type)

  is_mi <- inherits(x, "mira")

  models <- {
    if (is_mi) x[["analyses"]]
    else list(x)
  }

  m <- length(models)

  # get data: the predictors (model_data) and all variables in the original
  # data (eval_data), for subset, by, wts, and cluster
  if (is_mi) {
    model_data.complete <- process_model_data_mi(x, data)

    imp_split <- gsplit(seq_row(model_data.complete), model_data.complete[[".imp"]])

    eval_data.complete <- process_eval_data_mi(x, model_data.complete, data)

    pred_names <- setdiff(names(model_data.complete), ".imp")

    model_data.list <- lapply(imp_split, function(.i) {
      ss(model_data.complete, .i, pred_names)
    })

    eval_data.list <- lapply(imp_split, function(.i) {
      ss(eval_data.complete, .i)
    })
  }
  else {
    model_data.list <- list(process_model_data(x, data))

    eval_data.list <- list(process_eval_data(x, model_data.list[[1L]], data))
  }

  # Variables used in offsets supplied outside the formula are not predictors
  # but are needed to compute predictions
  model_data.list <- lapply(seq_len(m), function(.i) {
    add_offset_vars(models[[.i]], model_data.list[[.i]], eval_data.list[[.i]])
  })

  # treat
  treat_var <- lapply(model_data.list, process_treat, treat = treat) |>
    unlist()

  # wts
  wts.list <- process_wts_list(wts, models, eval_data.list)

  is_shift <- !is_null(shift_fun)

  # values; for a shift curve, the grid of shifts is set after `subset`
  if (!is_shift) {
    values <- process_range(range, n, treat_var, w = unlist(wts.list))
    n <- length(values)
  }

  # bayes
  are_bayes <- vapply(seq_len(m), function(.i) {
    .is_bayesian(models[[.i]], model_data.list[[.i]])
  }, logical(1L))

  is_bayes <- any(are_bayes)

  if (is_bayes && !all(are_bayes)) {
    arg::err("either all models fit to the imputed datasets must be Bayesian or no models must be")
  }

  # family
  .family <- process_family(models[[1L]], type)

  if (any_apply(models[-1L], function(mod) {
    !identical(.family, process_family(mod, type),
               ignore.environment = TRUE)
  })) {
    arg::err("all models must have the same family")
  }

  # vcov + cluster
  cl_list.list <- process_vcov_and_cluster_list(vcov, models, cluster,
                                                is_bayes, eval_data.list)

  vcov_type <- {
    if (identical(cl_list.list[[1L]][["vcov"]], "none")) "none"
    else if (is_bayes) "posterior"
    else if (any(cl_list.list[[1L]][["vcov"]] %in% c("boot", "fwb"))) "bootstrap"
    else if (identical(cl_list.list[[1L]][["vcov"]], "unconditional")) "unconditional"
    else "conditional"
  }

  # subset
  s.list <- lapply(seq_len(m), function(.i) {
    process_subset(subset, eval_data.list[[.i]], env = .adrf_env)
  })

  .shift <- NULL

  if (is_shift) {
    A.list <- lapply(seq_len(m), function(.i) {
      model_data.list[[.i]][[treat]][s.list[[.i]]]
    })

    .shift <- process_shift(shift, shift_fun, n = n, A.list = A.list,
                            w.list = lapply(seq_len(m), function(.i) {
                              wts.list[[.i]][s.list[[.i]]]
                            }),
                            treat = treat)

    values <- .shift[["values"]]
    n <- length(values)

    .shift[["values"]] <- NULL
  }

  # by, with subgroups defined across all models after subsetting
  processed_by_list <- lapply(seq_len(m), function(.i) {
    ss(eval_data.list[[.i]], s.list[[.i]])
  }) |>
    rowbind() |>
    process_by(by = by)

  by_grid <- processed_by_list$by_grid
  n_by <- get_n_by(.by_grid = by_grid)

  by_codes <- as.integer(processed_by_list$by_id) |>
    gsplit(rep.int(seq_len(m), vlengths(s.list)))

  if (any_apply(by_codes, function(b) !all(seq_len(n_by) %in% b))) {
    arg::err("every subgroup defined by {.arg by} must be present in every imputed dataset")
  }

  if (any_apply(seq_len(m), function(.i) {
    any(fsum(ss(wts.list[[.i]], s.list[[.i]]), g = by_codes[[.i]]) <= 0)
  })) {
    arg::err("the values of {.arg wts} must have a positive sum in every subgroup defined by {.arg subset} and {.arg by}")
  }

  # df
  .df <- process_df(x)

  # response
  .response <- process_response(models[[1L]])

  fits <- lapply(seq_len(m), function(.i) {
    .adrf_fit_one(x = models[[.i]], treat = treat, values = values,
                  model_data = model_data.list[[.i]], s = s.list[[.i]],
                  wts = wts.list[[.i]], by_id = by_codes[[.i]], n_by = n_by,
                  cl_list = cl_list.list[[.i]], vcov_type = vcov_type,
                  is_bayes = is_bayes, type = type, data = data,
                  fwb.args = fwb.args, shift_fun = .shift[["fun"]], ...)
  })

  .boot <- .draws <- .vcov <- NULL

  if (is_bayes) {
    # Posterior draws from all imputed datasets are stacked
    .draws <- lapply(fits, .subset2, "draws") |>
      do.call(what = "rbind")

    .est <- colMeans(.draws)

    if (identical(vcov_type, "none")) {
      .draws <- NULL
    }
    else {
      .vcov <- stats::cov(.draws)
    }
  }
  else if (m == 1L) {
    .est <- fits[[1L]][["est"]]
    .vcov <- fits[[1L]][["vcov"]]
    .boot <- fits[[1L]][["boot"]]
  }
  else {
    # Rubin's rules
    pooled <- pool_est(lapply(fits, .subset2, "est"),
                       if (vcov_type != "none") lapply(fits, .subset2, "vcov"))

    .est <- pooled[["est"]]
    .vcov <- pooled[["vcov"]]

    if (vcov_type != "none") {
      .df <- .mi_df(pooled[["vcov_w"]], pooled[["vcov_b"]], m = m, df_com = .df)
    }
  }

  attr(vcov_type, "clustered") <- !is_null(cl_list.list[[1L]][["cluster"]])

  # The shift function is not stored, so the object does not carry its
  # environment
  if (is_shift) {
    if (!all(is.finite(.est))) {
      arg::err("some estimates of the shift curve are not finite, possibly because some shifted treatment values are outside the values the model allows; try a narrower {.arg shift}")
    }

    .shift[["fun"]] <- NULL
  }

  # Density of the treatment in each subgroup, from all units used (not only
  # those within `range`), for weighting by the treatment distribution. A
  # shift curve is a function of the shift, which has no such distribution.
  .treat_density <- {
    if (is_shift) NULL
    else .get_treat_density(
      lapply(seq_len(m), function(.i) model_data.list[[.i]][[treat]][s.list[[.i]]]),
      lapply(seq_len(m), function(.i) wts.list[[.i]][s.list[[.i]]]),
      by_codes, n_by
    )
  }

  # `.treat` names the argument of the effect curve, which for a shift curve
  # is the shift; the treatment is then stored in `.shift`
  .make_fun(.est = drop(unname(.est)),
            .vcov = unname(.vcov),
            .values = values,
            .treat = if (is_shift) "shift" else treat,
            .shift = .shift,
            .vcov_type = vcov_type,
            .boot = .boot,
            .draws = .draws,
            .curve_type = "ADRF",
            .family = .family,
            .df = .df,
            .response = .response,
            .by_grid = by_grid,
            .treat_density = .treat_density)
}

# Computes the ADRF estimates for one fitted model: the point estimates and,
# depending on `vcov_type`, their covariance matrix, bootstrap replicates, or
# posterior draws. `by_id` contains integer subgroup codes for the units in
# `subset`, and `n_by` is the number of subgroups.
.adrf_fit_one <- function(x, treat, values, model_data, s, wts, by_id, n_by,
                          cl_list, vcov_type, is_bayes, type, data = NULL,
                          fwb.args = list(), shift_fun = NULL, ...) {
  vcov <- cl_list[["vcov"]]
  infl <- cl_list[["infl"]]
  cl_list[c("vcov", "infl")] <- NULL

  model_data <- ss(model_data, s)

  by_id <- qG(by_id)

  n <- length(values)

  # type is ignored for linear models
  if (inherits(x, "lm") && !inherits(x, "glm")) {
    type <- "response"
  }

  # The counterfactual treatment of every unit at the k-th grid value: the
  # grid value itself for an ADRF, or the shifted observed treatment for a
  # shift curve
  .cf_at <- {
    if (is_null(shift_fun)) {
      function(k) rep.int(values[k], fnrow(model_data))
    }
    else {
      A_obs <- model_data[[treat]]

      function(k) shift_fun(A_obs, values[k])
    }
  }

  # Speed up if treat is the only variable by treating dataset as a single
  # observation; not for a shift curve, where the counterfactual treatment
  # differs across units
  if (is_null(shift_fun) && allv(names(model_data), treat) && all_the_same(by_id) &&
      !isFALSE(...get(".quick"))) {
    model_data <- ss(model_data, 1L)
    wts <- 1.0
    s <- 1L
    by_id <- by_id[1L]

    .fmean <- function(x, ...) {x}
  }
  else {
    .fmean <- collapse::fmean
  }

  ni <- fnrow(model_data)
  ii <- seq_len(ni)

  vi <- seq_len(n)
  val_id <- rep(vi, each = ni)

  cross_id <- GRP(list(by_id = rep.int(by_id, n), val_id = val_id),
                  sort = TRUE)

  w_unit <- ss(wts, s)
  wi <- rep.int(w_unit, n)

  # The data with the treatment set to values[k] for every unit, and the grid
  # of every unit at every treatment value, with the rows ordered by treatment
  # value, then unit
  .data_at <- function(k) {
    d <- model_data
    d[[treat]] <- .cf_at(k)
    d
  }

  .data_grid <- function() {
    data_grid <- ss(model_data, rep.int(ii, n))
    data_grid[[treat]] <- {
      if (is_null(shift_fun)) values[val_id]
      else unlist(lapply(vi, .cf_at), use.names = FALSE)
    }
    data_grid
  }

  # A sample of up to 20 rows of the grid for checking the model type
  rows <- seq(1, ni * n, length.out = min(ni * n, 20L)) |>
    round() |>
    funique()

  rows_i <- (rows - 1L) %% ni + 1L
  rows_k <- (rows - 1L) %/% ni + 1L

  data_check <- ss(model_data, rows_i)
  data_check[[treat]] <- vapply(seq_along(rows), function(r) {
    .cf_at(rows_k[r])[rows_i[r]]
  }, numeric(1L))

  switch(.get_model_type(x, type = type, data_grid = data_check,
                         model_data = model_data, is_bayes = is_bayes, ...),
         "lm" = {
           linkinv <- .get_linkinv.lm(x)
           grad_fun <- .get_grad_fun.lm(x)
           pred_fun <- .get_pred_fun.lm(x)
         },
         "glm" = {
           linkinv <- .get_linkinv.glm(x)
           grad_fun <- .get_grad_fun.glm(x)
           pred_fun <- .get_pred_fun.glm(x, type = type)
         },
         {
           linkinv <- .get_linkinv.default(x)
           grad_fun <- .get_grad_fun.default(x)
           pred_fun <- .get_pred_fun.default(x, type = type)
         })

  if (!is_bayes) {
    # Aliased coefficients are dropped once, here, for every branch;
    # marginaleffects::get_vcov() drops them from the covariance matrix too
    beta0 <- marginaleffects::get_coef(x)
    aliased <- is.na(beta0)
    beta0 <- beta0[!aliased]
  }

  # Predictions from the model matrix and coefficients ("lm" and "glm" model
  # types) are computed one treatment value at a time, so that only one model
  # matrix the size of the data is held at once rather than one for the whole
  # grid
  fast <- !is_bayes && is_null(pred_fun)

  if (fast) {
    tt <- get_tt(x, model_data)

    # The variables of the model (the model frame) are evaluated on the whole
    # grid at once, as predict() would with the grid as new data, because a
    # variable can depend on other rows (e.g., `I(treat - mean(treat))` or
    # `rank(treat)`) and so would differ if evaluated one treatment value at
    # a time. The model matrix, which is computed row by row from the model
    # frame and is usually much larger, is then built for one treatment value
    # at a time from the rows of the model frame for that value.
    data_grid <- .data_grid()

    # Offsets (which model.matrix() leaves out) for each row of the grid, or 0
    # when the model has no offset
    off_grid <- .grid_offset(x, data_grid, model_data)

    mf <- stats::model.frame(tt, data = data_grid, xlev = .model_xlev(x))

    rm(data_grid)

    # nrow() rather than fnrow(), which counts the rows of the first column,
    # a matrix for terms like poly()
    if (nrow(mf) != ni * n) {
      if (is_null(shift_fun)) {
        arg::err("the variables in the model could not be computed for every unit at every treatment value")
      }
      else {
        arg::err("the variables in the model could not be computed for every unit at every shift, possibly because some shifted treatment values are outside the values the model allows; try a narrower {.arg shift}")
      }
    }

    .rows_at <- function(k) {
      (k - 1L) * ni + ii
    }

    .off_at <- function(k) {
      if (length(off_grid) == 1L) off_grid
      else off_grid[.rows_at(k)]
    }

    .mm_at <- function(k) {
      mm <- model.matrix(tt, data = mf[.rows_at(k), , drop = FALSE])

      # Subsetting copies the matrix, so only when needed
      if (any(aliased)) mm[, !aliased, drop = FALSE]
      else mm
    }

    # Coefficients of a refit model (e.g., in a bootstrap replicate) use the
    # original aliasing pattern so they match the columns of the model matrix;
    # a coefficient newly aliased in the refit is set to 0, as predict() does
    .refit_coef <- function(model) {
      b <- marginaleffects::get_coef(model)[!aliased]
      b[is.na(b)] <- 0
      b
    }

    # Predictions in the order of the grid and, if `grad = TRUE`, the gradient
    # of the weighted mean predictions in each subgroup and at each treatment
    # value with respect to the coefficients, in the order of `cross_id`
    # (subgroup, then treatment value). The gradient is computed as weighted
    # sums so that the rows of the model matrix are not multiplied in a copy.
    .fast_pass <- function(b, grad = FALSE) {
      phat <- vector("list", n)
      gr <- if (grad) vector("list", n)

      for (k in vi) {
        mm <- .mm_at(k)

        eta <- drop(mm %*% b) + .off_at(k)

        phat[[k]] <- linkinv(eta)

        if (grad) {
          gr[[k]] <- fsum(mm, g = by_id, w = w_unit * grad_fun(eta)) /
            fsum(w_unit, g = by_id)
        }
      }

      list(phat = unlist(phat, use.names = FALSE),
           gr = if (grad) {
             gr <- do.call("rbind", gr)
             ss(gr, .value_major_to_by_major(n, fnrow(gr) / n))
           })
    }

    est_fun <- function(model = NULL, b = NULL) {
      if (is_null(b)) {
        b <- .refit_coef(model)
      }

      .fast_pass(b)[["phat"]]
    }
  }
  else {
    # The Bayesian branch predicts one treatment value at a time with
    # .data_at(), as before
    data_grid <- if (!is_bayes) .data_grid()

    est_fun <- function(model = NULL, b = NULL) {
      pred_fun(data_grid = data_grid, beta0 = b, model = model, ...)
    }
  }

  out <- list(est = NULL, vcov = NULL, boot = NULL, draws = NULL)

  if (is_bayes) {
    # Predictions are computed one treatment value at a time to limit the size
    # of the posterior draws held at once
    out[["draws"]] <- lapply(vi, function(k) {
      pred_fun(data_grid = .data_at(k), model = x, ...) |>
        .attr("posterior_draws") |>
        .fmean(g = by_id, w = w_unit)
    }) |>
      do.call(what = "rbind") |>
      t() |>
      ss(j = .value_major_to_by_major(n, n_by))

    if (fnrow(out[["draws"]]) < 2L) {
      arg::err("at least two posterior draws are required")
    }

    if (anyNA(out[["draws"]])) {
      arg::err("the posterior draws of the predictions contain missing values")
    }

    out[["est"]] <- colMeans(out[["draws"]])
  }
  else if (identical(vcov_type, "none")) {
    out[["est"]] <- est_fun(x) |>
      .fmean(g = cross_id, w = wi)
  }
  else if (identical(vcov_type, "bootstrap")) {
    # The model is refit with its own weights times the bootstrap weights; the
    # target weights times the bootstrap weights are used to average the
    # predictions
    refit <- .make_boot_refit(x, data)

    if (fast) {
      # Predictions are needed in every replicate, so the model matrix is
      # stored rather than rebuilt. Columns that are identical at every
      # treatment value are stored once and the others once per treatment
      # value. Which columns are identical is found by comparing them rather
      # than from the formula, because a variable can depend on the treatment
      # without naming it, so the stored columns reproduce the model matrix at
      # each treatment value exactly.
      mm1 <- .mm_at(1L)

      tc <- rep.int(FALSE, fncol(mm1))

      for (k in vi[-1L]) {
        tc <- tc | .cols_differ(.mm_at(k), mm1)
      }

      mm_fixed <- mm1[, !tc, drop = FALSE]
      mm_treat <- lapply(vi, function(k) {
        (if (k == 1L) mm1 else .mm_at(k))[, tc, drop = FALSE]
      })

      rm(mm1)

      est_fun <- function(model = NULL, b = NULL) {
        if (is_null(b)) {
          b <- .refit_coef(model)
        }

        eta_fixed <- drop(mm_fixed %*% b[!tc])
        b_treat <- b[tc]

        lapply(vi, function(k) {
          linkinv(eta_fixed + drop(mm_treat[[k]] %*% b_treat) + .off_at(k))
        }) |>
          unlist(use.names = FALSE)
      }
    }

    bootfun <- function(data, w, ...) {
      x_boot <- withCallingHandlers(refit(w), warning = function(cnd) {
        if (grepl("non-integer", conditionMessage(cnd), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      })

      wi_boot <- {
        if (ni == 1L) alloc(1.0, n)
        else rep.int(ss(wts * w, s), n)
      }

      est_fun(x_boot) |>
        .fmean(g = cross_id, w = wi_boot)
    }

    if (is_null(fwb.args)) {
      fwb.args <- list()
    }

    arg::arg_list(fwb.args)

    if (identical(vcov, "boot")) {
      fwb.args[["wtype"]] <- "multinom"
    }

    if (!is_null(cl_list[["cluster"]])) {
      # Only one clustering variable is allowed with bootstrapping, which is
      # checked in process_vcov_and_cluster()
      fwb.args[["cluster"]] <- cl_list[["cluster"]][[1L]]
    }

    out[["boot"]] <- do.call(fwb::fwb, c(list(data = process_model_data(x, data),
                                              statistic = bootfun,
                                              drop0 = FALSE),
                                         fwb.args))

    out[["est"]] <- stats::coef(out[["boot"]])
    out[["vcov"]] <- stats::vcov(out[["boot"]])
  }
  else if (identical(vcov_type, "unconditional")) {
    if (fast) {
      fp <- .fast_pass(beta0, grad = TRUE)
      phat <- fp[["phat"]]
      gr <- fp[["gr"]]
    }
    else {
      phat <- est_fun(x)

      gr <- .gradient(function(b) {
        est_fun(x, b) |>
          .fmean(g = cross_id, w = wi) |>
          drop()
      }, .x = beta0)
    }

    out[["est"]] <- .fmean(phat, g = cross_id, w = wi)

    # `infl` (the influence of each unit on the coefficients) was computed by
    # .coef_influence() in process_vcov_and_cluster()

    if (!inherits(x, "svyglm")) {
      # bread() is scaled by the number of units with nonzero model weights,
      # but `S` is averaged over every row of estfun(), including units with
      # zero model weights, so the scales are matched
      n_bread <- try(stats::nobs(x), silent = TRUE)

      if (!null_or_error(n_bread) && is.numeric(n_bread) && length(n_bread) == 1L &&
          n_bread > 0 && n_bread < fnrow(infl)) {
        infl <- infl * fnrow(infl) / n_bread

        if (!is_null(cl_list[["infl_cl"]])) {
          cl_list[["infl_cl"]] <- lapply(cl_list[["infl_cl"]], `*`, fnrow(infl) / n_bread)
        }
      }
    }

    S <- tcrossprod(infl, gr)

    if (ni > 1L) {
      # Influence of sampling the covariates on the weighted mean of the
      # predictions in subgroup b, scaled like the rows of `S`
      N <- fnrow(S)

      # `phat` is ordered by treatment value, then unit, so its columns as a
      # matrix are the treatment values; the columns of `S` for subgroup g are
      # its n treatment values (in the order of `cross_id`)
      P <- matrix(phat, nrow = ni)

      b_split <- gsplit(ii, by_id)

      for (g in seq_len(n_by)) {
        u <- b_split[[g]]
        b <- s[u]
        k <- (g - 1L) * n + vi

        S[b, k] <- S[b, k] + wts[b] * TRA(ss(P, u), out[["est"]][k], "-") * N / sum(wts[b])
      }
    }

    out[["vcov"]] <- {
      if (inherits(x, "svyglm")) .svyglm_vcov(x, S)
      else if (is_null(cl_list)) crossprod(S / fnrow(S))
      else {
        V <- Reduce("+", lapply(seq_along(cl_list[["cluster"]]), function(j) {
          # With the CR2 adjustment, each clustering term has its own adjusted
          # influence on the coefficients
          S_j <- {
            if (is_null(cl_list[["infl_cl"]])) S
            else S + tcrossprod(cl_list[["infl_cl"]][[j]] - infl, gr)
          }

          crossprod(fsum(S_j, cl_list[["cluster"]][[j]]) / fnrow(S)) * cl_list[["adj"]][j]
        }))

        # Multiway clustering can give negative eigenvalues
        if (length(cl_list[["cluster"]]) > 1L) .clip_psd(V)
        else V
      }
    }
  }
  else {
    if (fast) {
      fp <- .fast_pass(beta0, grad = TRUE)
      phat <- fp[["phat"]]
      gr <- fp[["gr"]]
    }
    else {
      phat <- est_fun(x)

      gr <- .gradient(function(b) {
        est_fun(x, b) |>
          .fmean(g = cross_id, w = wi) |>
          drop()
      }, .x = beta0)
    }

    out[["est"]] <- .fmean(phat, g = cross_id, w = wi)

    V <- marginaleffects::get_vcov(x, vcov = vcov)

    if (fncol(gr) != fnrow(V)) {
      arg::err("the supplied covariance matrix does not have the right dimensions. It should have dimension {fncol(gr)} by {fncol(gr)} but actually has dimension {fnrow(V)} by {fncol(V)}")
    }

    out[["vcov"]] <- quad_mult(gr, V)
  }

  out
}
