process_by <- function(by, data) {
  if (is_null(by)) {
    return(list(by_id = qG(alloc(1L, nrow(data)))))
  }

  if (is.character(by)) {
    by <- reformulate(by)
  }
  else if (!rlang::is_formula(by)) {
    .err("{.arg by} must be a one-sided formula or character vector")
  }

  by_mf <- update(by, NULL ~ .) |>
    model.frame(data = data)

  attr(by_mf, "terms") <- NULL

  by_grid <- funique(by_mf) |>
    roworderv()

  by_grid$.by_id <- seq_row(by_grid)

  by_mf <- by_mf |>
    ftransform(.merge_id = seq_row(by_mf)) |>
    merge(by_grid,
          by = seq_col(by_mf),
          all.x = TRUE, all.y = FALSE,
          sort = FALSE) |>
    roworderv(".merge_id")

  list(by_grid = by_grid[-ncol(by_grid)],
       by_id = qG(by_mf[[length(by_mf)]]))
}

process_subset <- function(index.sub, data, env = parent.frame(2L)) {

  if (is_null(index.sub)) {
    return(seq_row(data))
  }

  subset <- eval(index.sub, data, env)

  if (!is.atomic(subset)) {
    .err("{.arg subset} must evaluate to an atomic vector")
  }

  if (is.logical(subset)) {
    if (length(subset) != nrow(data)) {
      .err("when {.arg subset} is logical, it must have the same length as the original dataset")
    }

    return(which(subset))
  }

  if (is.numeric(subset)) {
    s <- alloc(FALSE, nrow(data))
    if (all(subset > 0)) {
      s[subset] <- TRUE
    }
    else {
      s <- alloc(TRUE, nrow(data))
      s[subset] <- FALSE
    }
  }
  else if (is.character(subset)) {
    s <- setNames(alloc(FALSE, nrow(data)),
                  rownames(data))
    s[subset] <- TRUE
  }

  which(s)
}

process_subset_by_grid <- function(index.sub, .by_grid = NULL, .contrast = NULL) {

  if (is_null(index.sub)) {
    if (is_null(.by_grid)) {
      return(TRUE)
    }

    return(alloc(TRUE, nrow(.by_grid)))
  }

  if (is_not_null(.contrast)) {
    .err("{.arg subset} cannot be specified with {.fun curve_contast}")
  }

  if (is_null(.by_grid)) {
    .err("{.arg subset} cannot be specified when {.arg by} was not used in the original call to {.fun adrf}")
  }

  vars_in_subset <- get_varnames(index.sub)

  if (!all(vars_in_subset %in% names(.by_grid))) {
    .err("all variables named in {.arg subset} must be present in the original {.arg by} argument supplied to {.fun adrf}")
  }

  subset <- eval(index.sub, .by_grid)

  if (!is.atomic(subset) || !is.logical(subset)) {
    .err("{.arg subset} must evaluate to a logical vector")
  }

  if (!any(subset)) {
    .err("{.arg subset} must identify a real subset of units")
  }

  subset
}

process_range <- function(range = .95, n, treat_var, w = NULL, strict = FALSE) {
  arg_numeric(range)

  if (length(range) > 2L) {
    .err("{.arg range} must have length 1 or 2")
  }

  arg_count(n)
  arg_range(n, c(2L, 1000L))

  if (length(range) > 1L) {
    range <- sort(range)

    range_t <- .range(treat_var)

    if (strict) {
      if (any(range < range_t[1L]) || any(range > range_t[2L])) {
        r <- format(range_t, digits = 4L, drop0trailing = TRUE)
        .err("no values in {.arg range} can be outside the range of observed treatment values ({r[1L]} to {r[2L]})")
      }
    }
    else {
      range_diff <- range_t[2L] - range_t[1L]

      if (any(range < range_t[1L] - .1 * range_diff) ||
          any(range > range_t[2L] + .1 * range_diff)) {
        .err("some values in {.arg range} are outside the range of observed treatment values")
      }
    }
  }
  else if (range == 1) {
    range <- .range(treat_var)
  }
  else if (range > 0 && range < 1) {
    range <- .quantile(treat_var, probs = c(1 - range, 1 + range) / 2,
                       w = w)
  }
  else {
    .err("if {.arg range} is a single number, it must be a quantile (i.e., between 0 and 1)")
  }

  seq(range[1L], range[2L], length.out = n)
}

check_reference <- function(reference, values, strict = TRUE) {
  arg_number(reference)

  range_t <- .range(values)

  range_diff <- range_t[2L] - range_t[1L]

  if (strict) {
    if (any(reference < range_t[1L] - .01 * range_diff) ||
        any(reference > range_t[2L] + .01 * range_diff)) {
      r <- format(range_t, digits = 4L, drop0trailing = TRUE)
      .err("{.arg reference} must be within the range of treatment values supplied to {.fun adrf} ({r[1L]} to {r[2L]})")
    }
  }
  else {
    if (any(reference < range_t[1L] - .1 * range_diff) ||
        any(reference > range_t[2L] + .1 * range_diff)) {
      .wrn("{.arg reference} is outside the range of observed treatment values supplied to {.fun adrf}")
    }
  }
}

process_eps <- function(eps, values) {
  arg_number(eps)
  arg_gt(eps, 0)

  eps * diff1(frange(values)) / 4
}

process_model_data <- function(model, data = NULL) {
  if (is.environment(model[["data"]])) {
    vars <- insight::find_predictors(model, effects = "fixed", component = "all",
                                     flatten = TRUE, verbose = FALSE)

    dat <- model[["data"]] |> as.list() |> list2DF()

    return(ss(dat, j = intersect(vars, colnames(dat))))
  }

  if (utils::hasName(model, "survey.design")) {
    vars <- insight::find_predictors(model, effects = "fixed", component = "all",
                                     flatten = TRUE, verbose = FALSE)

    dat <- insight::get_data(model, verbose = FALSE,
                             source = "mf")

    return(ss(dat, j = intersect(vars, colnames(dat))))
  }

  dat <- try(insight::get_predictors(model, verbose = FALSE),
             silent = TRUE)

  if (!null_or_error(dat)) {
    return(dat)
  }

  if (is_null(data)) {
    if (insight::is_model_supported(model)) {
      .err("{.fun insight::get_predictors} could not extract the dataset used to fit the model. Please supply the dataset to {.arg data}")
    }
    else {
      .err("because the fitted model is not supported by {.pkg insight}, {.arg data} must be supplied")
    }
  }

  arg_data(data)

  data
}

process_model_data_mi <- function(model, data = NULL) {
  data_list <- get_data_mids(model)

  if (any(lengths(data_list) == 0L)) {
    if (is_null(data)) {
      .err("the dataset used to fit the model could not be extracted. Please supply the dataset to {.arg data}")
    }

    arg_is(data, "mids")

    m <- data[["m"]]

    if (!is.numeric(m) || length(m) != 1L || m != length(model[["analyses"]])) {
      .err("the {.cls mids} object supplied to {.arg data} must have the same number of imputed datasets as models fit in {.arg x}")
    }

    rlang::check_installed("mice")

    data_list[] <- mice::complete(data, action = "all")
  }

  rowbind(data_list, idcol = ".imp")
}

get_data_mids <- function(x, imp = NULL) {
  m <- length(x[["analyses"]])

  if (is_null(imp)) {
    imp <- seq_len(m)
  }

  data_list <- make_list(m)

  for (.i in imp) {
    xi <- x[["analyses"]][[.i]]

    vars <- insight::find_predictors(xi, effects = "fixed", component = "all",
                                     flatten = TRUE, verbose = FALSE)

    if (is.data.frame(xi[["data"]])) {
      dat <- xi[["data"]]

      data_list[[.i]] <- ss(dat, j = intersect(vars, names(dat)))
    }
    else {
      dat <- try({
        l <- {
          if (is.environment(xi[["data"]])) xi[["data"]]
          else insight::find_formula(xi)$conditional |>
            environment()
        } |>
          as.list()

        l[intersect(vars, names(l))] |> list2DF()
      }, silent = TRUE)

      if (!null_or_error(dat)) {
        data_list[[.i]] <- dat
      }
    }
  }

  if (length(imp) == 1L) {
    return(data_list[[imp]])
  }

  data_list
}

process_treat <- function(treat, data) {
  arg_string(treat)

  if (!any(names(data) == treat)) {
    .err(c("{.arg treat} must be the name of a variable in the original dataset.",
           "*" = "Supplied value: {.val {treat}}",
           "*" = "Available names: {.or {.val {names(data)}}}"))
  }

  treat_var <- data[[treat]]

  if (!is.numeric(treat_var)) {
    .err("{.arg treat} must be the name of a numeric variable in the original dataset")
  }

  treat_var
}

process_wts <- function(wts, model, data = NULL) {
  if (is_null(wts)) {
    if (endsWith(class(model)[1L], "_weightit") &&
        utils::hasName(model.frame(model), "(s.weights)")) {
      wts <- "(s.weights)"
    }
    else if (inherits(model, "svyglm")) {
      wts <- TRUE
    }
    else {
      wts <- FALSE
    }
  }
  else if (!rlang::is_string(wts) && !rlang::is_bool(wts) && !is.numeric(wts)) {
    .err("{.arg wts} must be {.val {list(NULL)}}, {.val {TRUE}}, {.val {FALSE}}, a numeric vector, or a string containing the name of the variable containing weights in the original dataset")
  }

  if (isTRUE(wts)) {
    wts <- insight::get_weights(model) %or% FALSE
  }

  if (rlang::is_string(wts)) {
    if (is_null(data)) {
      if (wts %in% c("(weights)", "(s.weights)") && is_not_null(model[["model"]]) &&
          utils::hasName(model[["model"]], wts)) {
        data <- model[["model"]]
      }
      else {
        data <- try(insight::get_data(model, verbose = FALSE),
                    silent = FALSE)
      }
    }

    if (!is.data.frame(data)) {
      .err("no dataset could be found to extract the weights named in {.arg wts}. Please supply a dataset to {.arg data}")
    }

    if (!utils::hasName(data, wts)) {
      .err("when supplied as a string, {.arg wts} must be the name of the variable in the original dataset")
    }

    wts <- data[[wts]]

    if (!is.numeric(wts)) {
      .err("when supplied as a string, {.arg wts} must be the name of a numeric variable")
    }

    return(as.vector(wts))
  }

  n <- insight::n_obs(model)

  if (is_null(n)) {
    if (is_null(data)) {
      data <- try(insight::get_data(model, verbose = FALSE),
                  silent = FALSE)

      if (!is.data.frame(data)) {
        .err("no dataset could be found to process the weights supplied to {.arg wts}. Please supply a dataset to {.arg data}")
      }
    }
    else {
      arg_data(data)
    }

    n <- nrow(data)
  }

  if (isFALSE(wts)) {
    return(alloc(1.0, n))
  }

  if (length(wts) != n) {
    .err("when supplied as a numeric vector, {.arg wts} must have as many values as the number of observations in the original model")
  }

  as.vector(wts)
}

process_wts_mi <- function(wts, model, data = NULL) {
  if (is_not_null(wts) && !rlang::is_string(wts) &&
      !rlang::is_bool(wts) && !is.numeric(wts)) {
    .err("{.arg wts} must be {.val {list(NULL)}}, {.val {TRUE}}, {.val {FALSE}}, a numeric vector, or a string containing the name of the variable containing weights in the original dataset")
  }

  m <- length(model[["analyses"]])

  wts_list <- make_list(m)

  if (is_null(wts)) {
    wts <- lapply(model[["analyses"]],
                  function(mod) {
                    if (endsWith(class(mod)[1L], "_weightit")) "(s.weights)"
                    else inherits(mod, "svyglm")
                  }) |>
      unlist()
  }

  if (!is.numeric(wts)) {
    if (length(wts) == 1L) {
      wts <- rep.int(wts, m)
    }

    for (.i in seq_len(m)) {
      if (isTRUE(wts[.i])) {
        test_w <- try(insight::get_weights(model[["analyses"]][[.i]]))

        if (is.numeric(test_w)) {
          wts_list[[.i]] <- as.vector(test_w)
          next
        }

        wts[.i] <- FALSE
      }

      if (isFALSE(wts[.i])) {
        n.i <- insight::n_obs(model[["analyses"]][[.i]])

        if (is_null(n.i)) {
          if (is_null(data)) {
            data.i <- try(get_data_mids(model, .i), silent = TRUE)

            if (!is.data.frame(data.i)) {
              .err("no dataset could be found to extract the weights named in {.arg wts}. Please supply a dataset to {.arg data}")
            }
          }
          else {
            arg_is(data, "mids")

            rlang::check_installed("mice")

            data.i <- mice::complete(data, action = .i)
          }

          n.i <- nrow(data.i)
        }

        wts_list[[.i]] <- alloc(1.0, n.i)
        next
      }

      if (is_null(data)) {
        if (wts[.i] %in% c("(weights)", "(s.weights)") && is_not_null(model[["analyses"]][[.i]][["model"]]) &&
            utils::hasName(model[["analyses"]][[.i]][["model"]], wts[.i])) {
          data.i <- model[["analyses"]][[.i]][["model"]]
        }
        else {
          data.i <- try(get_data_mids(model, .i), silent = TRUE)
        }

        if (!is.data.frame(data.i)) {
          .err("no dataset could be found to extract the weights named in {.arg wts}. Please supply a dataset to {.arg data}")
        }
      }
      else {
        arg_is(data, "mids")

        rlang::check_installed("mice")

        data.i <- mice::complete(data, action = .i)
      }

      if (!utils::hasName(data.i, wts[.i])) {
        .err("when supplied as a string, {.arg wts} must be the name of the variable in the original dataset")
      }

      test_w <- data.i[[wts[.i]]]

      if (!is.numeric(test_w)) {
        .err("when supplied as a string, {.arg wts} must be the name of a numeric variable")
      }

      wts_list[[.i]] <- as.vector(test_w)
    }

    return(unlist(wts_list))
  }

  ns <- integer(m)

  for (.i in seq_len(m)) {
    n.i <- try(insight::n_obs(model[["analyses"]][[.i]]),
               silent = TRUE)

    if (is_count(n.i)) {
      ns[.i] <- n.i
      next
    }

    if (is_not_null(data)) {
      arg_is(data, "mids")

      rlang::check_installed("mice")

      ns[] <- vapply(mice::complete(data, action = "all"),
                     nrow, integer(1L))
      break
    }

    data.i <- try(get_data_mids(model, .i), silent = TRUE)

    if (!is.data.frame(data.i)) {
      .err("no dataset could be found to process the weights supplied to {.arg wts}. Please supply a dataset to {.arg data}")
    }

    ns[.i] <- nrow(data.i)
  }

  if (length(wts) == sum(ns)) {
    return(as.vector(wts))
  }

  if (!all_the_same(ns)) {
    .err("when {.arg wts} is supplied as a numeric vector with multiply imputed data, it must have length equal to the sum of the number of observations across all imputed datasets, or the imputed datasets must have the same number of observations")
  }

  if (length(wts) != ns[1L]) {
    .err("when {.arg wts} is supplied as a numeric vector with multiply imputed data, it must have length equal to the sum of the number of observations across all imputed datasets or to the number of observations in each imputed dataset")
  }

  rep.int(as.vector(wts), m)
}

process_transform <- function(transform, x = NULL, .est = NULL) {
  inv_transform <- NULL
  d_transform <- NULL

  if (inherits(x, "effect_curve")) {
    if (!.is_pure_adrf(x)) {
      transform <- FALSE
    }
  }
  else if (!identical(get_curve_type(x), "ADRF") ||
           is_not_null(.attr(x, ".contrast")) ||
           is_not_null(.attr(x, ".reference"))) {
    transform <- FALSE
  }

  if (isFALSE(transform)) {
    .family <- stats::quasi("identity")
  }
  else if (isTRUE(transform)) {
    .family <- .attr(x, ".family") %or% stats::quasi("identity")
  }
  else if (inherits(transform, "family")) {
    .family <- transform
  }
  else if (inherits(transform, "link-glm")) {
    .family <- stats::quasi(transform)
  }
  else if (is.function(transform)) {
    if (identical(transform, log)) {
      .family <- stats::quasi("log")
    }
    else if (identical(transform, sqrt)) {
      .family <- stats::quasi("sqrt")
    }
    else if (identical(transform, qlogis)) {
      .family <- stats::quasi("logit")
    }
    else if (identical(transform, qnorm)) {
      .family <- stats::quasi("probit")
    }
    else if (identical(transform, exp)) {
      .family <- list(linkfun = exp,
                      linkinv = log,
                      dlinkfun = exp)
    }
    else {
      .family <- list(linkfun = transform)
    }
  }
  else if (is.list(transform) && is_not_null(names(transform)) &&
           utils::hasName(transform, "transform") &&
           all(names(transform) %in% c("transform", "inv_transform", "d_transform"))) {

    .family <- list(linkfun = transform[["transform"]],
                    linkinv = transform[["inv_transform"]],
                    dlinkfun = transform[["d_transform"]])
  }
  else {
    .err("{.arg transform} must be {.val {TRUE}}, {.val {FALSE}}, a family, or a function")
  }

  transform <- .family$linkfun

  inv_transform <- .family$linkinv %or% make_inverse(transform, .est)

  d_transform <- {
    if (is_not_null(.family$mu.eta)) {
      function(mu) {1 / .family$mu.eta(.family$linkfun(mu))}
    }
    else {
      .family$dlinkfun %or% function(mu) {
        eps <- 1e-5
        (transform(mu + eps) - transform(mu - eps)) / (2 * eps)
      }
    }
  }

  d_inv_transform <- .family$mu.eta %or% function(eta) {
    eps <- 1e-5
    (inv_transform(eta + eps) - inv_transform(eta - eps)) / (2 * eps)
  }

  list(transform = transform,
       inv_transform = inv_transform,
       d_transform = d_transform,
       d_inv_transform = d_inv_transform)
}

process_family <- function(model) {
  if (inherits(model, "bart")) {
    if (is_null(model[["sigma"]])) {
      .family <- stats::binomial("probit")
    }
    else {
      .family <- stats::gaussian()
    }
  }
  else {
    .family <- try(insight::get_family(model),
                   silent = TRUE)

    if (!inherits(.family, "family")) {
      .family <- stats::quasi("identity")
    }
  }

  .family[setdiff(names(.family), c("link", "linkfun", "linkinv", "mu.eta"))] <- NULL

  .family
}

process_response <- function(model) {
  insight::find_formula(model, verbose = FALSE)$conditional |>
    rlang::f_lhs() |>
    deparse1()
}

process_vcov_and_cluster <- function(vcov, model, cluster = NULL, is_bayes = FALSE,
                                     data = NULL) {
  if (is_null(vcov) || isFALSE(vcov) || identical(vcov, "none")) {
    if (is_not_null(cluster)) {
      .v <- if (is_null(vcov)) list(NULL) else vcov

      .wrn("{.arg cluster} is ignored when {.code vcov = {.val {(.v)}}}")
    }

    return(list(vcov = "none"))
  }

  if (is_bayes) {
    if (identical(vcov, "unconditional")) {
      vcov <- "posterior"
    }
    else {
      builtin_vcovs <- c("none", "posterior")

      vcov <- match_arg(vcov, builtin_vcovs, context = "with Bayesian models,")

      if (vcov == "none") {
        return(Recall(vcov = "none", model = model, cluster = cluster,
                      is_bayes = is_bayes, data = data))
      }
    }

    if (is_not_null(cluster)) {
      .wrn("{.arg cluster} is ignored with Bayesian models")
    }

    return(list(vcov = vcov))
  }

  builtin_vcovs <- c("none", "unconditional", "conditional", "boot", "fwb")

  if (rlang::is_string(vcov)) {
    p <- pmatch(vcov, builtin_vcovs, nomatch = 0L)

    if (p != 0L) {
      if (builtin_vcovs[p] == "none") {
        return(Recall(vcov = "none", model = model, cluster = cluster,
                      is_bayes = is_bayes, data = data))
      }

      if (builtin_vcovs[p] %in% c("boot", "fwb")) {
        rlang::check_installed("fwb (>= 0.5.0)")

        if (is_not_null(cluster)) {
          cl <- process_cluster(cluster, model, data)

          if (length(cl[["cluster"]]) > 1L) {
            .err("only one level of clustering is allowed with {.code vcov = {.val {builtin_vcovs[p]}}}")
          }

          return(c(list(vcov = builtin_vcovs[p]),
                   cl))
        }

        return(list(vcov = builtin_vcovs[p]))
      }

      if (builtin_vcovs[p] == "unconditional") {
        if (null_or_error(try(sandwich::bread(model), silent = TRUE)) ||
            null_or_error(try(sandwich::estfun(model), silent = TRUE))) {
          .wrn('{.arg vcov} cannot be {.val {"unconditional"}} with this type of model. Using {.val {"conditional"}}')

          return(Recall(vcov = "conditional", model = model, cluster = cluster,
                        is_bayes = is_bayes, data = data))
        }

        return(c(list(vcov = builtin_vcovs[p]),
                 process_cluster(cluster, model, data)))
      }

      # Use default get_vcov() below
      vcov <- NULL
    }
    else if (is_not_null(cluster) && vcov %in% c("HC0", "HC1", "HC2")) {
      cl <- process_cluster(cluster, model, data)

      vcov_try <- try(sandwich::vcovCL(model, type = vcov, cluster = cluster),
                      silent = TRUE)

      if (null_or_error(vcov_try)) {
        .err("the cluster covariance couldn't be calculated")
      }

      return(c(list(vcov = vcov_try),
               cl))
    }
  }

  vcov_try <- try(marginaleffects::get_vcov(model, vcov = cluster %or% vcov),
                  silent = TRUE)

  if (null_or_error(vcov_try)) {
    .err("{.arg vcov} must be either one of {.or {.val {builtin_vcovs}}} or one of the allowed arguments to {.fun marginaleffects::get_vcov}")
  }

  c(list(vcov = vcov_try),
    process_cluster(cluster, model, data))
}

process_vcov_and_cluster_mi <- function(vcov, models, cluster = NULL, is_bayes = FALSE,
                                        data_complete = NULL, imp_split) {
  if (is_null(vcov) || isFALSE(vcov) || identical(vcov, "none")) {
    if (is_not_null(cluster)) {
      .v <- if (is_null(vcov)) list(NULL) else vcov

      .wrn("{.arg cluster} is ignored when {.code vcov = {.val {(.v)}}}")
    }

    return(rep_with(list(vcov = "none"),
                    models))
  }

  if (is_bayes) {
    if (identical(vcov, "unconditional")) {
      vcov <- "posterior"
    }
    else {
      builtin_vcovs <- c("none", "posterior")

      vcov <- match_arg(vcov, builtin_vcovs, context = "with multiply imputed Bayesian models,")

      if (vcov == "none") {
        return(Recall(vcov = "none", models = models, cluster = cluster,
                      is_bayes = is_bayes, data_complete = data_complete,
                      imp_split = imp_split))
      }
    }

    if (is_not_null(cluster)) {
      .wrn("{.arg cluster} is ignored with Bayesian models")
    }

    return(rep_with(list(vcov = vcov),
                    models))
  }

  builtin_vcovs <- c("none", "unconditional", "conditional", "boot", "fwb")

  if (rlang::is_string(vcov)) {
    p <- pmatch(vcov, builtin_vcovs, nomatch = 0L)

    if (p != 0L) {
      if (builtin_vcovs[p] == "none") {
        return(Recall(vcov = "none", models = models, cluster = cluster,
                      is_bayes = is_bayes, data_complete = data_complete,
                      imp_split = imp_split))
      }

      if (builtin_vcovs[p] %in% c("boot", "fwb")) {
        .err("{.arg vcov} cannot be {.val {builtin_vcovs[p]}} with multiply imputed data")
      }

      if (builtin_vcovs[p] == "unconditional") {
        for (model in models) {
          if (null_or_error(try(sandwich::bread(model), silent = TRUE)) ||
              null_or_error(try(sandwich::estfun(model), silent = TRUE))) {
            .wrn('{.arg vcov} cannot be {.val {"unconditional"}} with this type of model. Using {.val {"conditional"}}')

            return(Recall(vcov = "conditional", models = models, cluster = cluster,
                          is_bayes = is_bayes, data_complete = data_complete,
                          imp_split = imp_split))
          }
        }

        out_list <- lapply(seq_along(models), function(.i) {
          c(list(vcov = builtin_vcovs[p]),
            process_cluster(cluster, models[[.i]], ss(data_complete, imp_split[[.i]])))
        })

        return(out_list)
      }

      # Use default get_vcov() below
      vcov <- NULL
    }
    else if (is_not_null(cluster) && vcov %in% c("HC0", "HC1", "HC2")) {
      vcov_try_list <- lapply(models, function(model) {
        try(sandwich::vcovCL(model, type = vcov, cluster = cluster),
            silent = TRUE)
      })

      if (any_apply(vcov_try_list, null_or_error)) {
        bad_imps <- which(vapply(vcov_try_list, null_or_error, logical(1L)))
        .err("the cluster covariance could not be calculated for imputations {.or {bad_imps}}")
      }

      out_list <- lapply(seq_along(models), function(.i) {
        c(list(vcov = vcov_try_list[[.i]]),
          process_cluster(cluster, models[[.i]], ss(data_complete, imp_split[[.i]])))
      })

      return(out_list)
    }
  }

  if (!is.list(vcov)) {
    vcov <- rep_with(list(vcov), models)
  }
  else if (length(vcov) != length(models)) {
    .err("{.arg vcov} must have length equal to the number of imputed datasets")
  }

  for (.i in seq_along(models)) {
    vcov[[.i]] <- try(marginaleffects::get_vcov(models[[.i]], vcov = cluster %or% vcov[[.i]]),
                      silent = TRUE)

    if (null_or_error(vcov[[.i]])) {
      .err("{.arg vcov} must be one of {.or {c(add_quotes(setdiff(builtin_vcovs, c('bs', 'fwb'))), 'one of the allowed arguments to')}} {.fun marginaleffects::get_vcov} or a list thereof")
    }
  }

  lapply(seq_along(models), function(.i) {
    c(list(vcov[[.i]]),
      process_cluster(cluster, models[[.i]], ss(data_complete, imp_split[[.i]])))
  })
}

process_cluster <- function(.cluster, model, data) {
  if (is_null(.cluster)) {
    return(NULL)
  }

  if (rlang::is_formula(.cluster)) {
    cluster <- try({
      cluster_tmp <- expand.model.frame(model, .cluster, na.expand = FALSE)

      model.frame(.cluster, data = cluster_tmp,
                  na.action = na.pass)},
      silent = TRUE)

    if (null_or_error(cluster)) {
      cluster <- try(model.frame(.cluster, data = data, na.action = na.pass),
                     silent = TRUE)
    }

    if (null_or_error(cluster)) {
      .err("clusters could not be extracted from the model; trying supplying {.arg cluster} as a data frame or list")
    }
  }
  else {
    cluster <- as.data.frame(.cluster)
  }

  if (nrow(cluster) != nrow(data)) {
    .err("the number of observations in {.arg cluster} must equal that in the original dataset")
  }

  if (anyNA(cluster, recursive = TRUE)) {
    .err("{.val {NA}}s are not allowed in the {cli::qty(ncol(cluster))} clustering variable{?s}")
  }

  cluster[] <- lapply(cluster, qF)

  p <- length(cluster)
  if (p > 1L) {
    clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE)) |>
      unlist(recursive = FALSE)
    sgn <- (-1L)^(lengths(clu) + 1L)
    cluster <- lapply(clu, function(i) finteraction(cluster[i]))
  }
  else {
    sgn <- 1
  }

  #For small sample adjustment (setting cadjust = TRUE in vcovCL)
  g <- vapply(cluster, fnlevels, numeric(1L))

  list(cluster = cluster,
       adj = sgn * g / (g - 1))
}

get_tt <- function(model, model_data = NULL) {
  tt <- terms(model) |>
    delete.response()

  if (is_null(.attr(tt, "dataClasses"))) {
    tt <- model.frame(tt, data = model_data) |>
      terms()
  }

  tt
}

process_conf_level <- function(conf_level = .95) {
  arg_number(conf_level)
  arg_range(conf_level, c(0, 1), inclusive = c(TRUE, FALSE))

  conf_level
}

process_ci.type <- function(ci.type, .vcov_type, simultaneous = NULL) {
  if (.vcov_type == "none") {
    return(NULL)
  }

  if (.vcov_type == "posterior") {
    arg_string(ci.type)

    ci.type <- match_arg(ci.type, c("perc", "wald"),
                         context = "with Bayesian models,")

    return(ci.type)
  }

  if (.vcov_type != "bootstrap") {
    return("wald")
  }

  arg_string(ci.type)

  if (!isTRUE(simultaneous)) {
    return(ci.type)
  }

  match_arg(ci.type, c("perc", "wald"),
            context = "when {.code simultaneous = TRUE},")
}

process_simultaneous <- function(simultaneous, .est) {
  if (length(.est) <= 1) {
    return(FALSE)
  }

  arg_flag(simultaneous)

  simultaneous
}

process_df <- function(fit) {
  if (inherits(fit, "mira")) {
    .df <- vapply(fit[["analyses"]], process_df, numeric(1L)) |>
      min()

    return(.df)
  }

  if (!insight::is_model_supported(fit) ||
      insight::is_bayesian_model(fit)) {
    return(Inf)
  }

  statistic <- insight::find_statistic(fit)

  if (identical(statistic, "chi-squared statistic")) {
    return(Inf)
  }

  insight::get_df(fit, type = "wald", statistic = statistic)
}

process_null <- function(null = NULL, x = NULL) {
  if (is_not_null(null)) {
    if (length(null) == 1L && allNA(null)) {
      null <- NA_real_
    }
    else {
      arg_number(null)
    }
  }
  else if (inherits(x, "effect_curve")) {
    null <- if (.is_pure_adrf(x)) NA_real_ else 0
  }
  else if (identical(get_curve_type(x), "ADRF") &&
           is_null(.attr(x, ".contrast")) &&
           is_null(.attr(x, ".reference"))) {
    null <- NA_real_
  }
  else {
    null <- 0
  }

  null
}

check_proj <- function(x, proj = NULL) {
  if (is_not_null(proj)) {
    arg_is(proj, "curve_projection")

    if (any_apply(c(".values", ".treat", ".vcov_type", ".curve_type", ".response",
                    ".by_grid", ".contrast", ".reference", ".family", ".df"),
                  function(a) !identical(.attr(x, a), .attr(proj, a)))) {
      .err("{.arg proj} must be the output of {.fun curve_projection} applied to the {.cls effect_curve} object")
    }
  }
}

get_curve_type <- function(x) {
  .curve_type <- .attr(x, ".curve_type")

  if (is_null(.curve_type)) {
    if (inherits(x, "adrf_curve")) {
      .curve_type <- "ADRF"
    }
    else if (inherits(x, "amef_curve")) {
      .curve_type <- "AMEF"
    }
  }

  .curve_type
}

check_effect_curve <- function(x, amef_ok = TRUE, reference_ok = TRUE,
                               contrast_ok = TRUE, projection_ok = TRUE) {
  nm <- deparse1(substitute(x))

  arg_is(x, "effect_curve")

  if (!amef_ok && inherits(x, "amef_curve")) {
    .err("{.arg {nm}} cannot be an {.cls amef_curve} object")
  }

  if (!reference_ok && inherits(x, "reference_curve")) {
    .err("{.arg {nm}} cannot be a {.cls reference_curve} object")
  }

  if (!contrast_ok && inherits(x, "contrast_curve")) {
    .err("{.arg {nm}} cannot be a {.cls contrast_curve} object")
  }

  if (!projection_ok && inherits(x, "curve_projection")) {
    .err("{.arg {nm}} cannot be a {.cls curve_projection} object")
  }
}

process_bw <- function(bw, v, constant = .5) {
  if (is_null(bw)) {
    spacing <- NULL
    diff_v <- diff1(v)

    if (all_the_same(diff_v)) {
      spacing <- mean(diff_v)
    }
    else {
      spacing_table <- table(round(diff_v, 8L))

      spacing <- as.numeric(names(spacing_table)[which.max(spacing_table)])

      if (mean(check_if_zero(diff_v - spacing)) < .7) {
        spacing <- NULL
      }
    }

    if (is_null(spacing)) {
      return(bw.nrd(v))
    }

    arg_number(constant)

    return(constant * spacing)
  }

  if (is.character(bw)) {
    bw <- sprintf("bw.%s", tolower(bw)) |>
      match.fun()
  }

  if (is.function(bw)) {
    bw <- match.fun(bw)

    bw <- bw(v)
  }

  if (!is.numeric(bw)) {
    .err("{.arg bw} must be a string, function, or number")
  }

  arg_number(bw)
  arg_gt(bw, 0)

  bw
}

get_n_by <- function(.contrast = NULL, .by_grid = NULL) {
  if (is_not_null(.contrast)) {
    return(length(.contrast))
  }

  if (is_not_null(.by_grid)) {
    return(nrow(.by_grid))
  }

  1L
}

add_est_labels <- function(res, .contrast, .by_grid, .values = NULL, est_name = NULL) {

  add_values <- is_not_null(.values)

  if (!add_values) {
    .values <- NA
  }

  # Apply by_grid and contrast labels if necessary
  if (is_not_null(.contrast)) {
    contrast <- factor(rep(.contrast, each = length(.values)),
                       levels = .contrast)

    if (add_values) {
      res[[est_name]] <- rep.int(.values, length(.contrast))
    }

    res <- cbind(contrast = contrast, res)
  }
  else if (is_not_null(.by_grid)) {
    res_names <- names(res)
    val_id <- seq_along(.values)
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
      res[[est_name]] <- .values[res[["..vi"]]]
    }

    res <- ss(res, j = c(names(.by_grid), res_names))
  }
  else if (add_values) {
    res[[est_name]] <- .values
  }

  res
}

# Make nearest positive semidefinite matrix; adapted from
# corpcor::make.positive.definite()
.to_psd <- function(m) {
  d <- nrow(m)

  es <- eigen(m, symmetric = TRUE)
  esv <- es$values

  tol <- d * max(abs(esv)) * .Machine$double.eps

  if (all(esv > tol)) {
    return(m)
  }

  delta <- 2 * tol

  tau <- pmax(0, delta - esv)

  m + tcrossprod(es$vectors %*% diag(sqrt(tau), d))
}

.is_pure_adrf <- function(x) {
  inherits(x, "effect_curve") &&
    inherits(x, "adrf_curve") &&
    !inherits(x, "amef_curve") &&
    !inherits(x, "contrast_curve") &&
    !inherits(x, "reference_curve")
}

.is_bayesian <- function(x, data = NULL) {
  if (insight::is_model_supported(x)) {
    return(insight::is_bayesian_model(x))
  }

  if (inherits(x, "bart")) {
    return(TRUE)
  }

  if (is_null(data)) {
    data <- insight::get_data(x, verbose = FALSE)
  }

  p <- try(marginaleffects::get_predict(x, newdata = ss(data, seq_len(min(nrow(data), 2L)))),
           silent = TRUE)

  is_not_null(.attr(p, "posterior_draws"))
}

# Make a function that yields its input when fun is applied to its output
make_inverse <- function(fun, original) {
  if (all(original == fun(original))) {
    return(identity)
  }

  function(x) {
    # length 1 x doesn't always play nice, so repeat it
    d <- length(x) == 1L
    if (d) {
      x <- x[c(1L, 1L)]
    }

    # Get starter closest to transformed x
    .start <- original[vapply(x, function(y) which.min(abs(fun(original) - y)),
                              integer(1L))]

    out <- optim(function(y) mean((fun(y) - x)^2),
                 par = .start,
                 control = list(maxit = 1e3,
                                reltol = 1e-12))$par
    if (d) {
      out <- mean(out)
    }

    out
  }
}

make_rmvt <- function(mu, Sigma, df = Inf, tol = 1e-7) {
  #Returns p x n matrix (easier to get col maxes)
  p <- length(mu)

  if (!allv(dim(Sigma), p)) {
    .err("incompatible arguments")
  }

  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values

  if (any(ev < -tol * abs(ev[1L]))) {
    .err("{.arg Sigma} must be positive definite")
  }

  mu <- drop(mu)
  scale_mat <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), p)

  if (is.finite(df)) {
    #Shifted t-distribution
    function(n) {
      X <- matrix(rnorm(n * p), nrow = n, ncol = p, byrow = TRUE)

      mu + (tcrossprod(scale_mat, X) %r/% sqrt(rchisq(n, df) / df))
    }
  }
  else {
    function(n) {
      X <- matrix(rnorm(n * p), nrow = n, ncol = p, byrow = TRUE)

      mu + tcrossprod(scale_mat, X)
    }
  }
}

#' Compute p-value of quadratic form by simulation
#' @noRd
#' @param stat value of `sum(x^2)` to test
#' @param v covariance matrix of `x` under H0
#' @param n number of iterations
#' @param max_size longest length allowed for a given vector
#' @param df degrees of freedom
.sim_p_value <- function(stat, v, n = 1e6, max_size = 1e7, df = Inf) {
  # Make multivariate t/normal random generator
  rmvt <- make_rmvt(alloc(0.0, ncol(v)), v, df = df)

  d <- nrow(v)

  # Max number of iterations in one run, based on max_size
  max_n <- max(1, max_size %/% d)

  if (n <= max_n) {
    return(mean(.colSums(rmvt(n)^2, d, n) >= stat))
  }

  iters <- n %/% max_n
  r <- n %% max_n

  p <- numeric(iters + 1L)
  for (i in seq_len(iters)) {
    p[i] <- sum(.colSums(rmvt(max_n)^2, d, max_n) >= stat)
  }

  if (r > 0) {
    p[iters + 1L] <- sum(.colSums(rmvt(r)^2, d, r) >= stat)
  }

  sum(p) / n
}

#' Get kernel weights for interpolation
#' @noRd
#' @param x the predictor values to impute
#' @param v the values on the existing grid
#' @param kernel the name of kernel; `"gaussian"`, `"epanechnikov"`, or `"triangular"`
#' @param constant when `bw` is `NULL` and `v` is on a regularly spaced grid, `constant` times the spacing is used as the bandwidth.
#' @param bw the bandwidth (function or value) to use. If `NULL` and the grid of `v` is not regularly space, `bw.nrd()` is used.
#' @param ... ignored.
get_kernel_w <- function(x, v = x, kernel = "gaussian", constant = .5, bw = NULL, ...) {

  n <- length(v)
  if (length(x) == n && all(check_if_zero(x - v))) {
    return(diag(1, n, n))
  }

  if (all(x %in% v)) {
    w <- outer(x, v, check_if_zero)

    return(w / rowSums(w))
  }

  bw <- process_bw(bw, v, constant)

  kernel <- match_arg(kernel, c("gaussian", "epanechnikov",
                                "triangular"))

  k <- switch(kernel,
              gaussian = function(a, b) {
                dnorm(a - b, sd = bw)
              },
              epanechnikov = function(a, b) {
                bw_a <- bw * sqrt(5)
                ax <- abs(a - b)
                o <- alloc(0, length(ax))
                in_window <- ax < bw_a
                o[in_window] <- .75 * (1 - (ax[in_window] / bw_a)^2) / bw_a
                o
              },
              triangular = function(a, b) {
                bw_a <- bw * sqrt(6)
                ax <- abs(a - b)
                o <- alloc(0, length(ax))
                in_window <- ax < bw_a
                o[in_window] <- (1 - (ax[in_window] / bw_a)) / bw_a
                o
              })

  w <- outer(x, v, k)

  rs <- rowSums(w)

  zeros <- check_if_zero(rs)

  if (any(zeros)) {
    w[zeros, ] <- 0
  }

  if (!all(zeros)) {
    w[!zeros, ] <- w[!zeros, ] / rs[!zeros]
  }

  w
}

#' Get weights for local polynomial interpolation
#' @noRd
#' @param x the predictor values to impute
#' @param v the values on the existing grid
#' @param w a matrix of length(x) by length(v) weights; if null, kernel weights are estimated
#' @param p the degree of the local polynomial. Set to 0 to just NW interpolation.
#' @param ... args passed to `get_kernel_w()` when `w` is `NULL`
get_locpoly_w <- function(x, v = x, w = NULL, p = 3, ...) {
  arg_count(p)

  if (is_null(w)) {
    w <- get_kernel_w(x = x, v = v, ...)
  }

  arg_equal(length(x), nrow(w))
  arg_equal(length(v), ncol(w))

  if (p == 0) {
    return(w)
  }

  n <- length(v)
  ind <- seq_len(n)

  v0 <- c(v, x)

  dv0 <- cbind(1, poly(v0, degree = p, raw = FALSE, simple = TRUE))

  dx <- ss(dv0, -ind)

  lw <- matrix(0, nrow = length(x), ncol = n)

  for (i in seq_along(x)) {
    wi <- w[i, ]
    pos <- which(wi > 1e-10)

    if (is_null(pos)) {
      .err("all kernel weights estimated as 0, indicating severe extrapolation")
    }

    wi <- wi[pos]
    wi <- wi / sum(wi)

    lw[i, pos] <- drop(ss(dx, i) %*% solve(.to_psd(crossprod(ss(dv0, pos) * sqrt(wi))), t(ss(dv0, pos) * wi)))
  }

  lw
}

# Pool multiply imputed estimates using Rubin's rules
pool_est <- function(est.list, vcov.list = NULL) {
  est.mat <- do.call("rbind", est.list)

  out <- list(est = colMeans(est.mat))

  if (is_not_null(vcov.list)) {
    m <- length(vcov.list)

    v_w <- Reduce("+", vcov.list) / m

    v_b <- cov(est.mat)

    out$vcov <- v_w + (1 + 1 / m) * v_b
  }

  out
}

# Get quadrature weights for trapezoidal Riemann sum along given grid
get_trapezoidal_w <- function(grid) {
  val_diff <- diff1(grid)

  w <- (c(0, val_diff) + c(val_diff, 0)) / 2

  w / sum(w)
}

get_by_grid_labels <- function(.by_grid, sep = ", ") {
  do.call(function(...) paste(..., sep = sep),
          lapply(names(.by_grid), function(i) {
            sprintf("%s = %s",
                    i,
                    add_quotes(.by_grid[[i]], is.character(.by_grid[[i]]) || is.factor(.by_grid[[i]])))
          }))
}

.print_estimate_table <- function(x, digits, topn, main = NULL, help.fn = NULL, bar = TRUE, rownames = FALSE, ...) {
  nr <- nrow(x)

  if (nr > 2 * topn + 1L) {
    head_ind <- seq_len(topn)
    tail_ind <- nr - rev(head_ind) + 1L
  }
  else {
    head_ind <- seq_len(nr)
    tail_ind <- integer()
  }

  if ("p.value" %in% colnames(x)) {
    x[["p.value"]] <- .format_p_value(x[["p.value"]],
                                      digits = digits)
  }

  if (is_not_null(head_ind)) {
    colnames(x) <- vapply(colnames(x),
                          function(nm)
                            switch(nm,
                                   term = "Term",
                                   comparison = "Comparison",
                                   contrast = "Contrast",
                                   estimate = "Estimate",
                                   std.error = "Std. Error",
                                   conf.low = "CI Low",
                                   conf.high = "CI High",
                                   p.value = "P-value",
                                   Chisq = "\u03C7\u00B2",
                                   nm),
                          character(1L))

    x <- ss(x, c(head_ind, tail_ind))

    for (i in which(vapply(x, is.numeric, logical(1L)))) {
      x[[i]] <- zapsmall(x[[i]], digits = digits)
    }

    tmp <- utils::capture.output({
      x |>
        round_df_char(digits = digits, na_vals = ".") |>
        print.data.frame(row.names = rownames,
                         ...)
    })

    out <- tmp[seq_along(c(1L, head_ind))]
  }
  else {
    tmp <- utils::capture.output({
      ss(x, 1L) |>
        round_df_char(digits = digits, na_vals = ".") |>
        print.data.frame(row.names = rownames,
                         ...)
    })

    out <- character(0L)
  }

  to_it <- NULL

  if (nr > 2L * topn + 1) {
    msg <- {
      if (is_null(help.fn))
        sprintf("--- %s rows omitted. ---",
                nr - 2L * topn)
      else
        sprintf("--- %s rows omitted. See ?%s ---",
                nr - 2L * topn, help.fn)
    }

    to_it <- length(out) + 1L

    out <- c(out, center_just(msg, wrt = tmp))

    if (is_not_null(tail_ind)) {
      out <- c(out, tmp[-seq_along(c(1L, head_ind))])
    }
  }

  nc <- max(nchar(out))

  if (is_not_null(to_it)) {
    out[to_it] <- .it(out[to_it])
  }

  if (is_not_null(main)) {
    #Add title
    out <- c(center_just(main, wrt = space(nc)),
             txtbar(nc),
             out)
  }

  #Add bar
  if (bar) {
    out <- c(out, txtbar(nc))
  }

  cat(out, sep = "\n")

  invisible(out)
}

.format_p_value <- function(p, digits = 4L) {
  p_text <- rep.int(NA_character_, length(p))

  min_p <- .1^digits

  nas <- is.na(p)

  p_above <- !nas & p >= min_p
  p_below <- !nas & p < min_p

  if (any(p_above)) {
    p_text[p_above] <- formatC(p[p_above], digits = digits, format = "f")
  }

  if (any(p_below)) {
    p_text[p_below] <- paste("<", formatC(min_p,
                                          digits = digits,
                                          format = "f"))
  }

  p_text
}
