process_by <- function(by, data) {
  if (is_null(by)) {
    return(list(by_grid = NULL,
                by_id = qG(alloc(1L, fnrow(data)))))
  }

  if (is.character(by)) {
    by <- reformulate(by)
  }
  else if (!rlang::is_formula(by)) {
    arg::err("{.arg by} must be a one-sided formula or character vector")
  }

  by_mf <- try(update(by, NULL ~ .) |>
                 model.frame(data = data, na.action = stats::na.pass),
               silent = TRUE)

  if (inherits(by_mf, "try-error")) {
    arg::err(c("{.arg by} could not be evaluated.",
               x = conditionMessage(attr(by_mf, "condition")),
               i = "Variables named in {.arg by} must be in the dataset used to fit the model. If that dataset cannot be found, supply it to {.arg data}."))
  }

  attr(by_mf, "terms") <- NULL

  if (anyNA(by_mf)) {
    arg::err("the variables in {.arg by} cannot contain missing values among the units used to compute the effect curve")
  }

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

  list(by_grid = by_grid[-fncol(by_grid)],
       by_id = qG(by_mf[[length(by_mf)]]))
}

process_subset <- function(index.sub, data, env = parent.frame(2L)) {

  if (is_null(index.sub)) {
    return(seq_row(data))
  }

  subset <- try(eval(index.sub, data, env), silent = TRUE)

  if (inherits(subset, "try-error")) {
    arg::err(c("{.arg subset} could not be evaluated.",
               x = conditionMessage(attr(subset, "condition")),
               i = "Variables named in {.arg subset} must be in the dataset used to fit the model. If that dataset cannot be found, supply it to {.arg data}."))
  }

  if (!is.atomic(subset)) {
    arg::err("{.arg subset} must evaluate to an atomic vector")
  }

  if (is.logical(subset)) {
    if (length(subset) != fnrow(data)) {
      arg::err("when {.arg subset} is logical, it must have one value for each unit used to fit the model")
    }

    return(.check_subset_nonempty(which(subset)))
  }

  if (is.factor(subset)) {
    arg::err(c("{.arg subset} cannot be a factor.",
               i = "Supply a logical vector (e.g., {.code subset = group == \"a\"}), row numbers, or row names."))
  }

  if (is.numeric(subset)) {
    subset <- subset[!is.na(subset)]

    if (!all(check_if_zero(subset - round(subset)))) {
      arg::err("when {.arg subset} is numeric, it must contain whole numbers corresponding to row numbers")
    }

    if (any(subset > 0) && any(subset < 0)) {
      arg::err("{.arg subset} cannot contain both positive and negative row numbers")
    }

    if (any(abs(subset) > fnrow(data))) {
      arg::err("{.arg subset} contains row numbers larger than the number of units used to fit the model ({fnrow(data)})")
    }

    # Indexing the row numbers handles positive and negative indices alike
    s <- seq_row(data) %in% seq_row(data)[subset]
  }
  else if (is.character(subset)) {
    subset <- subset[!is.na(subset)]

    bad <- setdiff(subset, rownames(data))

    if (!is_null(bad)) {
      arg::err(c("When {.arg subset} is a character vector, it must contain row names of units used to fit the model.",
                 x = "{cli::qty(length(bad))}The value{?s} {.val {bad}} {?is/are} not {?a /}row name{?s}."))
    }

    s <- rownames(data) %in% subset
  }
  else {
    arg::err("{.arg subset} must be a logical vector, row numbers, or row names")
  }

  .check_subset_nonempty(which(s))
}

.check_subset_nonempty <- function(s) {
  if (is_null(s)) {
    arg::err("{.arg subset} must select at least one unit")
  }

  s
}

process_subset_by_grid <- function(index.sub, .by_grid = NULL, .contrast = NULL) {

  if (is_null(index.sub)) {
    if (is_null(.by_grid)) {
      return(TRUE)
    }

    return(alloc(TRUE, fnrow(.by_grid)))
  }

  if (!is_null(.contrast)) {
    arg::err("{.arg subset} cannot be specified with {.fun curve_contrast}")
  }

  if (is_null(.by_grid)) {
    arg::err("{.arg subset} cannot be specified when {.arg by} was not used in the original call to {.fun adrf} or {.fun shift_curve}")
  }

  vars_in_subset <- get_varnames(index.sub)

  if (!all(vars_in_subset %in% names(.by_grid))) {
    arg::err("all variables named in {.arg subset} must be present in the original {.arg by} argument supplied to {.fun adrf} or {.fun shift_curve}")
  }

  subset <- eval(index.sub, .by_grid)

  if (!is.atomic(subset) || !is.logical(subset)) {
    arg::err("{.arg subset} must evaluate to a logical vector")
  }

  if (length(subset) == 1L) {
    subset <- rep.int(subset, fnrow(.by_grid))
  }
  else if (length(subset) != fnrow(.by_grid)) {
    arg::err("{.arg subset} must evaluate to a single value or one value for each subgroup")
  }

  subset[is.na(subset)] <- FALSE

  if (!any(subset)) {
    arg::err("{.arg subset} must identify a real subset of units")
  }

  subset
}

process_range <- function(range = .95, n, treat_var, w = NULL, strict = FALSE) {
  arg::arg_numeric(range)
  arg::arg_length(range, 1:2)


  arg::arg_count(n)
  arg::arg_between(n, c(2L, 1000L))

  if (length(range) > 1L) {
    range <- sort(range)

    if (range[1L] == range[2L]) {
      arg::err("the two values in {.arg range} must differ")
    }

    range_t <- .range(treat_var)

    if (strict) {
      if (any(range < range_t[1L]) || any(range > range_t[2L])) {
        r <- format(range_t, digits = 4L, drop0trailing = TRUE)
        arg::err("no values in {.arg range} can be outside the range of observed treatment values ({r[1L]} to {r[2L]})")
      }
    }
    else {
      range_diff <- range_t[2L] - range_t[1L]

      if (any(range < range_t[1L] - .1 * range_diff) ||
          any(range > range_t[2L] + .1 * range_diff)) {
        arg::err("some values in {.arg range} are outside the range of observed treatment values")
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
    arg::err("if {.arg range} is a single number, it must be a quantile (i.e., between 0 and 1)")
  }

  if (!(range[2L] > range[1L])) {
    arg::err("the treatment values in the range requested by {.arg range} are all equal, so an effect curve cannot be computed")
  }

  seq(range[1L], range[2L], length.out = n)
}

# Arguments of shift_curve() supplied to adrf(), or of adrf() supplied to
# shift_curve(), which would otherwise be passed through `...` and give
# confusing errors or be ignored
.check_misplaced_args <- function(dots_names, fun) {
  if (fun == "adrf") {
    bad <- intersect(c("shift", "shift_fun"), dots_names)

    if (!is_null(bad)) {
      arg::err("{.arg {bad}} can only be used with {.fun shift_curve}, which estimates the mean outcome when each unit's treatment is shifted from its observed value")
    }
  }
  else if (anyv(dots_names, "range")) {
    arg::err("{.arg range} cannot be used with {.fun shift_curve}; supply the range of shifts to {.arg shift}")
  }
}

# Processes the `shift_fun` argument of shift_curve(): a string naming a
# built-in shift or a function of the observed treatment values and a single
# shift value. `label` is how the shift is described in printed output.
process_shift_fun <- function(shift_fun, label = NULL) {
  if (is.function(shift_fun)) {
    # args() gives the formals of primitive functions too (or NULL, in which
    # case the function is checked only when it is called)
    fmls <- names(formals(args(shift_fun)))

    if (!is_null(args(shift_fun)) && length(fmls) < 2L && !anyv(fmls, "...")) {
      arg::err("when {.arg shift_fun} is a function, it must take two arguments: the observed treatment values and a single shift value")
    }

    return(list(fun = shift_fun, type = "custom",
                label = label %or% "custom"))
  }

  if (!rlang::is_string(shift_fun)) {
    arg::err("{.arg shift_fun} must be {.val additive}, {.val multiplicative}, or a function")
  }

  type <- arg::match_arg(shift_fun, c("additive", "multiplicative"))

  fun <- switch(type,
                additive = function(a, delta) a + delta,
                multiplicative = function(a, delta) a * (1 + delta))

  list(fun = fun, type = type, label = type)
}

# Processes the `shift` argument of shift_curve() and checks the shift
# function on the observed treatment values, `A.list`, which hold the treatment
# values of the units in `subset` in each imputed dataset (or in the single
# dataset), with their target weights in `w.list`. Returns the information
# stored in the `.shift` attribute of a shift curve, plus the grid of shifts in
# `values`.
process_shift <- function(shift, shift_fun, n, A.list, w.list, treat) {
  arg::arg_count(n)
  arg::arg_between(n, c(2L, 1000L))

  type <- shift_fun[["type"]]
  fun <- shift_fun[["fun"]]

  A <- unlist(A.list, use.names = FALSE)

  if (type == "multiplicative" && fmin(A) < 0) {
    arg::err("{.code shift_fun = \"multiplicative\"} requires the treatment to be nonnegative")
  }

  sd_A <- .w_sd(A, w = unlist(w.list, use.names = FALSE))

  if (is_null(shift)) {
    if (type != "additive") {
      arg::err("{.arg shift} must be supplied when {.arg shift_fun} is not {.val additive}")
    }

    shift <- .5
  }

  arg::arg_numeric(shift)
  arg::arg_length(shift, 1:2)

  if (!all(is.finite(shift))) {
    arg::err("{.arg shift} must contain finite values")
  }

  sd_units <- NULL

  if (length(shift) == 1L) {
    arg::arg_gt(shift, 0)

    bounds <- c(-shift, shift)

    if (type == "additive") {
      if (!is.finite(sd_A) || !(sd_A > 0)) {
        arg::err("the treatment does not vary, so {.arg shift} cannot be in standard deviations of the treatment; supply the bounds of the shift as two numbers")
      }

      sd_units <- bounds
      bounds <- bounds * sd_A
    }
  }
  else {
    bounds <- sort(shift)

    if (bounds[1L] == bounds[2L]) {
      arg::err("the two values in {.arg shift} must differ")
    }
  }

  if (type == "multiplicative" && bounds[1L] < -1) {
    arg::err("with {.code shift_fun = \"multiplicative\"}, {.arg shift} cannot be below -1, which would make the treatment negative")
  }

  values <- seq(bounds[1L], bounds[2L], length.out = n)

  # A shift of 0 computed from the bounds can differ from 0 by rounding error
  values[abs(values) < sqrt(.Machine$double.eps) * diff(bounds)] <- 0

  # The shift function must give a finite treatment value for every unit at
  # every shift
  for (A_i in A.list) {
    for (v in values) {
      a_new <- try(fun(A_i, v), silent = TRUE)

      if (inherits(a_new, "try-error")) {
        arg::err(c("{.arg shift_fun} produced an error:",
                   "!" = conditionMessage(attr(a_new, "condition"))))
      }

      if (!is.numeric(a_new) || length(a_new) != length(A_i) || !is_null(dim(a_new))) {
        arg::err("{.arg shift_fun} must return a numeric vector with one value for each value of the treatment")
      }

      if (!all(is.finite(a_new))) {
        arg::err("{.arg shift_fun} produced non-finite treatment values for some units at some values of {.arg shift}")
      }
    }
  }

  # Whether a shift of 0 leaves every unit's treatment unchanged, which is
  # required for ade()
  identity0 <- all_apply(A.list, function(A_i) {
    a0 <- try(fun(A_i, 0), silent = TRUE)

    is.numeric(a0) && length(a0) == length(A_i) &&
      all(abs(a0 - A_i) <= sqrt(.Machine$double.eps) * pmax(1, abs(A_i)))
  })

  list(values = values,
       fun = fun,
       type = type,
       label = shift_fun[["label"]],
       treat = treat,
       sd = sd_A,
       sd_units = sd_units,
       identity0 = identity0)
}

# The local polynomial used to interpolate between grid points cannot
# extrapolate reliably beyond the grid (just past its ends, the estimates are
# `NA`), so values outside the grid give an error. Values outside it by at most
# 0.1% of its width (e.g., an end point typed as printed, after rounding) are
# moved to the nearest end.
# The function that created the effect curve `x`, for messages
.curve_fun <- function(x) {
  if (is_null(.attr(x, ".shift"))) "adrf" else "shift_curve"
}

.snap_to_grid <- function(x, v, arg_name, shift = FALSE) {
  r <- .range(v)
  leeway <- .001 * (r[2L] - r[1L])

  if (any(x < r[1L] - leeway | x > r[2L] + leeway, na.rm = TRUE)) {
    rf <- format(r, digits = 10L, drop0trailing = TRUE)

    if (shift) {
      arg::err(c("All values of {.arg {arg_name}} must be within the range of shifts at which the shift curve was computed ({rf[1L]} to {rf[2L]}).",
                 i = "To estimate the shift curve outside this range, supply a wider range to {.arg shift} in the call to {.fun shift_curve}."))
    }

    arg::err(c("All values of {.arg {arg_name}} must be within the range of treatment values at which the effect curve was computed ({rf[1L]} to {rf[2L]}).",
               i = "To estimate the effect curve outside this range, supply a wider range to {.arg range} in the call to {.fun adrf}."))
  }

  squish(x, lo = r[1L], hi = r[2L])
}

check_reference <- function(reference, values, strict = TRUE, shift = FALSE) {
  arg::arg_number(reference)

  range_t <- .range(values)

  range_diff <- range_t[2L] - range_t[1L]

  if (strict) {
    return(.snap_to_grid(reference, values, "reference", shift = shift))
  }
  else {
    if (any(reference < range_t[1L] - .1 * range_diff) ||
        any(reference > range_t[2L] + .1 * range_diff)) {
      arg::wrn("{.arg reference} is outside the range of observed treatment values supplied to {.fun adrf}")
    }
  }

  reference
}

process_eps <- function(eps, values) {
  arg::arg_number(eps)
  arg::arg_gt(eps, 0)

  eps * diff1(frange(values)) / 4
}

process_model_data <- function(model, data = NULL) {
  if (utils::hasName(model, "survey.design")) {
    vars <- insight::find_predictors(model, effects = "fixed", component = "all",
                                     flatten = TRUE, verbose = FALSE)

    dat <- insight::get_data(model, verbose = FALSE,
                             source = "mf")

    # For designs that cannot drop rows (e.g., calibrated and PPS designs),
    # this includes the units the fit dropped for missing values
    na_act <- model[["na.action"]]

    if (!is_null(na_act) &&
        fnrow(dat) == nrow(stats::model.matrix(model)) + length(na_act)) {
      dat <- ss(dat, -na_act)
    }

    return(ss(dat, j = intersect(vars, colnames(dat))))
  }

  dat <- try(insight::get_predictors(model, verbose = FALSE),
             silent = TRUE)

  recovered <- !null_or_error(dat)

  if (recovered && .data_matches_model(model, dat)) {
    return(dat)
  }

  if (is_null(data)) {
    if (recovered) {
      arg::err(c("The dataset used to fit the model could not be recovered.",
                 i = "Only the model frame was found, and it contains transformed variables (e.g., from {.fun poly}) rather than the original ones. This can happen when the model is fit inside a function.",
                 i = "Supply the original dataset to {.arg data}."))
    }

    if (insight::is_model_supported(model)) {
      arg::err("{.fun insight::get_predictors} could not extract the dataset used to fit the model. Please supply the dataset to {.arg data}")
    }
    else {
      arg::err("because the fitted model is not supported by {.pkg insight}, {.arg data} must be supplied")
    }
  }

  arg::arg_data.frame(data)

  # Restrict the supplied dataset to the rows used in the fit (e.g., when the
  # model dropped rows with missing values)
  n_fit <- try(nrow(stats::model.frame(model)), silent = TRUE)

  if (!null_or_error(n_fit) && is.numeric(n_fit) && fnrow(data) != n_fit) {
    aligned <- .align_to_model_rows(model, data, n_fit)

    if (is_null(aligned)) {
      arg::err("the dataset supplied to {.arg data} has {fnrow(data)} rows, but the model was fit to {n_fit}, and the rows used could not be matched by row name")
    }

    data <- aligned
  }

  data
}

# The original dataset, with all of its variables, for the rows used in the fit.
# It is used to evaluate `subset`, `by`, `wts`, and `cluster`, which may refer
# to variables that are not in the model. Candidates are the survey design's
# variables, the user-supplied `data`, and the `data` argument of the model
# call; each must reproduce the model frame. Falls back to `model_data` (the
# predictors) when none does.
process_eval_data <- function(model, model_data, data = NULL) {
  candidates <- list(
    if (utils::hasName(model, "survey.design")) model[["survey.design"]][["variables"]],
    data,
    .recover_call_data(model)
  )

  for (cand in candidates) {
    if (!is.data.frame(cand)) {
      next
    }

    cand <- .align_to_model_rows(model, cand, fnrow(model_data))

    if (!is_null(cand) && .data_matches_model(model, cand)) {
      return(cand)
    }
  }

  model_data
}

# Evaluates the `data` argument of the model's call in the environment of its
# formula
.recover_call_data <- function(model) {
  model_call <- try(stats::getCall(model), silent = TRUE)

  if (null_or_error(model_call) || is_null(model_call[["data"]])) {
    return(NULL)
  }

  env <- try(environment(stats::formula(model)), silent = TRUE)

  if (!is.environment(env)) {
    env <- globalenv()
  }

  dat <- try(eval(model_call[["data"]], envir = env), silent = TRUE)

  if (!is.data.frame(dat)) {
    return(NULL)
  }

  dat
}

# Restricts `data` to the rows used in the model fit, matched by row name, when
# it has more rows than were used (e.g., because of missing values or `subset`)
.align_to_model_rows <- function(model, data, n) {
  if (fnrow(data) == n) {
    return(data)
  }

  mf <- try(stats::model.frame(model), silent = TRUE)

  # nrow() rather than fnrow(), which counts the columns of a matrix response
  # (e.g., cbind(successes, failures)) as rows
  if (null_or_error(mf) || nrow(mf) != n) {
    return(NULL)
  }

  idx <- match(rownames(mf), rownames(data))

  if (anyNA(idx)) {
    return(NULL)
  }

  # The original row names are kept so that `subset` can refer to them
  out <- ss(data, idx)
  attr(out, "row.names") <- rownames(data)[idx]

  out
}

# Adds variables used in the model's offsets that are missing from `model_data`
# (e.g., from an `offset` argument rather than an offset() term), taking them
# from `eval_data`
add_offset_vars <- function(model, model_data, eval_data) {
  off_vars <- try(insight::find_offset(model), silent = TRUE)

  if (null_or_error(off_vars) || !is.character(off_vars)) {
    return(model_data)
  }

  off_vars <- setdiff(intersect(off_vars, names(eval_data)), names(model_data))

  if (is_null(off_vars)) {
    return(model_data)
  }

  cbind(model_data, ss(eval_data, j = off_vars))
}

# Checks whether evaluating the model's terms on `data` reproduces the model
# frame. When insight cannot find the original dataset, it falls back to the
# model frame, where a term like poly(x, 2) is stored as its basis; using that
# basis as `x` gives silently wrong counterfactual predictions. Returns TRUE
# when the check cannot be performed.
.data_matches_model <- function(model, data) {
  mf <- try(stats::model.frame(model), silent = TRUE)

  if (null_or_error(mf) || !is.data.frame(mf)) {
    return(TRUE)
  }

  # The candidates are restricted to the rows used in the fit, so a different
  # number of rows means a different dataset (e.g., an object with the same
  # name as the model's `data` found in the environment of its formula)
  if (nrow(mf) != fnrow(data)) {
    return(FALSE)
  }

  tt <- try(stats::terms(mf) |> stats::delete.response(), silent = TRUE)

  if (null_or_error(tt)) {
    return(TRUE)
  }

  rebuilt <- try(stats::model.frame(tt, data = data, na.action = stats::na.pass),
                 silent = TRUE)

  if (null_or_error(rebuilt)) {
    return(FALSE)
  }

  for (nm in intersect(names(rebuilt), names(mf))) {
    a <- mf[[nm]]
    b <- rebuilt[[nm]]

    same <- {
      if (is.numeric(a) && is.numeric(b))
        isTRUE(all.equal(unname(as.matrix(a)), unname(as.matrix(b)),
                         check.attributes = FALSE))
      else identical(as.character(unlist(a)), as.character(unlist(b)))
    }

    if (!same) {
      return(FALSE)
    }
  }

  TRUE
}

process_model_data_mi <- function(model, data = NULL) {
  data_list <- get_data_mids(model)

  if (anyv(vlengths(data_list), 0L)) {
    if (is_null(data)) {
      arg::err("the dataset used to fit the model could not be extracted. Please supply the dataset to {.arg data}")
    }

    arg::arg_is(data, "mids")

    m <- data[["m"]]

    if (!is_number(m) || m != length(model[["analyses"]])) {
      arg::err("the {.cls mids} object supplied to {.arg data} must have the same number of imputed datasets as models fit in {.arg x}")
    }

    rlang::check_installed("mice")

    data_list[] <- mice::complete(data, action = "all")
  }

  rowbind(data_list, idcol = ".imp")
}

# All variables in each imputed dataset, stacked, for evaluating `subset`, `by`,
# and `cluster`; falls back to the stacked predictors when they cannot be
# recovered with the same rows and the same values of the predictors
process_eval_data_mi <- function(model, model_data, data = NULL) {
  data_list <- {
    if (inherits(data, "mids")) mice::complete(data, action = "all")
    else get_data_mids(model, all_vars = TRUE)
  }

  if (!is.list(data_list) || anyv(vlengths(data_list), 0L)) {
    return(model_data)
  }

  dat <- try(rowbind(data_list, idcol = ".imp"), silent = TRUE)

  if (null_or_error(dat) || fnrow(dat) != fnrow(model_data) ||
      !all(names(model_data) %in% names(dat))) {
    return(model_data)
  }

  for (v in setdiff(names(model_data), ".imp")) {
    if (!isTRUE(all.equal(dat[[v]], model_data[[v]], check.attributes = FALSE))) {
      return(model_data)
    }
  }

  dat
}

get_data_mids <- function(x, imp = NULL, all_vars = FALSE) {
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

      data_list[[.i]] <- {
        if (all_vars) dat
        else ss(dat, j = intersect(vars, names(dat)))
      }
    }
    else {
      dat <- try({
        l <- {
          if (is.environment(xi[["data"]])) xi[["data"]]
          else insight::find_formula(xi)$conditional |>
            environment()
        } |>
          as.list()

        l_vars <- l[intersect(vars, names(l))]

        if (all_vars) {
          # Keep every variable with as many values as the predictors
          n <- length(l_vars[[1L]])
          l_vars <- Filter(function(v) is.atomic(v) && length(v) == n, l)
        }

        list2DF(l_vars)
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
  arg::arg_string(treat)

  # The effect curve is a function whose first argument is named after the
  # treatment, so the treatment cannot share a name with its other argument or
  # with the objects it uses internally
  reserved <- c("subset", rlang::fn_fmls_names(.make_fun))

  if (anyv(reserved, treat)) {
    arg::err("the treatment cannot be named {.val {treat}}; rename it before fitting the model")
  }

  if (!anyv(names(data), treat)) {
    arg::err(c("{.arg treat} must be the name of a variable in the original dataset.",
               "*" = "Supplied value: {.val {treat}}",
               "*" = "Available names: {.or {.val {names(data)}}}"))
  }

  treat_var <- data[[treat]]

  if (!is.numeric(treat_var) || !is_null(dim(treat_var))) {
    arg::err("{.arg treat} must be the name of a numeric variable in the original dataset")
  }

  treat_var
}

process_wts <- function(wts, model, data = NULL) {
  if (is_null(wts)) {
    if (endsWith(class(model)[1L], "_weightit") &&
        utils::hasName(model.frame(model), "(s.weights)")) {
      wts <- "(s.weights)"
    }
    else {
      wts <- inherits(model, "svyglm")
    }
  }
  else if (!rlang::is_string(wts) && !rlang::is_bool(wts) && !is.numeric(wts)) {
    arg::err("{.arg wts} must be {.val {list(NULL)}}, {.val {TRUE}}, {.val {FALSE}}, a numeric vector, or a string containing the name of the variable containing weights in the original dataset")
  }

  if (isTRUE(wts)) {
    # insight::get_weights() returns NULL for svyglm models fit to replicate-weight
    # designs; prior.weights holds the (normalized) sampling weights for every svyglm
    wts <- {
      if (inherits(model, "svyglm")) model[["prior.weights"]]
      else insight::get_weights(model)
    } %or% FALSE
  }

  if (rlang::is_string(wts)) {
    if (wts %in% c("(weights)", "(s.weights)") && !is_null(model[["model"]]) &&
        utils::hasName(model[["model"]], wts)) {
      data <- model[["model"]]
    }
    else if (is_null(data)) {
      data <- try(insight::get_data(model, verbose = FALSE),
                  silent = TRUE)
    }

    if (!is.data.frame(data)) {
      arg::err("no dataset could be found to extract the weights named in {.arg wts}. Please supply a dataset to {.arg data}")
    }

    if (!utils::hasName(data, wts)) {
      arg::err("when supplied as a string, {.arg wts} must be the name of the variable in the original dataset")
    }

    wts <- data[[wts]]

    if (!is.numeric(wts)) {
      arg::err("when supplied as a string, {.arg wts} must be the name of a numeric variable")
    }

    return(.check_wts(as.vector(wts)))
  }

  # The rows of the evaluation data are the units used in the fit, including
  # any with zero model weights, which insight::n_obs() leaves out
  n <- {
    if (is.data.frame(data)) fnrow(data)
    else insight::n_obs(model)
  }

  if (is_null(n)) {
    if (is_null(data)) {
      data <- try(insight::get_data(model, verbose = FALSE),
                  silent = TRUE)

      if (!is.data.frame(data)) {
        arg::err("no dataset could be found to process the weights supplied to {.arg wts}. Please supply a dataset to {.arg data}")
      }
    }
    else {
      arg::arg_data.frame(data)
    }

    n <- fnrow(data)
  }

  if (isFALSE(wts)) {
    return(alloc(1.0, n))
  }

  if (length(wts) != n) {
    arg::err("when supplied as a numeric vector, {.arg wts} must have as many values as the number of observations in the original model")
  }

  .check_wts(as.vector(wts))
}

.check_wts <- function(wts) {
  if (!all(is.finite(wts))) {
    arg::err("{.arg wts} cannot contain missing or non-finite values")
  }

  if (any(wts < 0)) {
    arg::err("{.arg wts} cannot contain negative values")
  }

  if (sum(wts) <= 0) {
    arg::err("at least one value of {.arg wts} must be positive")
  }

  wts
}

# Target weights for each model: one model, or the models fit to each imputed
# dataset. Numeric weights for multiply imputed data may have one value per row
# of the stacked imputed datasets or one value per row of each imputed dataset.
process_wts_list <- function(wts, models, eval_data.list) {
  m <- length(models)

  # The design-based variance of svyglm models linearizes the estimates over
  # the population the design represents, so the target weights must be the
  # design weights
  if (!is_null(wts) && !isTRUE(wts) && inherits(models[[1L]], "svyglm")) {
    arg::err("{.arg wts} cannot be supplied with {.cls svyglm} models, for which the survey weights are always used")
  }

  if (m == 1L) {
    return(list(process_wts(wts, models[[1L]], eval_data.list[[1L]])))
  }

  if (!is.numeric(wts)) {
    return(lapply(seq_len(m), function(.i) {
      process_wts(wts, models[[.i]], eval_data.list[[.i]])
    }))
  }

  ns <- vapply(eval_data.list, fnrow, integer(1L))

  if (length(wts) == sum(ns)) {
    return(as.vector(wts) |>
             split(rep.int(seq_len(m), ns)) |>
             unname())
  }

  if (!all_the_same(ns)) {
    arg::err("when {.arg wts} is supplied as a numeric vector with multiply imputed data, it must have length equal to the sum of the number of observations across all imputed datasets, or the imputed datasets must have the same number of observations")
  }

  if (length(wts) != ns[1L]) {
    arg::err("when {.arg wts} is supplied as a numeric vector with multiply imputed data, it must have length equal to the sum of the number of observations across all imputed datasets or to the number of observations in each imputed dataset")
  }

  rep.int(list(as.vector(wts)), m)
}

# `null` is a numeric null value to be tested on the transformed scale, named
# `null_arg` in messages, and `check_draws` is whether the bootstrap or
# posterior draws will be transformed (they are not for percentile intervals).
# A transformation requested explicitly that cannot be applied to every value
# it will be applied to gives an error, as does a null value outside its
# domain; the link function of the outcome model (`transform = TRUE`) is
# replaced by the identity with a warning instead, since the user did not
# choose it.
process_transform <- function(transform, x = NULL, .est = NULL, null = NULL,
                              null_arg = "null", check_draws = TRUE) {
  if (inherits(x, "effect_curve")) {
    if (!.is_pure_adrf(x)) {
      transform <- FALSE
    }
  }
  else if (!identical(get_curve_type(x), "ADRF") ||
           !is_null(.attr(x, ".contrast")) ||
           !is_null(.attr(x, ".reference"))) {
    transform <- FALSE
  }

  if (isFALSE(transform)) {
    .family <- stats::quasi("identity")
  }
  else if (isTRUE(transform)) {
    .family <- .attr(x, ".family") %or% stats::quasi("identity")
  }
  else if (inherits(transform, c("family", "link-glm")) &&
           all(c("linkfun", "linkinv", "mu.eta") %in% names(transform))) {
    .family <- transform
  }
  else if (is.function(transform) && !is_null(.known_transform(transform))) {
    .family <- .known_transform(transform)
  }
  else {
    arg::err("{.arg transform} must be {.val {TRUE}}, {.val {FALSE}}, a family, a link object (e.g., from {.fun make.link}), or one of {.fun log}, {.fun sqrt}, {.fun exp}, {.fun qlogis}, or {.fun qnorm}")
  }

  if (!is_null(.est) && !.is_identity_family(.family)) {
    .family <- .check_transform_values(.family, x, .est, null, null_arg,
                                       check_draws, explicit = !isTRUE(transform))
  }

  linkfun <- .family$linkfun
  mu.eta <- .family$mu.eta

  list(transform = linkfun,
       inv_transform = .family$linkinv,
       d_transform = function(mu) 1 / mu.eta(linkfun(mu)),
       d_inv_transform = mu.eta)
}

.is_identity_family <- function(.family) {
  identical(.family[["link"]] %or% .family[["name"]], "identity")
}

# Checks that the transformation can be applied to and reversed for the null
# value, the estimates, and, when they will be transformed, the bootstrap or
# posterior draws, which between grid points can fall outside the domain of
# the transformation. Returns the family to use.
.check_transform_values <- function(.family, x, .est, null, null_arg,
                                    check_draws, explicit) {
  null <- null[!is.na(null)]

  if (!is_null(null) && !.transform_is_valid(.family, null)) {
    arg::err(c("The null value supplied to {.arg {null_arg}} ({null}) is outside the domain of the transformation, so it cannot be tested on the transformed scale.",
               i = "Set {.code transform = FALSE} to perform the test on the original scale."))
  }

  draws <- {
    if (check_draws) .attr(x, ".boot")[["t"]] %or% .attr(x, ".draws")
  }

  what <- {
    if (!.transform_is_valid(.family, .est)) "estimates"
    else if (!is_null(draws) && !.transform_is_valid(.family, draws)) {
      if (is_null(.attr(x, ".boot"))) "posterior draws"
      else "bootstrap estimates"
    }
  }

  if (is_null(what)) {
    return(.family)
  }

  if (explicit) {
    arg::err(c("The transformation supplied to {.arg transform} cannot be applied to all of the {what}, some of which are outside its domain.",
               i = "Supply a different transformation or set {.code transform = FALSE}{if (what != \"estimates\") \", or use percentile intervals where available\" else \"\"}."))
  }

  arg::wrn("some of the {what} are outside the domain of the link function of the outcome model, so they are not transformed. Set {.code transform = FALSE} to silence this warning")

  stats::quasi("identity")
}

# Families for the transformations that can be supplied as functions; each has
# a closed-form inverse and derivative
.known_transform <- function(transform) {
  if (identical(transform, log)) {
    stats::quasi("log")
  }
  else if (identical(transform, sqrt)) {
    stats::quasi("sqrt")
  }
  else if (identical(transform, qlogis)) {
    stats::quasi("logit")
  }
  else if (identical(transform, qnorm)) {
    stats::quasi("probit")
  }
  else if (identical(transform, exp)) {
    list(linkfun = exp,
         linkinv = log,
         mu.eta = function(eta) 1 / eta)
  }
}

.transform_is_valid <- function(.family, .est) {
  if (!all(is.finite(.est))) {
    return(FALSE)
  }

  .t_est <- try(suppressWarnings(.family$linkfun(.est)), silent = TRUE)

  if (null_or_error(.t_est) || length(.t_est) != length(.est) ||
      !all(is.finite(.t_est))) {
    return(FALSE)
  }

  .back <- try(suppressWarnings(.family$linkinv(.t_est)), silent = TRUE)

  if (null_or_error(.back) || length(.back) != length(.est) ||
      !all(is.finite(.back))) {
    return(FALSE)
  }

  # Tolerance relative to the scale of the estimates
  tol <- sqrt(.Machine$double.eps) * max(abs(.est), .Machine$double.eps)

  all(abs(.back - .est) <= tol)
}

process_family <- function(model, type = "response") {
  if (type %!in% c("response", "probs")) {
    .family <- stats::quasi("identity")
  }
  else if (inherits(model, "bart")) {
    .family <- {
      if (is_null(model[["sigma"]])) stats::binomial("probit")
      else stats::gaussian()
    }
  }
  else {
    .family <- try(insight::get_family(model),
                   silent = TRUE)

    if (!inherits(.family, "family")) {
      .family <- stats::quasi("identity")
    }
  }

  # Some families (e.g., from bartisan) store only the name of the link;
  # rebuild the link functions from it when possible
  if (!all(c("linkfun", "linkinv", "mu.eta") %in% names(.family))) {
    link_obj <- try(stats::make.link(.family[["link"]]), silent = TRUE)

    if (inherits(link_obj, "link-glm")) {
      .family[c("linkfun", "linkinv", "mu.eta")] <- link_obj[c("linkfun", "linkinv", "mu.eta")]
    }
    else {
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

# The influence of each unit on the coefficients, estfun() %*% t(bread()),
# used for the unconditional variance, or NULL when the model has no estfun()
# or bread() method. For svyglm models, it is computed from the survey design
# and scaled by the full sample size.
.coef_influence <- function(model, hc2 = TRUE, cluster = NULL) {
  if (inherits(model, "svyglm")) {
    X <- stats::model.matrix(model)
    infl <- .attr(model, "influence")

    if (is_null(infl)) {
      infl <- tcrossprod(X * stats::resid(model, "working") * model$weights,
                         model$naive.cov)
    }
    else if (!is_null(model[["na.action"]]) &&
             fnrow(infl) == nrow(X) + length(model[["na.action"]])) {
      # svyglm(influence = TRUE) pads the units dropped for missing values with
      # rows of 0
      infl <- infl[-model[["na.action"]], , drop = FALSE]
    }

    return(fnrow(infl) * infl)
  }

  eb <- {
    if (inherits(model, "gam")) .gam_estfun_bread(model)
    else .sandwich_estfun_bread(model)
  }

  if (is_null(eb)) {
    return(NULL)
  }

  ef <- eb[["ef"]]
  br <- eb[["br"]]

  infl <- tcrossprod(ef, br)

  if (hc2) {
    infl <- .hc2_adjust(infl, model)
  }

  if (!is_null(cluster)) {
    cr2 <- lapply(cluster, function(f) .cr2_estfun(model, ef, f))

    if (!any_apply(cr2, is_null)) {
      attr(infl, "cr2") <- lapply(cr2, tcrossprod, br)
    }
  }

  infl
}

# estfun() and bread() from sandwich, as a list, or NULL when either fails
.sandwich_estfun_bread <- function(model) {
  ef <- try({
    # A WeightIt fit with `vcov = "HC0"` treats the weights as fixed, so the
    # estimation of the weights is left out of the influence function
    if (inherits(model, "glm_weightit") && identical(model[["vcov_type"]], "HC0"))
      sandwich::estfun(model, asympt = FALSE)
    else
      sandwich::estfun(model)
  }, silent = TRUE)

  if (null_or_error(ef)) {
    return(NULL)
  }

  # With na.action = na.exclude, estfun() pads the units dropped from the fit
  # with NA rows; the influence values are for the units used in the fit
  na_act <- model[["na.action"]]

  if (inherits(na_act, "exclude") && fnrow(ef) > length(na_act)) {
    ef <- ef[-na_act, , drop = FALSE]
  }

  br <- try(sandwich::bread(model), silent = TRUE)

  if (null_or_error(br)) {
    return(NULL)
  }

  list(ef = ef, br = br)
}

# The equivalents of estfun() and bread() for mgcv's gam() and bam() fits, as
# a list, or NULL when they cannot be computed. In sandwich 3.1-3 and
# earlier, the methods do not agree for these fits: bread() has a gam method,
# which omits the dispersion, but estfun() falls back to the glm method, which
# divides by it, so with an estimated scale the standard errors are too small
# by a factor of the dispersion. This is fixed in the development version
# (3.1-3.99), with which the two give the same influence values for the
# families supported here; this function can be replaced by
# .sandwich_estfun_bread() once adrftools requires the fixed version.
# mgcv::vcov.gam(sandwich = TRUE) cannot be used instead
# because it omits the prior weights. The penalized coefficients solve
# X'W(z - X beta) - S beta = 0, so the influence of unit i on them is
# (X'WX + S)^(-1) x_i w_i r_i = Vp x_i w_i r_i / scale, with w_i and r_i the
# working weights (which include the prior weights) and residuals. The
# uncertainty in the smoothing parameters is not included. Extended families
# (e.g., nb(), tw(), ocat()) estimate other parameters (theta) jointly with the coefficients,
# which this leaves out, so they are not supported.
.gam_estfun_bread <- function(model) {
  if (!inherits(model[["family"]], "family") ||
      inherits(model[["family"]], c("extended.family", "general.family"))) {
    return(NULL)
  }

  X <- try(stats::model.matrix(model), silent = TRUE)

  if (null_or_error(X)) {
    return(NULL)
  }

  Vp <- model[["Vp"]]
  scale <- model[["sig2"]]
  w <- model[["weights"]]
  r <- model[["residuals"]]
  n <- fnrow(X)

  if (!is.matrix(Vp) || fncol(X) != fncol(Vp) || length(w) != n || length(r) != n ||
      !is_number(scale) || !is.finite(scale) || scale <= 0) {
    return(NULL)
  }

  ef <- X * (w * r / scale)
  attr(ef, "lpi") <- NULL
  attr(ef, "model.offset") <- NULL

  # Scaled like sandwich::bread(), by the number of units with nonzero weight
  n_bread <- try(stats::nobs(model), silent = TRUE)

  if (null_or_error(n_bread) || !is_number(n_bread) || !is.finite(n_bread)) {
    n_bread <- sum(w > 0)
  }

  list(ef = ef, br = n_bread * Vp)
}

# The CR2 adjustment for the clustering factor `f`, as described by
# Pustejovsky and Tipton (2018): within each cluster g, the vector of
# (weighted) residuals u_g is replaced by (I - H_gg)^(-1/2) u_g, with H_gg the
# block of the hat matrix for the cluster, which makes the cluster-robust
# variance unbiased when the working model's errors are correct. For GLMs,
# the residuals and hat matrix are those of the final weighted least squares
# fit; sandwich::vcovCL(type = "HC2") instead applies the nonsymmetric
# (I - X_g (X'WX)^(-1) X_g' W_g)^(-1/2) to W_g r_g, which differs when the
# working weights vary within a cluster. Returns the adjusted rows of estfun(), whose
# sums within clusters are the adjusted cluster scores, or NULL when the model
# is not an lm() or glm() fit whose estfun() is the usual product of the
# weighted residuals and model matrix (e.g., a WeightIt fit that accounts for
# the estimation of the weights). With clusters of one unit, it is the HC2
# adjustment.
.cr2_estfun <- function(model, ef, f) {
  q <- model[["qr"]]

  if (!inherits(model, "lm") || !inherits(q, "qr")) {
    return(NULL)
  }

  n <- fnrow(ef)

  w <- model[["weights"]] %or% rep.int(1, n)
  res <- model[["residuals"]]
  X <- try(stats::model.matrix(model), silent = TRUE)

  if (null_or_error(X) || length(w) != n || length(res) != n || fnrow(X) != n ||
      length(f) != n) {
    return(NULL)
  }

  aliased <- is.na(stats::coef(model))
  X <- X[, !aliased, drop = FALSE]

  pos <- which(w > 0)

  if (length(pos) != nrow(q[["qr"]]) || fncol(X) != fncol(ef)) {
    return(NULL)
  }

  # estfun() should be a scalar multiple (e.g., 1 / dispersion) of z * u, with
  # z = sqrt(w) x and u = sqrt(w) r for working residuals r
  sw <- sqrt(w[pos])
  Z <- X[pos, , drop = FALSE] * sw
  u <- sw * res[pos]
  ef_std <- Z * u

  k <- sum(ef[pos, ] * ef_std) / sum(ef_std^2)

  tol <- sqrt(.Machine$double.eps) * max(abs(ef), .Machine$double.xmin)

  if (!is.finite(k) || any(abs(ef[pos, , drop = FALSE] - k * ef_std) > tol) ||
      any(abs(ef[-pos, , drop = FALSE]) > tol)) {
    return(NULL)
  }

  # Rows of the orthonormal factor of the weighted model matrix, so that
  # H_gg = Q_g Q_g'
  Q <- qr.qy(q, diag(1, length(pos), q[["rank"]]))

  # With the singular value decomposition Q_g = U D V', I - H_gg = I - U D^2 U',
  # so (I - H_gg)^(-1/2) u_g = u_g + U ((1 - d^2)^(-1/2) - 1) U'u_g; a direction
  # with leverage 1, in which the residual is 0, is dropped
  Au <- u

  for (idx in split(seq_along(pos), f[pos], drop = TRUE)) {
    sv <- svd(Q[idx, , drop = FALSE], nv = 0L)
    d2 <- sv$d^2
    a <- rep.int(-1, length(d2))
    ok <- d2 < 1 - sqrt(.Machine$double.eps)
    a[ok] <- 1 / sqrt(1 - d2[ok]) - 1

    Au[idx] <- u[idx] + drop(sv$u %*% (a * crossprod(sv$u, u[idx])))
  }

  out <- ef
  out[pos, ] <- k * Z * Au
  out[-pos, ] <- 0

  out
}

# The HC2 adjustment: each unit's influence on the coefficients is divided by
# sqrt(1 - h), with h its leverage. A residual is shrunk toward 0 by the fit,
# with E[e^2] = (1 - h) sigma^2 under homoskedasticity, so the HC0 variance is
# biased downward, most for the contrasts that depend on units with high
# leverage (e.g., an ADRF at the extremes of the treatment); the HC2
# adjustment removes this bias under homoskedasticity (MacKinnon and White,
# 1985). When the leverages cannot be found, the HC1 adjustment is used
# instead, which removes the bias for the average leverage, p / n, and the
# influence values get the attribute `hc2_failed`.
.hc2_adjust <- function(infl, model) {
  h <- .hat_values(model, fnrow(infl))

  if (is_null(h)) {
    infl <- infl * sqrt(.hc1_factor(model, infl))
    attr(infl, "hc2_failed") <- TRUE
    return(infl)
  }

  # A unit with leverage 1 has a residual of 0, so its row is left unscaled
  # rather than divided by 0
  adj <- rep.int(1, length(h))
  ok <- h < 1 - sqrt(.Machine$double.eps)
  adj[ok] <- 1 / sqrt(1 - h[ok])

  infl * adj
}

# The HC1 factor n / (n - p), with n the number of units with nonzero weight
# and p the number of estimated parameters (columns of estfun(), including
# ancillary parameters such as thresholds, or the effective degrees of
# freedom of a penalized fit); the leverages of any model sum to p, so this
# corrects the bias for the average leverage. It is 1 when n <= p.
.hc1_factor <- function(model, infl) {
  n <- .n_fit_units(model, infl)
  p <- .n_params(model, infl)

  if (n <= p) 1
  else n / (n - p)
}

# The CR1 factor (n - 1) / (n - p) that multiplies G / (G - 1), as in
# sandwich::vcovCL() with its default type. It is 1 when n <= p.
.cr1_factor <- function(model, infl) {
  n <- .n_fit_units(model, infl)
  p <- .n_params(model, infl)

  if (n <= p) 1
  else (n - 1) / (n - p)
}

# The number of estimated parameters: the columns of estfun(), or the
# effective degrees of freedom of a gam() or bam() fit, whose penalized
# coefficients are not free parameters
.n_params <- function(model, infl) {
  if (inherits(model, "gam") && is.numeric(model[["edf"]]) &&
      length(model[["edf"]]) == fncol(infl)) {
    return(sum(model[["edf"]]))
  }

  fncol(infl)
}

# The number of units used to fit the model with nonzero weight
.n_fit_units <- function(model, infl) {
  n <- try(stats::nobs(model), silent = TRUE)

  if (null_or_error(n) || !is_number(n) || !is.finite(n)) {
    n <- sum(rowSums(abs(infl)) > 0)
  }

  n
}

# Leverages (hat values) of the units, one per row of estfun(), or NULL when
# they cannot be found. For lm() and glm() fits (and models inheriting from
# them), they are computed from the QR decomposition stored in the fit, which
# covers only the units with positive (working) weights; the others get a
# leverage of 0. gam() fits store their leverages. For other models,
# hatvalues() is used.
.hat_values <- function(model, n) {
  q <- model[["qr"]]

  if (inherits(model, "lm") && inherits(q, "qr")) {
    r <- q[["rank"]]
    h0 <- rowSums(qr.qy(q, diag(1, nrow(q[["qr"]]), r))^2)

    w <- model[["weights"]]

    if (is_null(w)) {
      w <- rep.int(1, n)
    }

    pos <- which(w > 0)

    if (length(w) == n && length(pos) == length(h0)) {
      h <- numeric(n)
      h[pos] <- h0
      return(h)
    }
  }

  # gam() fits store the leverages of the penalized fit, which sum to the
  # effective degrees of freedom; bam() fits store something else in `hat`
  if (inherits(model, "gam") && !inherits(model, "bam")) {
    h <- model[["hat"]]

    if (is.numeric(h) && length(h) == n && !anyNA(h) && all(h >= 0) &&
        all(h <= 1 + sqrt(.Machine$double.eps))) {
      return(unname(h))
    }
  }

  h <- try(stats::hatvalues(model), silent = TRUE)

  if (null_or_error(h) || !is.numeric(h) || length(h) != n ||
      anyNA(h) || any(h < 0) || any(h > 1 + sqrt(.Machine$double.eps))) {
    return(NULL)
  }

  unname(h)
}

process_vcov_and_cluster <- function(vcov, model, cluster = NULL, is_bayes = FALSE,
                                     data = NULL) {
  if (is_null(vcov) || isFALSE(vcov) || identical(vcov, "none")) {
    if (!is_null(cluster)) {
      .v <- if (is_null(vcov)) list(NULL) else vcov

      arg::wrn("{.arg cluster} is ignored when {.code vcov = {.val {(.v)}}}")
    }

    return(list(vcov = "none"))
  }

  if (is_bayes) {
    if (identical(vcov, "unconditional")) {
      vcov <- "posterior"
    }
    else {
      builtin_vcovs <- c("none", "posterior")

      vcov <- arg::match_arg(vcov, builtin_vcovs, .context = "with Bayesian models,")

      if (vcov == "none") {
        return(Recall(vcov = "none", model = model, cluster = cluster,
                      is_bayes = is_bayes, data = data))
      }
    }

    if (!is_null(cluster)) {
      arg::wrn("{.arg cluster} is ignored with Bayesian models")
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
        if (inherits(model, "svyglm")) {
          arg::err(c("Bootstrapping cannot be used with {.cls svyglm} models because resampling units ignores the survey design.",
                     i = "Use {.code vcov = \"unconditional\"}, which uses the survey design, or fit the model to a replicate-weight design."))
        }

        rlang::check_installed("fwb (>= 1.0.0)")

        if (!is_null(cluster)) {
          cl <- process_cluster(cluster, model, data)

          if (length(cl[["cluster"]]) > 1L) {
            arg::err(c("Only one clustering variable can be supplied to {.arg cluster} with {.code vcov = {.val {builtin_vcovs[p]}}}.",
                       i = "For multiway clustering, use {.code vcov = \"unconditional\"} with {.arg cluster}, or supply a multiway cluster bootstrap covariance matrix of the model coefficients (e.g., {.code vcov = fwb::vcovFWB(fit, cluster = ~a + b)}), which uses Wald inference."))
          }

          return(c(list(vcov = builtin_vcovs[p]),
                   cl))
        }

        return(list(vcov = builtin_vcovs[p]))
      }

      if (builtin_vcovs[p] == "unconditional") {
        if (inherits(model, "svyglm") && !is_null(cluster)) {
          arg::wrn("{.arg cluster} is ignored with {.cls svyglm} models when {.code vcov = \"unconditional\"}; the survey design is used instead")
          cluster <- NULL
        }

        cl <- process_cluster(cluster, model, data)

        # Computed here, where it also checks that the model is supported, and
        # passed on so that estfun() and bread() are not computed twice. The
        # HC2 adjustment is used without clustering and the CR2 adjustment
        # with it (survey designs are handled in .coef_influence())
        infl <- .coef_influence(model, hc2 = is_null(cl), cluster = cl[["cluster"]])

        if (is_null(infl)) {
          arg::wrn('{.arg vcov} cannot be {.val unconditional} with this type of model. Using {.val conditional}')

          return(Recall(vcov = "conditional", model = model, cluster = cluster,
                        is_bayes = is_bayes, data = data))
        }

        if (isTRUE(.attr(infl, "hc2_failed"))) {
          arg::msg("the leverages of the units could not be found for this type of model, so the unconditional variance uses the HC1 small-sample adjustment in place of HC2")
        }

        if (!is_null(cl)) {
          if (is_null(.attr(infl, "cr2"))) {
            arg::msg("the CR2 small-sample adjustment is available only for {.fun lm} and {.fun glm} models, so the cluster-robust unconditional variance uses the CR1 adjustment in place of CR2")

            # CR1: G / (G - 1), already in `adj`, times (n - 1) / (n - p)
            cl[["adj"]] <- cl[["adj"]] * .cr1_factor(model, infl)
          }
          else {
            # One adjusted influence matrix per clustering term; the CR2
            # adjustment replaces the G / (G - 1) factor
            cl[["infl_cl"]] <- .attr(infl, "cr2")
            cl[["adj"]] <- sign(cl[["adj"]])
          }
        }

        attr(infl, "hc2_failed") <- NULL
        attr(infl, "cr2") <- NULL

        return(c(list(vcov = builtin_vcovs[p], infl = infl), cl))
      }

      if (!is_null(cluster)) {
        # The default type of sandwich::vcovCL()
        return(.cluster_vcov(model, cluster, data, type = NULL))
      }

      # Use default get_vcov() below
      vcov <- NULL
    }
    else if (!is_null(cluster)) {
      cl_types <- c("HC0", "HC1", "HC2", "HC3")

      if (vcov %!in% cl_types) {
        arg::err(c("{.arg vcov} cannot be {.val {vcov}} when {.arg cluster} is supplied.",
                   i = "With {.arg cluster}, {.arg vcov} can be {.or {.val {c(builtin_vcovs[2:3], cl_types)}}}, or a cluster-robust covariance matrix of the model coefficients."))
      }

      return(.cluster_vcov(model, cluster, data, type = vcov))
    }
  }

  if ((is.function(vcov) || rlang::is_formula(vcov)) && !is_null(cluster)) {
    arg::err(c("{.arg cluster} cannot be supplied when {.arg vcov} is a function or formula.",
               i = "Supply the clustering variables to {.arg cluster} with a string for {.arg vcov}, or include them in the function or formula."))
  }

  if (is.numeric(vcov) && !is_null(cluster)) {
    arg::wrn("{.arg cluster} is ignored when {.arg vcov} is supplied as a numeric matrix")
    cluster <- NULL
  }

  vcov_try <- try(marginaleffects::get_vcov(model, vcov = cluster %or% vcov),
                  silent = TRUE)

  if (null_or_error(vcov_try)) {
    arg::err("{.arg vcov} must be either one of {.or {.val {builtin_vcovs}}} or one of the allowed arguments to {.fun marginaleffects::get_vcov}")
  }

  c(list(vcov = vcov_try),
    process_cluster(cluster, model, data))
}

# Cluster-robust covariance matrix of the model coefficients from
# sandwich::vcovCL(), using the clustering variables as processed by
# process_cluster(), so formulas with transformations and data with missing
# values work. `type = NULL` uses vcovCL()'s default type.
# Sets negative eigenvalues of a symmetric matrix to 0, as
# sandwich::vcovCL(fix = TRUE) does for multiway clustering
.clip_psd <- function(V) {
  e <- eigen(V, symmetric = TRUE)

  if (all(e$values >= 0)) {
    return(V)
  }

  tcrossprod(e$vectors %*% diag(pmax(e$values, 0), nrow = length(e$values)), e$vectors)
}

.cluster_vcov <- function(model, cluster, data, type = NULL) {
  cl <- process_cluster(cluster, model, data)

  # With multiway clustering, negative eigenvalues are set to 0
  vcov_try <- try(sandwich::vcovCL(model, type = type, cluster = cl[["vars"]],
                                   fix = fncol(cl[["vars"]]) > 1L),
                  silent = TRUE)

  if (null_or_error(vcov_try)) {
    .type <- type %or% "conditional"

    arg::err(c("The cluster-robust covariance matrix could not be computed with {.code vcov = {.val {(.type)}}}.",
               x = conditionMessage(attr(vcov_try, "condition")),
               i = if (.type %in% c("HC2", "HC3")) "{.val {(.type)}} with {.arg cluster} is available only for {.cls lm} and {.cls glm} models."))
  }

  c(list(vcov = vcov_try), cl)
}

# The covariance specification for each model: one model, or the models fit to
# each imputed dataset. With multiply imputed data, `vcov` may be a list with
# one entry per imputed dataset, and bootstrapping is not allowed.
process_vcov_and_cluster_list <- function(vcov, models, cluster = NULL, is_bayes = FALSE,
                                          data_list) {
  m <- length(models)

  if (m == 1L) {
    return(list(process_vcov_and_cluster(vcov, models[[1L]], cluster,
                                         is_bayes, data_list[[1L]])))
  }

  if (rlang::is_string(vcov) && !is_bayes) {
    builtin_vcovs <- c("none", "unconditional", "conditional", "boot", "fwb")
    p <- pmatch(vcov, builtin_vcovs, nomatch = 0L)

    if (p %in% 4:5) {
      arg::err("{.arg vcov} cannot be {.val {builtin_vcovs[p]}} with multiply imputed data")
    }
  }

  vcov.list <- {
    if (is.list(vcov)) vcov
    else rep.int(list(vcov), m)
  }

  if (length(vcov.list) != m) {
    arg::err("{.arg vcov} must have length equal to the number of imputed datasets")
  }

  # Warnings and messages would be the same for every imputed dataset, so
  # they are given only for the first
  lapply(seq_len(m), function(.i) {
    out <- function() {
      process_vcov_and_cluster(vcov.list[[.i]], models[[.i]], cluster,
                               is_bayes, data_list[[.i]])
    }

    if (.i == 1L) out()
    else suppressMessages(suppressWarnings(out()))
  })
}

process_cluster <- function(.cluster, model, data) {
  if (is_null(.cluster)) {
    return(NULL)
  }

  if (rlang::is_formula(.cluster)) {
    # Retrieve the variables named in the formula for exactly the rows used in
    # the fit (na.expand = TRUE), then evaluate its terms (e.g.,
    # interaction(a, b)) using them
    cluster <- try({
      model |>
        expand.model.frame(extras = reformulate(all.vars(.cluster)),
                           na.expand = TRUE) |>
        model.frame(formula = .cluster, na.action = na.pass)
    }, silent = TRUE)

    if (null_or_error(cluster)) {
      cluster <- try(model.frame(.cluster, data = data, na.action = na.pass),
                     silent = TRUE)
    }

    if (null_or_error(cluster)) {
      arg::err("clusters could not be extracted from the model; try supplying {.arg cluster} as a data frame or list")
    }
  }
  else {
    cluster <- qDF(.cluster)
  }

  if (fnrow(cluster) != fnrow(data)) {
    arg::err("the number of observations in {.arg cluster} must equal that in the original dataset")
  }

  if (anyNA(cluster, recursive = TRUE)) {
    arg::err("{.val {NA}}s are not allowed in the {cli::qty(fncol(cluster))} clustering variable{?s}")
  }

  # Levels with no units (e.g., because of missing values or `subset`) would
  # count as clusters
  cluster[] <- lapply(cluster, qF, drop = TRUE)

  if (any(vapply(cluster, fnlevels, numeric(1L)) < 2)) {
    arg::err("each clustering variable must have at least 2 clusters among the units used to fit the model")
  }

  # The clustering variables themselves, for sandwich::vcovCL()
  cluster_vars <- cluster

  p <- length(cluster)
  if (p > 1L) {
    clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE)) |>
      unlist(recursive = FALSE)
    sgn <- (-1L)^(vlengths(clu) + 1L)
    cluster <- lapply(clu, function(i) fdroplevels(finteraction(cluster[i])))
  }
  else {
    sgn <- 1
  }

  #For small sample adjustment (setting cadjust = TRUE in vcovCL)
  g <- vapply(cluster, fnlevels, numeric(1L))

  list(cluster = cluster,
       adj = sgn * g / (g - 1),
       vars = cluster_vars)
}

# The levels of the factor and character predictors used to fit the model, so
# that model.matrix() keeps the columns of levels absent from a subset of the
# data (as it already does for factors)
.model_xlev <- function(model) {
  xlev <- try(model[["xlevels"]], silent = TRUE)

  if (is.list(xlev)) xlev
  else NULL
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
  arg::arg_number(conf_level)
  arg::arg_between(conf_level, c(0, 1), inclusive = c(TRUE, FALSE))

  conf_level
}

process_ci.type <- function(ci.type, .vcov_type, simultaneous = NULL) {
  if (.vcov_type == "none") {
    return(NULL)
  }

  if (.vcov_type == "posterior") {
    ci.type <- arg::match_arg(ci.type, c("perc", "wald"),
                              .context = "with Bayesian models,")

    return(ci.type)
  }

  if (.vcov_type != "bootstrap") {
    return("wald")
  }

  arg::arg_string(ci.type)

  if (isTRUE(simultaneous)) {
    ci.type <- arg::match_arg(ci.type, c("perc", "wald"),
                              .context = "when {.code simultaneous = TRUE},")
  }

  ci.type
}

process_simultaneous <- function(simultaneous, .est) {
  if (length(.est) <= 1) {
    return(FALSE)
  }

  arg::arg_flag(simultaneous)

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
  if (!is_null(null)) {
    if (length(null) == 1L && allNA(null)) {
      return(NA_real_)
    }

    arg::arg_number(null)
    return(null)
  }

  if (inherits(x, "effect_curve")) {
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
  if (!is_null(proj)) {
    arg::arg_is(proj, "curve_projection")

    if (any_apply(c(".values", ".treat", ".vcov_type", ".curve_type", ".response",
                    ".by_grid", ".contrast", ".reference", ".family", ".df", ".shift"),
                  function(a) !identical(.attr(x, a), .attr(proj, a)))) {
      arg::err("{.arg proj} must be the output of {.fun curve_projection} applied to the {.cls effect_curve} object")
    }
  }
}

get_curve_type <- function(x) {
  .curve_type <- .attr(x, ".curve_type")

  if (!is_null(.curve_type)) {
    return(.curve_type)
  }

  if (inherits(x, "adrf_curve")) {
    return("ADRF")
  }

  if (inherits(x, "amef_curve")) {
    return("AMEF")
  }

  NULL
}

# The name of the curve type in printed output and plot titles. A shift curve
# has the curve type of an ADRF (its estimates are mean outcomes), and its
# derivative that of an AMEF, but they are described as shift curves.
.curve_type_label <- function(x, title = FALSE) {
  .curve_type <- get_curve_type(x)

  if (is_null(.attr(x, ".shift"))) {
    return(.curve_type)
  }

  label <- switch(.curve_type,
                  ADRF = "shift curve",
                  AMEF = "shift curve derivative",
                  return(.curve_type))

  if (title) {
    label <- gsub("\\b([a-z])", "\\U\\1", label, perl = TRUE)
  }

  label
}

check_effect_curve <- function(x, amef_ok = TRUE, reference_ok = TRUE,
                               contrast_ok = TRUE, projection_ok = TRUE,
                               nm = rlang::caller_arg(x)) {

  arg::arg_is(x, "effect_curve")

  if (!amef_ok && inherits(x, "amef_curve")) {
    arg::err("{.arg {nm}} cannot be an {.cls amef_curve} object")
  }

  if (!reference_ok && inherits(x, "reference_curve")) {
    arg::err("{.arg {nm}} cannot be a {.cls reference_curve} object")
  }

  if (!contrast_ok && inherits(x, "contrast_curve")) {
    arg::err("{.arg {nm}} cannot be a {.cls contrast_curve} object")
  }

  if (!projection_ok && inherits(x, "curve_projection")) {
    arg::err("{.arg {nm}} cannot be a {.cls curve_projection} object")
  }
}

process_bw <- function(bw, v, constant = .5) {
  if (is_null(bw)) {
    spacing <- NULL
    diff_v <- diff1(v)

    # Spacings are compared relative to the largest one, so that the check does
    # not depend on the scale of the treatment
    max_diff <- max(abs(diff_v))
    tol <- sqrt(.Machine$double.eps)

    if (max_diff > 0 && all(abs(diff_v - mean(diff_v)) <= tol * max_diff)) {
      spacing <- mean(diff_v)
    }
    else if (max_diff > 0) {
      spacing_table <- table(round(diff_v / max_diff, 8L))

      spacing <- as.numeric(names(spacing_table)[which.max(spacing_table)]) * max_diff

      if (mean(abs(diff_v - spacing) <= tol * max_diff) < .7) {
        spacing <- NULL
      }
    }

    if (is_null(spacing)) {
      return(bw.nrd(v))
    }

    arg::arg_number(constant)

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
    arg::err("{.arg bw} must be a string, function, or number")
  }

  arg::arg_number(bw)
  arg::arg_gt(bw, 0)

  bw
}

get_n_by <- function(.contrast = NULL, .by_grid = NULL) {
  if (!is_null(.contrast)) {
    return(length(.contrast))
  }

  if (!is_null(.by_grid)) {
    return(fnrow(.by_grid))
  }

  1L
}

# Returns a function that refits `x` with bootstrap weights `w` (one per unit
# used in the fit). The model is refit through its own call, on the rows used in
# the original fit, so that missing values, `subset`, and the environment the
# model was fit in do not matter. The model's own weights are multiplied by the
# bootstrap weights. For models fit with a weightit object, the weights are
# re-estimated in each replicate by refitting the weighting model with its
# sampling weights multiplied by the bootstrap weights, as WeightIt does.
.make_boot_refit <- function(x, data = NULL) {
  fit_data <- try(insight::get_data(x, verbose = FALSE), silent = TRUE)

  if (null_or_error(fit_data) || !is.data.frame(fit_data) ||
      !.data_matches_model(x, fit_data)) {
    if (is_null(data)) {
      arg::err("the dataset used to fit the model could not be recovered for bootstrapping. Supply it to {.arg data}")
    }

    arg::arg_data.frame(data)
    fit_data <- data
  }

  model_call <- stats::getCall(x)

  if (is_null(model_call)) {
    arg::err("bootstrapping requires a model object that stores the call used to fit it")
  }

  genv <- environment(stats::formula(x)) %or% globalenv()

  fn <- try(eval(model_call[[1L]], envir = genv), silent = TRUE)

  if (is.function(fn)) {
    model_call <- match.call(fn, model_call)
  }

  if (utils::hasName(as.list(model_call), "formula")) {
    model_call$formula <- stats::formula(x)
  }

  model_call$data <- quote(.adrftools_data)
  model_call$subset <- NULL

  weightit <- x[["weightit"]]

  if (!is_null(weightit)) {
    wcall <- weightit[["call"]]

    # weightit() stores its calling environment only for top-level calls
    wenv <- weightit[["env"]] %or% environment(weightit[["formula"]]) %or% globalenv()

    if (is_null(wcall) || !is.environment(wenv)) {
      arg::err("bootstrapping cannot be used because the {.cls weightit} object supplied to the model cannot be re-estimated")
    }

    sw <- weightit[["s.weights"]] %or% 1

    model_call$weightit <- quote(.adrftools_weightit)
    model_call$vcov <- "none"

    return(function(w) {
      wcall$s.weights <- sw * w

      weightit_boot <- suppressMessages(suppressWarnings(eval(wcall, wenv)))

      list(.adrftools_data = fit_data,
           .adrftools_weightit = weightit_boot) |>
        list2env(parent = genv) |>
        eval(expr = model_call)
    })
  }

  v <- try(stats::weights(x), silent = TRUE)

  if (null_or_error(v) || !is.numeric(v)) {
    v <- 1
  }
  else if (!is_null(x[["na.action"]]) &&
           length(v) == fnrow(fit_data) + length(x[["na.action"]])) {
    # With na.action = na.exclude, weights() pads the dropped rows with NAs
    v <- v[-x[["na.action"]]]
  }

  # Weights are passed as a column of the data so that model.frame() finds them
  model_call$weights <- quote(.adrftools_w)

  function(w) {
    fit_data$.adrftools_w <- v * w

    list(.adrftools_data = fit_data) |>
      list2env(parent = genv) |>
      eval(expr = model_call)
  }
}

# Design-based variance of estimates from a svyglm model, given the matrix `S`
# of their influence values (one row per unit, scaled so that
# crossprod(S / nrow(S)) is the i.i.d. variance). Following Binder (1983), the
# estimates are linearized as totals of the influence values; dividing by the
# sampling weights and totaling them with svytotal() lets survey apply the
# variance estimator appropriate to the design (clusters, strata,
# finite-population corrections, calibration, two-phase sampling, or replicate
# weights). The target weights are always the design weights (see
# process_wts_list()), which is what makes the calibration adjustment
# svytotal() applies to the whole influence value correct.
.svyglm_vcov <- function(x, S) {
  des <- x[["survey.design"]]

  w_d <- {
    if (inherits(des, "svyrep.design")) stats::weights(des, "sampling")
    else stats::weights(des)
  }

  n <- fnrow(S)
  na_act <- x[["na.action"]]

  # Designs that cannot drop rows (e.g., calibrated and PPS designs) keep the
  # units the fit dropped for missing values, with weights of 0
  keep <- {
    if (length(w_d) == n) seq_len(n)
    else if (!is_null(na_act) && length(w_d) == n + length(na_act) &&
             all(w_d[na_act] == 0)) -na_act
  }

  if (is_null(keep)) {
    arg::err("the survey design and the fitted model have different numbers of observations")
  }

  z <- matrix(0, nrow = length(w_d), ncol = fncol(S))
  z[keep, ] <- S / n / w_d[keep]

  # Units with weights of 0 (e.g., outside a subset() of the design) have
  # influence values of 0
  z[w_d == 0, ] <- 0

  # The matrix is passed directly rather than added to the design's variables,
  # which two-phase designs do not read and which may already have a column of
  # the same name
  v <- survey::svytotal(z, design = des) |>
    stats::vcov()

  matrix(v, nrow = fncol(z), ncol = fncol(z))
}

# Column order that turns estimates stored value-major (every subgroup at the
# first treatment value, then every subgroup at the second, ...) into the
# subgroup-major order used for `.est` throughout (every treatment value in the
# first subgroup, then the second, ...)
.value_major_to_by_major <- function(n_values, n_by = 1L) {
  matrix(seq_len(n_values * n_by), nrow = n_by) |>
    t() |>
    as.vector()
}

add_est_labels <- function(res, .contrast, .by_grid, .values = NULL, est_name = NULL) {

  add_values <- !is_null(.values)

  if (!add_values) {
    .values <- NA
  }

  # Apply by_grid and contrast labels if necessary
  if (!is_null(.contrast)) {
    contrast <- factor(rep(.contrast, each = length(.values)),
                       levels = .contrast)

    if (add_values) {
      res[[est_name]] <- rep.int(.values, length(.contrast))
    }

    res <- cbind(contrast = contrast, res)
  }
  else if (!is_null(.by_grid)) {
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

# Prepares the correlation matrix `m` of a set of estimates for
# mvtnorm::qmvt() and pmvt(). A matrix of full rank is returned unchanged. The
# correlation matrix of estimates along an effect curve often has low rank
# (e.g., rank 3 for the AMEF of a cubic polynomial), and qmvt() and pmvt() fail
# on such a matrix or on one nudged only slightly away from singularity. In
# that case, the matrix is rebuilt from its eigenvectors with eigenvalues
# above `rel_tol` times the largest, rescaled to have a unit diagonal, and
# combined with the identity matrix with weight `ridge`. With this
# regularization, critical values agree with those simulated from the
# low-rank factor to within Monte Carlo error.
.to_psd <- function(m, rel_tol = 1e-8, ridge = 1e-6) {
  es <- eigen(m, symmetric = TRUE)
  esv <- es$values

  keep <- esv > rel_tol * max(esv)

  if (all(keep)) {
    return(m)
  }

  L <- es$vectors[, keep, drop = FALSE] %*% diag(sqrt(esv[keep]), sum(keep))

  out <- (1 - ridge) * stats::cov2cor(tcrossprod(L)) + ridge * diag(fnrow(m))

  dimnames(out) <- dimnames(m)

  out
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

  if (inherits(x, c("bart", "stan4bartFit", "bartisan_fit"))) {
    return(TRUE)
  }

  if (is_null(data)) {
    data <- insight::get_data(x, verbose = FALSE)
  }

  p <- try(marginaleffects::get_predict(x, newdata = ss(data, seq_len(min(fnrow(data), 2L)))),
           silent = TRUE)

  !is_null(.attr(p, "posterior_draws"))
}

make_rmvt <- function(mu, Sigma, df = Inf, tol = 1e-7) {
  #Returns p x n matrix (easier to get col maxes)
  p <- length(mu)

  if (!allv(dim(Sigma), p)) {
    arg::err("incompatible arguments")
  }

  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values

  if (any(ev < -tol * abs(ev[1L]))) {
    arg::err("{.arg Sigma} must be positive definite")
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
# P(Q >= stat) for Q = sum_k ev[k] * Z_k^2, with Z_k independent standard
# normal, by simulation; with finite `df`, Q is divided by an independent
# chi^2_df / df, as for the quadratic form of a multivariate t vector. Each
# term is drawn directly, which needs only the eigenvalues, and the draws are
# made in chunks of `chunk` so memory does not grow with `n`.
.sim_p_value <- function(stat, ev, n = 1e6, df = Inf, chunk = 1e5) {
  ev <- ev[ev > 0]

  sizes <- c(rep.int(chunk, n %/% chunk), n %% chunk)
  sizes <- sizes[sizes > 0]

  hits <- vapply(sizes, function(m) {
    q <- numeric(m)

    for (k in seq_along(ev)) {
      q <- q + ev[k] * stats::rnorm(m)^2
    }

    if (is.finite(df)) {
      q <- q / (stats::rchisq(m, df) / df)
    }

    sum(q >= stat)
  }, numeric(1L))

  sum(hits) / n
}

#' Get kernel weights for interpolation
#' @noRd
#' @param x the predictor values to impute
#' @param v the values on the existing grid
#' @param kernel the name of kernel; `"gaussian"`, `"epanechnikov"`, or `"triangular"`
#' @param constant when `bw` is `NULL` and `v` is on a regularly spaced grid, `constant` times the spacing is used as the bandwidth.
#' @param bw the bandwidth (function or value) to use. If `NULL` and the grid of `v` is not regularly spaced, `bw.nrd()` is used.
#' @param ... ignored.
get_kernel_w <- function(x, v = x, kernel = "gaussian", constant = .5, bw = NULL, ...) {

  n <- length(v)
  tol <- .grid_tol(v)

  if (length(x) == n && all(abs(x - v) <= tol)) {
    return(diag(1, n, n))
  }

  if (all(x %in% v)) {
    w <- outer(x, v, function(a, b) abs(a - b) <= tol)

    return(w / rowSums(w))
  }

  bw <- process_bw(bw, v, constant)

  kernel <- arg::match_arg(kernel, c("gaussian", "epanechnikov",
                                     "triangular"))

  # Kernels are evaluated on standardized distances, and the normalizing
  # constants are omitted because each row is normalized below, so the weights
  # do not depend on the scale of the treatment
  k <- switch(kernel,
              gaussian = function(a, b) {
                dnorm((a - b) / bw)
              },
              epanechnikov = function(a, b) {
                az <- abs(a - b) / (bw * sqrt(5))
                o <- alloc(0, length(az))
                in_window <- az < 1
                o[in_window] <- .75 * (1 - az[in_window]^2)
                o
              },
              triangular = function(a, b) {
                az <- abs(a - b) / (bw * sqrt(6))
                o <- alloc(0, length(az))
                in_window <- az < 1
                o[in_window] <- 1 - az[in_window]
                o
              })

  w <- outer(x, v, k)

  rs <- rowSums(w)

  zeros <- rs <= 0

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
#' @param p the degree of the local polynomial. Set to 0 to use Nadaraya-Watson interpolation.
#' @param ... args passed to `get_kernel_w()` when `w` is `NULL`
get_locpoly_w <- function(x, v = x, w = NULL, p = 3, ...) {
  arg::arg_count(p)

  if (is_null(w)) {
    w <- get_kernel_w(x = x, v = v, ...)
  }

  arg::arg_equal(length(x), fnrow(w))
  arg::arg_equal(length(v), fncol(w))

  if (p == 0) {
    return(w)
  }

  n <- length(v)

  # A polynomial basis requires fewer degrees than unique points
  p <- min(p, length(funique(c(v, x))) - 1L)

  # Each local polynomial is fit in coordinates centered at the point being
  # interpolated and scaled by the width of the window, so its basis is well
  # conditioned however narrow the window is relative to the grid, and the
  # interpolated value is the intercept. A basis orthogonal over the whole grid
  # is nearly collinear within a narrow window (e.g., with 101 grid points),
  # which made the decomposition below rank deficient.
  lw <- matrix(0, nrow = length(x), ncol = n)

  for (i in seq_along(x)) {
    wi <- w[i, ]
    pos <- which(wi > 1e-10)

    if (is_null(pos)) {
      arg::err("all kernel weights estimated as 0, indicating severe extrapolation")
    }

    wi <- wi[pos]
    wi <- wi / sum(wi)

    # A local polynomial of degree p needs at least p + 1 grid points
    p_i <- min(p, length(pos) - 1L)

    z <- v[pos] - x[i]
    z <- z / max(abs(z), .Machine$double.eps)

    Z <- outer(z, 0:p_i, `^`)

    # Weighted least squares through the QR decomposition of sqrt(w) Z. The
    # kernel weights can span many orders of magnitude, so solving the normal
    # equations (Z'WZ)^-1 Z'W loses precision, which numerical derivatives of
    # the interpolated curve amplify
    sw <- sqrt(wi)
    qr_i <- qr(Z * sw)

    lw[i, pos] <- qr.coef(qr_i, diag(sw, nrow = length(pos)))[1L, ]
  }

  lw
}

# Pool multiply imputed estimates using Rubin's rules
pool_est <- function(est.list, vcov.list = NULL) {
  est.mat <- do.call("rbind", est.list)

  out <- list(est = colMeans(est.mat))

  if (!is_null(vcov.list)) {
    m <- length(vcov.list)

    v_w <- Reduce("+", vcov.list) / m

    v_b <- cov(est.mat)

    out$vcov <- v_w + (1 + 1 / m) * v_b
    out$vcov_w <- v_w
    out$vcov_b <- v_b
  }

  out
}

# Denominator degrees of freedom for joint inference on estimates pooled from
# m imputed datasets: Reiter (2007) when the complete-data degrees of freedom
# `df_com` are finite, and Li, Raghunathan, and Rubin (1991) otherwise. These
# assume the fraction of missing information is the same for every estimate,
# so the same degrees of freedom apply to any linear combination of them.
.mi_df <- function(v_w, v_b, m, df_com = Inf) {
  # The within-imputation covariance of points along a curve is rank deficient
  # when there are more points than model coefficients, so the average relative
  # increase in variance is computed in its span
  e <- eigen(v_w, symmetric = TRUE)
  keep <- e$values > sqrt(.Machine$double.eps) * max(e$values, 0)

  if (!any(keep)) {
    return(df_com)
  }

  Q <- e$vectors[, keep, drop = FALSE]
  k <- sum(keep)

  r <- (1 + 1 / m) * sum(diag(crossprod(Q, v_b %*% Q)) / e$values[keep]) / k

  if (r <= sqrt(.Machine$double.eps)) {
    return(df_com)
  }

  t <- k * (m - 1)

  df_li <- {
    if (t > 4) 4 + (t - 4) * (1 + (1 - 2 / t) / r)^2
    else t * (1 + 1 / k) * (1 + 1 / r)^2 / 2
  }

  if (!is.finite(df_com)) {
    return(df_li)
  }

  v_star <- (df_com + 1) / (df_com + 3) * df_com

  if (t > 4) {
    a <- r * t / (t - 2)
    c0 <- 1 / (t - 4)
    c1 <- v_star - 2 * (1 + a)
    c2 <- v_star - 4 * (1 + a)

    z <- 1 / c2 +
      c0 * (a^2 * c1 / ((1 + a)^2 * c2)) +
      c0 * (8 * a^2 * c1 / ((1 + a) * c2^2) + 4 * a^2 / ((1 + a) * c2)) +
      c0 * (4 * a^2 / (c2 * c1) + 16 * a^2 * c1 / c2^3) +
      c0 * (8 * a^2 / c2^2)

    df_reiter <- 4 + 1 / z

    if (is.finite(df_reiter) && df_reiter > 0) {
      return(df_reiter)
    }
  }

  # With few estimates or imputations, combine the large-sample degrees of
  # freedom with those of the observed data, as in Barnard and Rubin (1999)
  1 / (1 / df_li + (1 + r) / v_star)
}

# Get quadrature weights for trapezoidal Riemann sum along given grid
# Kernel density estimates of the treatment in each subgroup, for weighting
# the grid by the distribution of the treatment. `A_list`, `w_list`, and
# `by_list` hold, for each imputed dataset, the treatment values, target
# weights, and subgroup codes of the units used; with multiply imputed data,
# each dataset's units get weight 1 / m. Each density uses stats::density()
# with the bandwidth from Scott's rule (bw.nrd()) for the subgroup's treatment values, and all are
# evaluated on a common fine grid spanning the treatment values, so a density
# can be interpolated at any grid point (including points added later, such as
# a reference value). Returns the fine grid `x`, a matrix `y` with one column
# per subgroup, and each subgroup's total target weight `tot`, used to pool
# subgroups for contrasts.
.get_treat_density <- function(A_list, w_list, by_list, n_by, n_fine = 512L) {
  m <- length(A_list)

  A <- unlist(A_list)
  w <- unlist(lapply(w_list, function(wi) wi / m))
  b <- unlist(by_list)

  A_split <- gsplit(A, b)
  w_split <- gsplit(w, b)
  tot <- fsum(w, b, use.g.names = FALSE)

  bws <- vapply(A_split, function(Ak) {
    if (fnunique(Ak) > 1L) bw.nrd(Ak) else 1
  }, numeric(1L))

  x <- seq(fmin(A) - 3 * max(bws), fmax(A) + 3 * max(bws), length.out = n_fine)

  y <- vapply(seq_len(n_by), function(k) {
    stats::density(A_split[[k]], weights = w_split[[k]] / tot[k], bw = bws[k],
                   from = x[1L], to = x[n_fine], n = n_fine)$y
  }, numeric(n_fine))

  list(x = x,
       y = matrix(y, nrow = n_fine),
       tot = tot)
}

# Weights of the grid points for omnibus tests and projections, one vector
# per block of estimates (subgroup or contrast) in `blocks`: the trapezoidal
# weights of the grid ("uniform"), or those times the density of the
# treatment at each grid point ("empirical"), normalized to sum to 1
get_grid_w <- function(values, integrate = "uniform", treat_density = NULL,
                       blocks = 1L) {
  w0 <- get_trapezoidal_w(values)

  if (integrate == "uniform") {
    return(rep.int(list(w0), length(blocks)))
  }

  lapply(blocks, function(k) {
    f <- stats::approx(treat_density[["x"]], treat_density[["y"]][, k],
                       xout = values, rule = 2L)$y
    w <- w0 * pmax(f, 0)
    w / sum(w)
  })
}

# Processes the `integrate` argument of summary.effect_curve() and
# curve_projection()
process_integrate <- function(integrate, x) {
  integrate <- arg::match_arg(integrate, c("uniform", "empirical"))

  if (integrate == "empirical") {
    if (!is_null(.attr(x, ".shift"))) {
      arg::err("{.code integrate = \"empirical\"} cannot be used with a shift curve, because the shift has no distribution to weight by")
    }

    if (is_null(.attr(x, ".treat_density"))) {
      arg::err("{.code integrate = \"empirical\"} requires the distribution of the treatment, which is stored by {.fun adrf} in versions of {.pkg adrftools} after 0.1.0; rerun {.fun adrf} to use it")
    }
  }

  integrate
}

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
  nr <- fnrow(x)

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

  if (is_null(head_ind)) {
    tmp <- utils::capture.output({
      ss(x, 1L) |>
        round_df_char(digits = digits, na_vals = ".") |>
        print.data.frame(row.names = rownames,
                         ...)
    })

    out <- character(0L)
  }
  else {
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

    if (!is_null(tail_ind)) {
      out <- c(out, tmp[-seq_along(c(1L, head_ind))])
    }
  }

  nc <- max(nchar(out))

  if (!is_null(to_it)) {
    out[to_it] <- .it(out[to_it])
  }

  if (!is_null(main)) {
    #Add title
    out <- c(center_just(main, wrt = space(nc)),
             txtbar(nc),
             out)
  }

  #Add bar
  if (bar) {
    out <- c(out, txtbar(nc))
  }

  cli::cat_line(out)

  invisible(out)
}

.format_p_value <- function(p, digits = 4L) {
  p_text <- alloc(NA_character_, length(p))

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
