process_by <- function(by, data) {
  if (is_null(by)) {
    return(list(by_id = qG(alloc(1L, nrow(data)))))
  }

  if (is.character(by)) {
    by <- reformulate(by)
  }
  else if (!inherits(by, "formula")) {
    .err("`by` must be a one-sided formula or character vector")
  }

  by_mf <- model.frame(update(by, NULL ~ .), data = data)
  attr(by_mf, "terms") <- NULL
  by_name <- names(by_mf)

  by_grid <- funique(by_mf)
  by_grid <- ss(by_grid, do.call("order", as.list(by_grid)))

  by_grid$.by_id <- seq_row(by_grid)

  by_mf <- by_mf |>
    cbind(.merge_id = seq_row(by_mf)) |>
    merge(by_grid,
          by = seq_along(by_name),
          all.x = TRUE, all.y = FALSE,
          sort = FALSE)

  by_mf <- ss(by_mf, order(by_mf$.merge_id))

  list(by_grid = by_grid[-ncol(by_grid)],
       by_id = qG(by_mf[[length(by_mf)]]))
}

process_subset <- function(index.sub, data, env = parent.frame(2L)) {

  if (is_null(index.sub)) {
    return(seq_row(data))
  }

  subset <- eval(index.sub, data, env)

  if (!chk::vld_atomic(subset)) {
    .err("`subset` must evaluate to an atomic vector")
  }

  if (is.logical(subset) && length(subset) != nrow(data)) {
    .err("when `subset` is logical, it must have the same length as the original dataset")
  }

  if (is.logical(subset)) {
    if (length(subset) != nrow(data)) {
      .err("when `subset` is logical, it must have the same length as the original dataset")
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

process_values <- function(values = NULL, n, treat_var, q = .95, strict = FALSE) {
  if (is_null(values)) {
    chk::chk_count(n)
    chk::chk_gte(n, 2)

    if (q == 1) {
      range_t <- .range(treat_var)
    }
    else {
      range_t <- .quantile(treat_var, probs = c(1 - q, 1 + q) / 2)
    }

    return(seq(range_t[1L], range_t[2L], length.out = n))
  }

  chk::chk_numeric(values)
  values <- sort(values)

  range_t <- .range(treat_var)

  if (strict) {
    if (any(values < range_t[1L]) || any(values > range_t[2L])) {
      r <- format(range_t, digits = 4L, drop0trailing = TRUE)
      .err(sprintf("`values` cannot be outside the range of observed treatment values (%s to %s)",
                   r[1L],
                   r[2L]))
    }
  }
  else {
    range_diff <- range_t[2L] - range_t[1L]

    if (any(values < range_t[1L] - .1 * range_diff) ||
        any(values > range_t[2L] + .1 * range_diff)) {
      .wrn("some `values` are outside the range of observed treatment values")
    }
  }

  values
}

process_reference <- function(reference = NULL, treat_var, strict = FALSE) {
  if (is_null(reference)) {
    return(NULL)
  }

  chk::chk_number(reference)

  range_t <- .range(treat_var)

  if (strict) {
    if (any(reference < range_t[1L]) || any(reference > range_t[2L])) {
      r <- format(range_t, digits = 4L, drop0trailing = TRUE)
      .err(sprintf("`reference` cannot be outside the range of observed treatment values (%s to %s)",
                   r[1L],
                   r[2L]))
    }
  }
  else {
    range_diff <- range_t[2L] - range_t[1L]

    if (any(reference < range_t[1L] - .1 * range_diff) ||
        any(reference > range_t[2L] + .1 * range_diff)) {
      .wrn("`reference` is outside the range of observed treatment values")
    }
  }

  reference
}

process_eps <- function(eps, treat_var) {
  chk::chk_number(eps)
  chk::chk_gt(eps, 0)

  eps * sd(treat_var)
}

process_treat <- function(treat, data) {
  chk::chk_string(treat)
  chk::chk_subset(treat, names(data))

  treat_var <- data[[treat]]

  if (!chk::vld_numeric(treat_var)) {
    .err("`treat` must be the name of a numeric variable in the original dataset")
  }

  treat_var
}

process_wts <- function(wts, model) {
  if (is_null(wts)) {
    if (endsWith(class(model)[1L], "_weightit")) {
      wts <- "(s.weights)"
    }
    else if (inherits(model, "svyglm")) {
      wts <- TRUE
    }
    else {
      wts <- FALSE
    }
  }

  n <- insight::n_obs(model)

  if (isFALSE(wts)) {
    return(alloc(1.0, n))
  }

  if (isTRUE(wts)) {
    return(process_wts(insight::get_weights(model), model))
  }

  if (is.numeric(wts)) {
    if (length(wts) != n) {
      .err("when supplied as a numeric vector, `wts` must have as many values as the number of observations in the original model")
    }

    return(as.vector(wts))
  }

  if (!chk::vld_string(wts)) {
    .err("`wts` must be `NULL`, `TRUE`, `FALSE`, a numeric vector, or a string containing the name of the variable containing weights in the original dataset")
  }

  if (wts %in% c("(weights)", "(s.weights)") && is_not_null(model[["model"]]) &&
      utils::hasName(model[["model"]], wts)) {
    data <- model[["model"]]
  }
  else {
    data <- insight::get_data(model, verbose = FALSE)
  }

  if (!utils::hasName(data, wts)) {
    .err("when supplied as a string, `wts` must be the name of the variable in the original dataset")
  }

  wts <- data[[wts]]

  if (!is.numeric(wts)) {
    .err("when supplied as a string, `wts` must be the name of a numeric variable")
  }

  wts
}

process_transform <- function(transform = NULL, family = NULL) {
  inv_transform <- NULL
  d_transform <- NULL

  if (isTRUE(transform) && is_not_null(family) && !identical(family$link, "identity")) {
    transform <- family$linkfun
    inv_transform <- family$linkinv
    d_transform <- switch(family$link,
                          "log" = function(x) {1 / x},
                          "logit" = function(x) {1 / ((1 - x) * x)},
                          "probit" = function(x) {1 / dnorm(qnorm(x))},
                          "cauchit" = function(x) {1 / dcauchy(qcauchy(x))},
                          "sqrt" = function(x) {.5 / sqrt(x)},
                          "cloglog" = function(x) {-1 / ((1 - x) * log(1 - x))},
                          "inverse" = function(x) {-1 / sqrt(x)},
                          "1/mu^2" = function(x) {-2 / x^3},
                          NULL)
  }
  else if (is_null(transform) || isFALSE(transform) || is_null(family) || identical(family$link, "identity")) {
    transform <- identity
    inv_transform <- identity
    d_transform <- function(x) {
      alloc(1.0, length(x))
    }
  }
  else if (is.list(transform) && is_not_null(names(transform)) &&
           utils::hasName(transform, "transform") &&
           all(names(transform) %in% c("transform", "inv_transform", "d_transform"))) {

    if (utils::hasName(transform, "inv_transform")) {
      inv_transform <- transform[["inv_transform"]]
    }

    if (utils::hasName(transform, "d_transform")) {
      d_transform <- transform[["d_transform"]]
    }

    transform <- transform[["transform"]]
  }

  if (!is.function(transform)) {
    .err("`transform` must be `TRUE`, `FALSE`, or a function")
  }

  if (is_null(inv_transform)) {
    if (identical(transform, exp)) {
      inv_transform <- log
      d_transform <- exp
    }
    else if (identical(transform, log)) {
      inv_transform <- exp
      d_transform <- function(x) {1 / x}
    }
    else if (identical(transform, sqrt)) {
      inv_transform <- function(x) {x^2}
      d_transform <- function(x) {.5 / sqrt(x)}
    }
    else if (identical(transform, qlogis)) {
      inv_transform <- plogis
      d_transform <- function(x) {
        1 / ((1 - x) * x)
      }
    }
    else if (identical(transform, qnorm)) {
      inv_transform <- pnorm
      d_transform <- function(x) {
        1 / dnorm(qnorm(x))
      }
    }
  }

  if (is_null(d_transform)) {
    d_transform <- function(x) {
      eps <- 1e-6
      (transform(x + eps / 2) - transform(x - eps / 2)) / eps
    }
  }

  list(transform = transform,
       inv_transform = inv_transform,
       d_transform = d_transform)
}

process_vcov <- function(vcov, model) {
  if (is_null(vcov) || isFALSE(vcov)) {
    return("none")
  }

  builtin_vcovs <- c("none", "unconditional", "conditional", "boot", "fwb")

  if (chk::vld_string(vcov)) {
    p <- pmatch(vcov, builtin_vcovs, nomatch = 0L)

    if (p != 0L) {
      if (builtin_vcovs[p] == "unconditional") {
        if (null_or_error(try(sandwich::bread(model), silent = TRUE)) ||
            null_or_error(try(sandwich::estfun(model), silent = TRUE))) {
          .wrn('`vcov` cannot be "unconditional" with this type of model. Using "conditional"')

          return(Recall(vcov = "conditional", model = model))
        }

        return("unconditional")
      }

      if (builtin_vcovs[p] == "conditional") {
        vcov <- marginaleffects::get_vcov
      }
      else if (builtin_vcovs[p] %in% c("boot", "fwb")) {
        rlang::check_installed("fwb (>= 0.5.0)")
        return(builtin_vcovs[p])
      }
      else {
        return(builtin_vcovs[p])
      }
    }
  }

  vcov_try <- try(marginaleffects::get_vcov(model, vcov = vcov),
                  silent = TRUE)

  if (null_or_error(vcov_try)) {
    .err(sprintf('`vcov` must be %s or one of the allowed arguments to `marginaleffects::get_vcov()`',
                 word_list(builtin_vcovs, and.or = "or", quotes = TRUE)))
  }

  vcov_try
}

process_vcov_mi <- function(vcov, models) {
  if (is_null(vcov) || isFALSE(vcov)) {
    return("none")
  }

  builtin_vcovs <- c("none", "unconditional", "conditional", "boot", "fwb")

  if (chk::vld_string(vcov)) {
    p <- pmatch(vcov, builtin_vcovs, nomatch = 0L)

    if (p != 0L) {
      if (builtin_vcovs[p] %in% c("boot", "fwb")) {
        .err(sprintf("`vcov` cannot be %s with multiply imputed data",
                     add_quotes(builtin_vcovs[p])))
      }

      if (builtin_vcovs[p] == "unconditional") {
        for (model in models) {
          if (null_or_error(try(sandwich::bread(model), silent = TRUE)) ||
              null_or_error(try(sandwich::estfun(model), silent = TRUE))) {
            .wrn('`vcov` cannot be "unconditional" with this type of model. Using "conditional"')

            return(Recall(vcov = "conditional", models = models))
          }
        }
        vcov <- "unconditional"
      }
      else if (builtin_vcovs[p] == "conditional") {
        vcov <- marginaleffects::get_vcov
      }
      else {
        vcov <- builtin_vcovs[p]
      }

      return(rep_with(list(vcov), models))
    }
  }

  if (!is.list(vcov)) {
    vcov <- rep_with(list(vcov), models)
  }
  else if (length(vcov) != length(models)) {
    .err("`vcov` must have length equal to the number of imputed datasets")
  }

  for (.i in seq_along(models)) {
    vcov[[.i]] <- try(marginaleffects::get_vcov(models[[.i]], vcov = vcov[[.i]]),
                      silent = TRUE)

    if (null_or_error(vcov[[.i]])) {
      .err("`vcov` must be one of %s or one of the allowed arguments to `marginaleffects::get_vcov()` or a list thereof",
           word_list(setdiff(builtin_vcovs, c("bs", "fwb")),
                     and.or = "or", quotes = TRUE))
    }
  }

  vcov
}

process_call <- function(definition = sys.function(sys.parent()),
                         call = sys.call(sys.parent()),
                         expand.dots = TRUE) {

  mcall <- match.call(definition = definition, call = call,
                      expand.dots = expand.dots,
                      envir = parent.frame(3L))

  generic <- get0(".Generic", parent.frame(1L), inherits = FALSE)

  if (is_not_null(generic)) {
    mcall[[1L]] <- str2lang(generic)
  }

  mcall
}

process_data <- function(model) {
  if (is.environment(model[["data"]])) {
    vars <- insight::find_predictors(model, effects = "fixed", component = "all",
                                     flatten = TRUE, verbose = FALSE)

    dat <- model[["data"]] |> as.list() |> list2DF()
    dat <- dat[, intersect(vars, colnames(dat)), drop = FALSE]

    return(dat)
  }

  insight::get_predictors(model, verbose = FALSE)
}

chk_proj <- function(x, proj = NULL) {
  if (is_not_null(proj)) {
    .chk_is(proj, "curve_proj")

    if ((inherits(x, c("adrf_est", "amef_est")) && !isTRUE(all.equal(x, .attr(proj, "obj")))) ||
        inherits(x, c("summary.adrf_est", "summary.amef_est")) && !isTRUE(all.equal(.attr(x, "obj"), .attr(proj, "obj")))) {
      .err(sprintf("`proj` must be the output of `curve_proj()` applied to the original %s object",
                   class(x)[1L]))
    }
  }
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
  dm <- tcrossprod(es$vectors %*% diag(sqrt(tau), d))

  m + dm
}

#Make a function that yields its input when fun is applied to its output
make_inverse <- function(fun, original) {
  function(x) {

    # length 1 x doesn't always play nice, so repeat it
    d <- length(x) == 1L
    if (d) {
      x <- x[c(1L, 1L)]
    }

    # Get starter closest to transformed x
    start <- original[vapply(x, function(y) which.min(abs(fun(original) - y)),
                             integer(1L))]

    out <- optim(function(y) mean((fun(y) - x)^2),
                 par = start,
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

  if (!all(dim(Sigma) == p)) {
    .err("incompatible arguments")
  }

  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values

  if (any(ev < -tol * abs(ev[1L]))) {
    .err("`Sigma` is not positive definite")
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

.sim_p_value <- function(stat, v, n = 1e6, max_size = 1e7, df = Inf) {
  rmvt <- make_rmvt(alloc(0.0, ncol(v)), v, df = df)

  d <- nrow(v)

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

get_df <- function(fit) {
  if (inherits(fit, "mira")) {
    df <- vapply(fit[["analyses"]], get_df, numeric(1L)) |>
      mean()

    return(df)
  }

  if (!insight::is_model_supported(fit)) {
    return(Inf)
  }

  statistic <- insight::find_statistic(fit)

  if (identical(statistic, "chi-squared statistic")) {
    return(Inf)
  }

  insight::get_df(fit, type = "wald", statistic = statistic)
}

get_kernel_w <- function(x, v = x, kernel = "gaussian", constant = 1, bw = NULL, ...) {

  n <- length(v)
  if (length(x) == n && all(check_if_zero(x - v))) {
    return(diag(1, n, n))
  }

  if (all(x %in% v)) {
    w <- outer(x, v, "==")
    return(w / rowSums(w))
  }

  if (is_not_null(bw)) {
    chk::chk_number(bw)
    chk::chk_gt(bw, 0)
  }
  else if (all_the_same(diff1(v))) {
    bw <- constant * mean(diff1(v))
  }
  else {
    bw <- bw.nrd(v)
  }

  kernel <- match_arg(kernel, c("gaussian", "epanechnikov"))

  k <- switch(kernel,
              gaussian = function(a, b) dnorm(a - b, sd = bw),
              epanechnikov = function(a, b) {
                bw_a <- bw * sqrt(5)
                ax <- abs(a - b)
                o <- alloc(0, length(ax))
                in_window <- ax < bw_a
                o[in_window] <- 3/4 * (1 - (ax[in_window] / bw_a)^2) / bw_a
                o
              })

  w <- outer(x, v, k)

  rs <- rowSums(w)

  zeros <- check_if_zero(rs)

  if (any(zeros)) {
    w[zeros,] <- 0
  }

  if (!all(zeros)) {
    w[!zeros,] <- w[!zeros,] / rs[!zeros]
  }

  w
}

get_locpoly_w <- function(x, v = x, w = NULL, p = 3, ...) {
  chk::chk_count(p)

  if (is_null(w)) {
    w <- get_kernel_w(x = x, v = v, ...)
  }

  chk::chk_equal(length(x), nrow(w))
  chk::chk_equal(length(v), ncol(w))

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
    wi <- w[i,]
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

get_trapezoidal_w <- function(grid) {
  val_diff <- diff1(grid)

  (c(0, val_diff) + c(val_diff, 0)) / 2
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

  out
}
