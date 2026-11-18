# Determines how predictions are computed: from the model matrix and the
# coefficients with an identity link ("lm") or an inverse link ("glm"), which is
# fast and gives analytic gradients, or with marginaleffects::get_predict()
# ("default"). Rather than relying on the class of the model, the model-matrix
# predictions are compared with get_predict() on a sample of rows of the
# counterfactual data, so that models whose predictions involve more than the
# model matrix, the coefficients, and any offsets (auxiliary parameters, fixed
# effects, thresholds, random effects) take the general path.
.get_model_type <- function(x, type = "response", data_grid, model_data = NULL,
                            is_bayes = FALSE, ...) {
  if (is_bayes) {
    return("default")
  }

  .family <- process_family(x, type)

  rows <- seq(1, fnrow(data_grid), length.out = min(fnrow(data_grid), 20L)) |>
    round() |>
    funique()

  nd <- ss(data_grid, rows)

  fast_ok <- try({
    beta <- marginaleffects::get_coef(x)
    aliased <- is.na(beta)

    X <- get_tt(x, model_data) |>
      model.matrix(data = nd, xlev = .model_xlev(x))

    if (fncol(X) != length(beta) ||
        (!is_null(names(beta)) && !identical(colnames(X), names(beta)))) {
      FALSE
    }
    else {
      off <- .grid_offset(x, nd, model_data)

      fast <- (drop(ss(X, j = !aliased) %*% beta[!aliased]) + off) |>
        .family$linkinv()

      p <- .get_predict(x, newdata = nd, type = type, ...)

      length(p) == length(fast) &&
        isTRUE(all.equal(unname(p), unname(fast), tolerance = 1e-8,
                         check.attributes = FALSE))
    }
  }, silent = TRUE)

  if (!isTRUE(fast_ok)) {
    return("default")
  }

  if (identical(.family[["link"]], "identity")) "lm"
  else "glm"
}

# The offset for each row of `data`: offset() terms in the formula and an
# `offset` argument in the model call, evaluated on `data`. Returns 0 when the
# model has no offset.
.grid_offset <- function(x, data, model_data = NULL) {
  off <- 0

  tt <- try(get_tt(x, model_data), silent = TRUE)

  # The model frame of the whole grid is built only when the formula has an
  # offset() term
  if (!null_or_error(tt) && !is_null(.attr(tt, "offset"))) {
    mf <- try(stats::model.frame(tt, data = data, na.action = stats::na.pass),
              silent = TRUE)

    if (!null_or_error(mf)) {
      o <- stats::model.offset(mf)

      if (!is_null(o)) {
        off <- off + o
      }
    }
  }

  call_off <- try(stats::getCall(x)[["offset"]], silent = TRUE)

  if (!null_or_error(call_off) && !is_null(call_off)) {
    env <- try(environment(stats::formula(x)), silent = TRUE)

    if (!is.environment(env)) {
      env <- globalenv()
    }

    o <- eval(call_off, data, env)

    if (!is.numeric(o) || length(o) != fnrow(data)) {
      stop("the offset could not be evaluated")
    }

    off <- off + o
  }

  off
}

# Which columns differ between two model matrices of the same dimensions; a
# missing value in either counts as a difference
.cols_differ <- function(a, b) {
  d <- a != b
  d[is.na(d)] <- TRUE

  colSums(d) > 0
}

.get_linkinv <- function(.x, ...) {
  UseMethod(".get_linkinv")
}

# The function that maps the linear predictor to the predictions
#' @exportS3Method NULL
.get_linkinv.lm <- function(.x, ...) {
  function(eta) {
    eta
  }
}

#' @exportS3Method NULL
.get_linkinv.glm <- function(.x, ...) {
  process_family(.x)$linkinv
}

#' @exportS3Method NULL
.get_linkinv.default <- function(.x, ...) {
  NULL
}


.get_grad_fun <- function(.x, ...) {
  UseMethod(".get_grad_fun")
}

#' @exportS3Method NULL
# The gradient of each prediction with respect to the coefficients is its row
# of the model matrix times the returned value (the derivative of the inverse
# link at the linear predictor); returning the multiplier rather than the
# product avoids a second copy of the model matrix
.get_grad_fun.lm <- function(.x, ...) {
  function(eta) {
    1
  }
}

#' @exportS3Method NULL
.get_grad_fun.glm <- function(.x, ...) {
  .mu.eta <- process_family(.x)$mu.eta

  function(eta) {
    .mu.eta(eta)
  }
}

#' @exportS3Method NULL
.get_grad_fun.default <- function(.x, ...) {
  NULL
}


.get_pred_fun <- function(.x, ...) {
  UseMethod(".get_pred_fun")
}

#' @exportS3Method NULL
.get_pred_fun.lm <- function(.x, type = NULL, ...) {
  NULL
}

#' @exportS3Method NULL
.get_pred_fun.glm <- function(.x, type = NULL, ...) {
  NULL
}

#' @exportS3Method NULL
.get_pred_fun.default <- function(.x, type = NULL, ...) {
  if (inherits(.x, "stan4bartFit")) {
    rlang::check_installed("stan4bart")

    function(data_grid, beta0 = NULL, model, ...) {
      p <- predict(object = model, newdata = data_grid,
                   type = "ev", ...)

      out <- rowMeans(p)
      attr(out, "posterior_draws") <- p

      out
    }
  }
  else {
    function(data_grid, beta0 = NULL, model, ...) {
      if (!is_null(beta0)) {
        model <- marginaleffects::set_coef(model, beta0)
      }

      .get_predict(model, newdata = data_grid, type = type, ...)
    }
  }
}

# Predictions from marginaleffects::get_predict(), with the posterior draws (if
# any) in the "posterior_draws" attribute. For svyglm models,
# stats::predict.glm() is used instead, because survey::predict.svyglm() (which
# get_predict() calls) rebuilds the model terms without their stored
# `predvars`, so terms whose basis depends on the data, such as ns(), poly(),
# or scale(), are recomputed from `newdata` and give wrong predictions on
# counterfactual data.
.get_predict <- function(model, newdata, type = "response", ...) {
  if (inherits(model, "svyglm") && anyv(c("response", "link"), type)) {
    out <- stats::predict.glm(model, newdata = newdata, type = type) |>
      as.numeric()

    return(out)
  }

  rlang::with_options(marginaleffects_posterior_center = "mean", {
    p <- marginaleffects::get_predict(model = model, newdata = newdata,
                                      type = type, ...)
  })

  out <- p[["estimate"]]
  attr(out, "posterior_draws") <- .attr(p, "posterior_draws")

  out
}
