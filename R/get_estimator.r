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
      model.matrix(data = nd)

    if (fncol(X) != length(beta) ||
        (!is_null(names(beta)) && !identical(colnames(X), names(beta)))) {
      FALSE
    }
    else {
      off <- .grid_offset(x, nd, model_data)

      fast <- (drop(ss(X, j = !aliased) %*% beta[!aliased]) + off) |>
        .family$linkinv()

      p <- marginaleffects::get_predict(x, newdata = nd, type = type, ...)[["estimate"]]

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

  mf <- try(stats::model.frame(get_tt(x, model_data), data = data,
                               na.action = stats::na.pass),
            silent = TRUE)

  if (!null_or_error(mf)) {
    o <- stats::model.offset(mf)

    if (!is_null(o)) {
      off <- off + o
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

.get_estimator <- function(.x, ...) {
  UseMethod(".get_estimator")
}

#' @exportS3Method NULL
.get_estimator.lm <- function(.x, ...) {
  function(mm, beta0 = NULL, model = NULL, offset = 0, ...) {
    if (is_null(beta0)) {
      beta0 <- marginaleffects::get_coef(model)
    }

    drop(tcrossprod(na_rm(beta0), mm)) + offset
  }
}

#' @exportS3Method NULL
.get_estimator.glm <- function(.x, ...) {
  .linkinv <- process_family(.x)$linkinv

  function(mm, beta0 = NULL, model = NULL, offset = 0, ...) {
    if (is_null(beta0)) {
      beta0 <- marginaleffects::get_coef(model)
    }

    (drop(tcrossprod(na_rm(beta0), mm)) + offset) |>
      .linkinv()
  }
}

#' @exportS3Method NULL
.get_estimator.default <- function(.x, ...) {
  NULL
}


.get_grad_fun <- function(.x, ...) {
  UseMethod(".get_grad_fun")
}

#' @exportS3Method NULL
.get_grad_fun.lm <- function(.x, ...) {
  function(mm, eta = NULL) {
    mm
  }
}

#' @exportS3Method NULL
.get_grad_fun.glm <- function(.x, ...) {
  .mu.eta <- process_family(.x)$mu.eta

  function(mm, eta = NULL) {
    .mu.eta(eta) * mm
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

      rlang::with_options(marginaleffects_posterior_center = "mean", {
        p <- marginaleffects::get_predict(model = model, newdata = data_grid,
                                          type = type, ...)
      })

      out <- p$estimate
      attr(out, "posterior_draws") <- .attr(p, "posterior_draws")

      out
    }
  }
}
