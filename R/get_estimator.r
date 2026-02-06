.get_model_type <- function(x, type = "response") {
  if (.is_bayesian(x)) {
    return("default")
  }

  if (class(x)[1L] %in% "ordinal_weightit" && type %in% c("link", "lp")) {
    return("lm")
  }

  if (class(x)[1L] %in% c("lm", "lm_robust")) {
    return("lm")
  }

  if (class(x)[1L] %in% c("glm", "glm_weightit", "svyglm", "negbin", "fixest")) {
    .family <- insight::get_family(x)

    if (identical(.family[["link"]], "identity")) {
      return("lm")
    }

    if (type == "response") {
      return("glm")
    }
  }

  "default"
}

.get_estimator <- function(.x, ...) {
  UseMethod(".get_estimator")
}

#' @exportS3Method NULL
.get_estimator.lm <- function(.x, ...) {
  function(mm, beta0 = NULL, model = NULL, ...) {
    if (is_null(beta0)) {
      beta0 <- marginaleffects::get_coef(model)
    }

    tcrossprod(na_rm(beta0), mm) |>
      drop()
  }
}

#' @exportS3Method NULL
.get_estimator.glm <- function(.x, ...) {
  .linkfun <- process_family(.x)$linkinv

  function(mm, beta0 = NULL, model = NULL, ...) {
    if (is_null(beta0)) {
      beta0 <- marginaleffects::get_coef(model)
    }

    tcrossprod(na_rm(beta0), mm) |>
      drop() |>
      .linkfun()
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
      if (is_not_null(beta0)) {
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
