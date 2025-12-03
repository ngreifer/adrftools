#' Contrast point estimates along an effect curve
#'
#' @description
#' `point_contrast()` computes pairwise contrasts of estimates from an effect curve.
#'
#' @inheritParams summary.curve_est
#' @param object for `point_contrast()`, a [`curve_est`] object; the output of a an `effect_curve` object. For `summary()`, a `curve_est_contrast` object; the output of a call to `point_contrast()`.
#' @param null the null value for hypothesis tests. Default is 0. Set to `NA` to omit tests.
#'
#' @returns
#' `point_contrast()` returns an object of class `curve_est_contrast`, which is like a `curve_est` object but with its own `summary()` method.
#'
#' @details
#' `point_contrast()` computes all pairwise contrasts between effect curve estimates. Because pairwise contrasts are a linear operation over the original estimates, the delta method can be used to perform Wald inference for the contrasts. When `by` was specified in the original call to [adrf()] or the effect curve is a `contrast_curve` object resulting from [curve_contrast()], pairwise contrasts occur only within subgroups or within subgroup contrasts, respectively. To compare points on an effect curve to a single point, use [reference_curve()].
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [reference_curve()] for comparing points on an effect curve to a single point
#' * [summary.curve_est()] for inference on individual points on an effect curve
#' * [marginaleffects::hypotheses()] for general hypotheses on `curve_est` (and other) objects
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- lm(Math ~ poly(logBLL, 5) *
#'             (Male + Age + Race + PIR +
#'                Enough_Food),
#'           data = nhanes3lead)
#'
#' # ADRF of logBLL on Math, unconditional
#' # inference
#' adrf1 <- adrf(fit, treat = "logBLL")
#'
#' # Differences among ADRF estimates at given points
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   point_contrast() |>
#'   summary()

#' @export
point_contrast <- function(object) {
  .chk_is(object, "curve_est")

  .treat <- attr(object, ".treat")
  .vcov_type <- .attr(object, ".vcov_type")
  .contrast <- attr(object, ".contrast")
  .by_grid <- attr(object, ".by_grid")

  .est <- stats::coef(object)

  values <- .attr(object, "values")

  n <- length(values)

  if (n < 2L) {
    fn_name <- rlang::current_call() |> rlang::call_name()
    .err("{.fun {fn_name}} can only be used when the `effect_curve` object was called with two or more values")
  }

  combos <- utils::combn(seq_len(n), 2L, simplify = TRUE)

  nc <- ncol(combos)

  contr_mat <- matrix(0, nrow = nc, ncol = n)

  contr_mat[cbind(seq_len(nc), combos[1L, ])] <- -1
  contr_mat[cbind(seq_len(nc), combos[2L, ])] <- 1

  n_by <- get_n_by(.contrast, .by_grid)

  contr_mat <- block_diag(contr_mat, n_by)

  if (.vcov_type == "bootstrap") {
    .boot <- .attr(object, ".boot")

    .boot[["t"]] <- .boot[["t"]] |>
      tcrossprod(contr_mat)

    .boot[["t0"]] <- .boot[["t0"]] |>
      tcrossprod(contr_mat) |>
      drop()

    attr(object, ".vcov") <- stats::vcov(.boot)
    attr(object, ".boot") <- .boot
  }
  else if (.vcov_type == "posterior") {
    .draws <- .attr(object, ".draws") |>
      tcrossprod(contr_mat)

    attr(object, ".vcov") <- stats::cov(.draws)
    attr(object, ".draws") <- .draws
  }
  else if (.vcov_type != "none") {
    .vcov <- stats::vcov(object)
    attr(object, ".vcov") <- quad_mult(contr_mat, .vcov)
  }


  .contrast_names <- sprintf("[%s = %s] - [%s = %s]",
                             .treat, values[combos[2L, ]],
                             .treat, values[combos[1L, ]])

  .est0 <- drop(contr_mat %*% .est)

  object[setdiff(names(object), c(c(.treat, "estimate")))] <- NULL
  setrename(object, setNames("term", .treat))

  object <- object |>
    ss(rep.int(1L, length(.est0))) |>
    ftransform(estimate = .est0) |>
    add_est_labels(.contrast, .by_grid, .contrast_names, "term")

  attr(object, ".contrast_names") <- .contrast_names

  class(object) <- unique(c("curve_est_contrast", class(object)))

  object
}

#' @exportS3Method summary curve_est_contrast
#' @rdname point_contrast
summary.curve_est_contrast <- function(object, conf_level = 0.95, simultaneous = TRUE, null = 0,
                                       ci.type = "perc", df = NULL, ...) {

  .vcov_type <- .attr(object, ".vcov_type")

  if (.vcov_type == "none") {
    conf_level <- 0
    null <- NA

    do_inference <- FALSE
  }
  else {
    conf_level <- process_conf_level(conf_level)
    null <- process_null(null)

    do_inference <- (conf_level > 0 || !anyNA(null)) &&
      is_not_null(stats::vcov(object))
  }

  if (!do_inference) {
    class(object) <- unique(c("summary.curve_est_contrast",
                              "summary.curve_est",
                              class(object)))

    return(object)
  }

  .boot <- .attr(object, ".boot")
  .draws <- .attr(object, ".draws")

  est <- stats::coef(object)

  simultaneous <- process_simultaneous(simultaneous, est)

  ci.type <- process_ci.type(ci.type, .vcov_type, simultaneous)

  if (ci.type == "wald") {
    df <- df %or% .attr(object, ".df")

    chk::chk_number(df)
    chk::chk_gt(df, 0)

    stat <- if (is.finite(df)) "t" else "z"

    vcov <- stats::vcov(object)

    se <- sqrt(diag(vcov))
    zeros <- se < 1e-10

    res_names <- c("std.error",
                   stat[!anyNA(null)],
                   "p.value"[!anyNA(null)],
                   "conf.low"[conf_level > 0],
                   "conf.high"[conf_level > 0])

    res <- make_df(res_names, length(est))

    res$std.error <- se
  }
  else {
    res_names <- c("p.value"[!anyNA(null)],
                   "conf.low"[conf_level > 0],
                   "conf.high"[conf_level > 0])

    res <- make_df(res_names, length(est))
  }

  t_crit <- vp <- NULL

  # Compute p-values
  if (!anyNA(null)) {
    if (ci.type == "wald") {
      # Compute test statistic on transformed estimates if null is specified
      res[[stat]][!zeros] <- (est[!zeros] - null) / se[!zeros]

      if (simultaneous) {
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
      }
      else {
        pvals <- switch(stat,
                        "z" = 2 * pnorm(-abs(res[[stat]][!zeros])),
                        "t" = 2 * pt(-abs(res[[stat]][!zeros]), df))
      }

      res$p.value[!zeros] <- pvals
    }
    else if (.vcov_type == "bootstrap") {
      s <- summary(.boot,
                   conf = 0,
                   parm = seq_along(est),
                   ci.type = ci.type,
                   p.value = TRUE,
                   null = null,
                   simultaneous = simultaneous)

      res$p.value <- s[, ncol(s)]
    }
    else if (.vcov_type == "posterior") {
      res$p.value <- posterior_p_value(.draws, parm = seq_along(est),
                                       null = null, simultaneous)
    }
  }

  # Compute CI
  if (conf_level > 0) {
    if (ci.type == "wald") {
      if (simultaneous) {
        if (is_null(vp)) {
          vp <- ss(vcov, !zeros, !zeros) |>
            cov2cor() |>
            .to_psd()
        }

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
      res$conf.low[] <- est - t_crit * se
      res$conf.high[] <- est + t_crit * se

    }
    else if (.vcov_type == "bootstrap") {
      boot_ci <- confint(.boot, level = conf_level,
                         parm = seq_along(est),
                         ci.type = ci.type,
                         simultaneous = simultaneous)

      res$conf.low[] <- boot_ci[, 1L]
      res$conf.high[] <- boot_ci[, 2L]
    }
    else if (.vcov_type == "posterior") {
      draw_ci <- posterior_ci(.draws, level = conf_level,
                              parm = seq_along(est),
                              simultaneous = simultaneous)

      res$conf.low[] <- draw_ci[, 1L]
      res$conf.high[] <- draw_ci[, 2L]
    }
  }

  for (i in names(res)) {
    object[[i]] <- res[[i]]
  }

  attr(object, "simultaneous") <- simultaneous
  attr(object, "conf_level") <- conf_level

  if (conf_level > 0 && is_not_null(t_crit)) {
    attr(t_crit, "df") <- df
    attr(object, "crit") <- t_crit
  }

  attr(object, "ci.type") <- ci.type

  if (!anyNA(null)) {
    attr(object, "null") <- null
  }

  class(object) <- unique(c("summary.curve_est_contrast",
                            "summary.curve_est",
                            class(object)))

  object
}
