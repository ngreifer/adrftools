#' Contrast Points on Effect Curves
#'
#' @description
#' `point_contrast()` computes point estimates and confidence intervals for the difference between specified points along an effect curve. Interpolation is used to compare points that may not have been specified in the original call to [adrf()] or [amef()].
#'
#' @inheritParams summary.adrf_est
#' @param values a vector of values of the treatment at which estimates are to be compared.
#' @param reference optional; a single number corresponding to the value of the treatment to which all estimates will be compared. If `NULL` (the default), all pairwise comparisons between estimates at `values` will be reported.
#' @param simultaneous `logical`; whether the computed confidence bands and p-values should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known as uniform) bands cover all contrasts at the desired confidence level, whereas pointwise confidence bands only cover each contrast at the desired level. Default is `TRUE`. See [summary.adrf_est()] for details.
#' @param ... ignored.
#'
#' @returns
#' `point_contrast()` returns an object of class `"point_contrast"`, which is a data frame containing the estimates.
#'
#' @details
#'
#' @seealso
#' * [adrf()], [amef()], and [curve_contrast()] for estimating effect curves
#' * [curve_test()] for testing omnibus null hypotheses about effect curves, such as whether they are flat
#' * [curve_proj()] for projecting the effect curve onto a simpler linear model
#'
#' @examples
#'

#' @export
point_contrast <- function(object, values, reference = NULL, conf_level = 0.95, simultaneous = TRUE, null = 0,
                           df = NULL, ci.type = "perc", ...) {

  .chk_is(object, c("adrf_est", "amef_est"))

  # values
  chk::chk_not_missing(values, "`values`")
  values <- process_values(values, n = NULL, treat_var = object$values,
                           strict = TRUE)
  reference <- process_reference(reference, treat_var = object$values,
                                 strict = TRUE)

  vcov_type <- object[["vcov_type"]]

  if (is_null(vcov_type)) {
    conf_level <- 0
    null <- NULL
    simultaneous <- FALSE
  }

  chk::chk_number(conf_level)
  chk::chk_range(conf_level, c(0, 1))

  if (is_not_null(null) && !allNA(null)) {
    chk::chk_number(null)
  }
  else {
    null <- NULL
  }

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  if (do_inference) {
    # Process simultaneous
    if (length(coef(object)) > 1L) {
      chk::chk_flag(simultaneous)
    }
    else {
      simultaneous <- FALSE
    }

    if (vcov_type == "bootstrap") {
      chk::chk_string(ci.type)
      if (simultaneous) {
        ci.type <- match_arg(ci.type, c("perc", "wald"))
      }
    }
    else {
      ci.type <- "wald"
    }

    # Process df
    if (is_null(df)) {
      df <- get_df(.attr(object, "model"))
    }
    else {
      chk::chk_number(df)
      chk::chk_gt(df, 2)
    }
  }

  res <- .point_contrast_internal(object, values, reference,
                                  conf_level, null = null,
                                  simultaneous = simultaneous,
                                  df = df,
                                  ci.type = ci.type)

  class(res) <- c("point_contrast", class(res))

  res
}

#' @exportS3Method print point_contrast
#' @rdname point_contrast
print.point_contrast <- function(x, digits = max(4L, getOption("digits") - 3L), topn = 10L, ...) {
  chk::chk_whole_number(digits)
  chk::chk_count(topn)

  .print_summary_internal(x = x, digits = digits, topn = topn,
                          ...)

  invisible(x)
}

#' @exportS3Method summary point_contrast
summary.point_contrast <- function(object, ...) {
  object
}

.point_contrast_internal <- function(object, values, reference, conf_level = 0.95, null = 0, simultaneous = TRUE,
                                     df = Inf, ci.type = "perc") {

  if (is_null(reference)) {
    values <- sort(unique(values))
    n <- length(values)
    combos <- combn(n, 2L, simplify = FALSE)
  }
  else {
    values <- sort(unique(c(values, reference)))
    n <- length(values)
    ref_ind <- match(reference, values)
    combos <- lapply(seq_len(n)[-ref_ind], c, ref_ind)
  }

  nc <- length(combos)

  # Inference
  vcov_type <- object[["vcov_type"]]

  do_inference <- is_not_null(vcov_type) &&
    is_not_null(vcov(object)) &&
    (conf_level > 0 || is_not_null(null))

  # Prepare interpolation if needed
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

  # contrast_mat0 encodes contrasts between estimates at values
  contrast_mat0 <- matrix(0, nrow = nc, ncol = n)
  for (i in seq_along(combos)) {
    contrast_mat0[i, combos[[i]]] <- c(1, -1)
  }

  by_list <- {
    if (is_null(object$by_grid) && is_null(object$contrast)) list(seq_len(nc))
    else if (is_not_null(object$contrast)) split(seq_len(nc * n), rep(seq_along(object$contrast), each = nc))
    else split(seq_len(nc * n), rep(seq_row(object$by_grid), each = nc))
  }

  contrast_mat <- block_diag(contrast_mat0, length(by_list))

  # est is the contrast estimates
  est <- drop(tcrossprod(est0, contrast_mat))

  if (do_inference) {
    stat <- if (is.finite(df)) "t" else "z"

    if (is_null(null)) {
      res_names <- c("contrast", "estimate", "std.error",
                     stat,
                     "conf.low", "conf.high",
                     "comp.vals")
    }
    else {
      res_names <- c("contrast", "estimate", "std.error",
                     stat, "p.value",
                     "conf.low", "conf.high",
                     "comp.vals")
    }

    res <- make_df(res_names, length(est))

    vcov0 <- quad_mult(val_mat, object$vcov)

    v <- quad_mult(contrast_mat, vcov0)

    se <- sqrt(diag(v))
    zeros <- se < 1e-10

    if (vcov_type == "bootstrap" && ci.type != "wald") {
      rlang::check_installed(c("fwb", "generics"))

      boot <- object$boot
      boot[["t"]] <- tcrossprod(boot[["t"]][, seq_along(coef(object)), drop = FALSE], val_mat) |>
        tcrossprod(contrast_mat)
      boot[["t0"]] <- tcrossprod(boot[["t0"]][seq_along(coef(object))], val_mat) |>
        tcrossprod(contrast_mat)|>
        drop()

      s <- summary(boot,
                   conf = conf_level,
                   ci.type = ci.type,
                   p.value = TRUE,
                   null = 0,
                   simultaneous = simultaneous) |>
        generics::tidy()

      res$p.value <- s$p.value
      res$conf.low <- s$conf.low
      res$conf.high <- s$conf.high
    }
    else {
      res[[stat]][!zeros] <- est[!zeros] / se[!zeros]

      if (simultaneous) {
        vp <- ss(v, !zeros, !zeros) |>
          cov2cor() |>
          .to_psd()

        t_crit <- try(mvtnorm::qmvt(conf_level,
                                    tail = "both.tails",
                                    df = df, corr = vp,
                                    keepAttr = FALSE,
                                    maxiter = 1e5,
                                    abseps = 1e-5)$quantile,
                      silent = TRUE)

        pvals <- try(vapply(abs(res[[stat]][!zeros]),
                            function(zi) {
                              1 - mvtnorm::pmvt(lower = rep.int(-zi, sum(!zeros)),
                                                upper = rep.int(zi, sum(!zeros)),
                                                df = df, corr = vp,
                                                keepAttr = FALSE,
                                                abseps = 1e-5)
                            }, numeric(1L)),
                     silent = TRUE)

        if (null_or_error(pvals) || null_or_error(t_crit)) {
          .err("there was an error computing simultaneous confidence intervals and p-values")
          # use_mvtnorm <- FALSE
        }

      }
      else {
        t_crit <- switch(stat,
                         "z" = abs(qnorm((1 - conf_level) / 2)),
                         "t" = abs(qt((1 - conf_level) / 2, df)))

        pvals <- switch(stat,
                        "z" = 2 * pnorm(-abs(res[[stat]][!zeros])),
                        "t" = 2 * pt(-abs(res[[stat]][!zeros]), df))
      }

      res$conf.low[] <- est - t_crit * se
      res$conf.high[] <- est + t_crit * se
      res$p.value[!zeros] <- pvals
    }

    res$std.error <- se
  }
  else {
    res_names <- c("contrast", "estimate",
                   "comp.vals")

    res <- make_df(res_names, length(est))
  }

  res$estimate <- est

  nm <- vapply(combos, function(combo) {
    sprintf("[%s = %s] - [%s = %s]",
            object[["treat"]], round(values[combo[1L]], 4),
            object[["treat"]], round(values[combo[2L]], 4))
  }, character(1L))

  comp.vals <- lapply(combos, function(combo) values[combo])

  # Apply by_grid and contrast labels if necessary
  if (is_null(object$by_grid) && is_null(object$contrast)) {
    res$contrast <- nm
    res$comp.vals <- comp.vals

    attr(res, "treat") <- object$treat
  }
  else if (is_not_null(object$by_grid)) {
    res_names <- names(res)
    comp_id <- seq_len(nc)
    by_id <- seq_row(object$by_grid)

    res$.ci <- rep.int(comp_id, length(by_id))
    res$.by_id <- rep(by_id, each = nc)

    res <- merge(res,
                 cbind(object$by_grid, .by_id = by_id),
                 by = ".by_id",
                 all.x = TRUE, all.y = FALSE,
                 sort = FALSE)

    res <- ss(res, order(res$.by_id, res$.ci))

    res$contrast <- nm[res$.ci]
    res$comp.vals <- comp.vals[res$.ci]

    by_names <- names(object$by_grid)
    res <- res[c(by_names, res_names)]

    attr(res, "treat") <- object$treat
    attr(res, "by") <- by_names
  }
  else {
    contrast <- factor(rep(object$contrast, each = nc),
                       levels = object$contrast)

    res$contrast <- rep.int(nm, length(object$contrast))
    res$comp.vals <- rep.int(comp.vals, length(object$contrast))

    res <- cbind(contrast = contrast, res)

    attr(res, "treat") <- object$treat
    attr(res, "contrast") <- object$contrast
  }

  attr(res, "obj") <- object
  attr(res, "values") <- values

  if (do_inference) {
    attr(res, "simultaneous") <- simultaneous
    attr(res, "conf_level") <- conf_level

    if (is_not_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(res, "crit") <- t_crit
    }

    if (vcov_type == "bootstrap") {
      attr(res, "ci.type") <- ci.type
    }
  }

  res
}
