#' Compute points on an effect curve
#'
#' @description
#' `summary()` computes estimates and confidence intervals for specified points on the supplied effect curve.
#'
#' @param object a [`curve_est`] object; the output of a call to an [`effect_curve`] object.
#' @param conf_level the desired confidence level. Set to 0 to omit confidence intervals. Default is .95.
#' @param simultaneous `logical`; whether the computed p-values and confidence intervals should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known as uniform) intervals jointly cover all specified estimates at the desired confidence level, whereas pointwise confidence intervals only cover each estimate at the desired level. Simultaneous p-values are inversions of the simultaneous confidence intervals. Default is `TRUE`. See Details.
#' @param null the null value for the hypothesis tests. Default is to use a null value of 0 when the effect curve is an AMEF, a curve contrast, or a reference effect curve, and to omit hypothesis tests otherwise. Set to `NA` to manually omit hypothesis tests.
#' @param transform whether to compute intervals and perform tests on the transformed estimates. Allowable options include `TRUE`, `FALSE`, a family object (e.g., `binomial()`), a link object (e.g., the output of [make.link()]), or one of the functions `log`, `sqrt`, `exp`, `qlogis`, or `qnorm`. Ignored unless `object` is an ADRF and Wald inference is used (i.e., ignored with percentile intervals, the default with bootstrapping and Bayesian models). See Details.
#' @param df the "denominator" degrees of freedom to use for the tests and critical test statistics for confidence intervals. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise, or, for models fit to multiply imputed data, degrees of freedom that account for the number of imputations (see [adrf()]).
#' @param ci.type string; when bootstrapping or Bayesian inference is used in the original effect curve, which type of confidence interval is to be computed. For bootstrapping, allowable options include `"perc"` for percentile intervals, `"wald"` for Wald intervals, and other options allowed by \pkgfun{fwb}{summary.fwb}. When `simultaneous = TRUE`, only `"perc"` and `"wald"` are allowed. For Bayesian models, allowable options include `"perc"` for equi-tailed intervals and `"wald"` for Wald intervals. Default is `"perc"`. Ignored when bootstrapping is not used and the model is not Bayesian.
#' @param ... ignored.
#'
#' @returns
#' `summary()` returns an object of class `summary.curve_est`, which inherits from [`curve_est`]. This is a `data.frame` with columns for the treatment, estimates, and uncertainty measures (p-values, confidence intervals, etc.).
#'
#' @details
#' ## Transform
#'
#' The usual confidence intervals and tests assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. `transform` transforms the estimates to ones that are unbounded and computes the corresponding distribution of transformed estimates using the delta method. By default, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simply a transformation of the estimated points along the curve already computed. Confidence intervals are computed using the transformed estimates before being back-transformed to ensure they are within the bounds of the outcome. When `null` is a number, that number is also transformed. When a transformation other than the identity is used, including the default link function of a generalized linear model, standard errors are not reported, because the standard errors used for tests and confidence intervals are those of the transformed estimates. `transform` is only used with Wald inference; percentile intervals, the default with bootstrapping and Bayesian models, already respect the bounds of the outcome because they are invariant to transformation.
#'
#' When `transform` is `TRUE`, the link function of the outcome model is used; when a family or link object is supplied, its link function is used. Only transformations with a known inverse and derivative are allowed. A value of `null` outside the domain of the transformation (e.g., a null value of 0 or 1 with a logit transformation) gives an error; set `transform = FALSE` to test it on the original scale. A transformation supplied to `transform` that cannot be applied to every estimate, or, with Wald intervals, to every bootstrap or posterior draw, also gives an error. When `transform = TRUE` and the link function of the outcome model cannot be applied to every draw (e.g., a draw between grid points just outside 0 and 1), the estimates are left untransformed with a warning.
#'
#' ## Simultaneous confidence intervals and tests
#'
#' Simultaneous confidence intervals ensure all estimates, not just a given individual point, are contained within their confidence intervals at the given confidence level. These are wider than pointwise intervals to reflect that they are covering multiple estimates, which otherwise would decrease the true coverage rate from that specified. `summary()` uses the "sup-t" simultaneous confidence interval, which is the smallest one-parameter interval that covers all estimates at the desired rate. Simultaneous hypothesis tests are performed by inverting the simultaneous confidence intervals; the p-value for each test is the complement of the smallest confidence level for which a simultaneous confidence interval accounting for the other tests contains the null value. The widths of the confidence intervals and the p-values depend on how many and which estimates are computed. With Wald inference, the simultaneous critical value and p-values are computed from the multivariate t (or normal) distribution of the estimates using \pkgfun{mvtnorm}{qmvt} and \pkgfun{mvtnorm}{pmvt}, which use a randomized algorithm, so they can differ slightly between calls; use [set.seed()] for reproducible results. With bootstrapping or Bayesian models and `ci.type = "perc"`, the simultaneous intervals are the equal-tailed percentile intervals with the tail probability calibrated so that they jointly contain the desired proportion of the bootstrap or posterior draws, as described by Montiel Olea and Plagborg-Møller (2019).
#'
#' @references
#' Montiel Olea, J. L., & Plagborg-Møller, M. (2019). Simultaneous confidence bands: Theory, implementation, and an application to SVARs. *Journal of Applied Econometrics*, 34(1), 1–17. \doi{10.1002/jae.2656}
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [`curve_est`] for information on the output of an effect curve
#' * [plot.effect_curve()] for plotting an effect curve
#' * [summary.effect_curve()] for testing omnibus hypotheses about an effect curve
#'
#' @examples
#' data("nhanes3lead")
#'
#' fit <- glm(Block >= 12 ~ poly(logBLL, 3) *
#'              Male * (Age + Race + PIR + NICU +
#'                        Smoke_Pregnant),
#'            data = nhanes3lead,
#'            family = binomial)
#'
#' # ADRF of logBLL on P(Block >= 12)
#' adrf1 <- adrf(fit, treat = "logBLL")
#'
#' # Estimates along ADRF with simultaneous CIs computed
#' # from transformed estimates
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary()
#'
#' # Estimates along ADRF with pointwise CIs computed
#' # from transformed estimates
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary(simultaneous = FALSE)
#'
#' # Estimates along ADRF with simultaneous CIs computed
#' # from original estimates
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary(transform = FALSE)
#'
#' # Estimates along ADRF with simultaneous CIs computed
#' # from transformed estimates, hypothesis tests against
#' # null of .1
#' adrf1(logBLL = c(0, 1, 2)) |>
#'   summary(null = .1)

#' @exportS3Method summary curve_est
summary.curve_est <- function(object, conf_level = 0.95, simultaneous = TRUE, null = NULL,
                              transform = TRUE, ci.type = "perc", df = NULL, ...) {

  .vcov_type <- attr(object, ".vcov_type")

  if (.vcov_type == "none") {
    do_inference <- FALSE
  }
  else {
    conf_level <- process_conf_level(conf_level)
    null <- process_null(null, object)

    do_inference <- (conf_level > 0 || !anyNA(null)) &&
      !is_null(stats::vcov(object))
  }

  if (!do_inference) {
    class(object) <- unique(c("summary.curve_est", class(object)))

    return(object)
  }

  .boot <- .attr(object, ".boot")
  .draws <- .attr(object, ".draws")

  est0 <- stats::coef(object)

  simultaneous <- process_simultaneous(simultaneous, est0)

  ci.type <- process_ci.type(ci.type, .vcov_type, simultaneous)

  # transform; only Wald inference uses the transformation, so only then are
  # the null value and the bootstrap or posterior draws transformed
  transform_list <- process_transform(transform, object, est0,
                                      null = if (ci.type == "wald") null,
                                      check_draws = ci.type == "wald")

  # First, apply transforms if necessary
  transform <- transform_list$transform
  inv_transform <- transform_list$inv_transform

  est <- transform(est0)

  do_transform <- ci.type != "perc" &&
    !all(check_if_zero(est - est0))

  if (ci.type == "wald") {
    df <- df %or% .attr(object, ".df")

    arg::arg_number(df)
    arg::arg_gt(df, 0)

    stat <- if (is.finite(df)) "t" else "z"

    if (.vcov_type == "bootstrap") {
      vcov <- .boot[["t"]] |>
        ss(j = seq_along(est)) |>
        transform() |>
        cov()
    }
    else if (.vcov_type == "posterior") {
      vcov <- .draws |>
        ss(j = seq_along(est)) |>
        transform() |>
        cov()
    }
    else if (do_transform) {
      d_transform <- transform_list$d_transform
      # nrow is given because diag() of a single number makes an identity matrix
      vcov <- quad_mult(diag(d_transform(est0), nrow = length(est0)),
                        stats::vcov(object))
    }
    else {
      vcov <- stats::vcov(object)
    }

    se <- sqrt(diag(vcov))
    zeros <- se < 1e-10

    res_names <- c("std.error"[!do_transform],
                   stat[!anyNA(null)],
                   "p.value"[!anyNA(null)],
                   "conf.low"[conf_level > 0],
                   "conf.high"[conf_level > 0])

    res <- make_df(res_names, length(est0))

    if (!do_transform) {
      res$std.error <- se
    }
  }
  else {
    res_names <- c("p.value"[!anyNA(null)],
                   "conf.low"[conf_level > 0],
                   "conf.high"[conf_level > 0])

    res <- make_df(res_names, length(est0))
  }

  t_crit <- vp <- NULL

  # Compute p-values
  if (!anyNA(null)) {
    if (ci.type == "wald") {
      # Compute test statistic on transformed estimates and null
      res[[stat]][!zeros] <- (est[!zeros] - transform(null)) / se[!zeros]

      if (simultaneous) {
        if (is_null(vp)) {
          vp <- ss(vcov, !zeros, !zeros) |>
            cov2cor() |>
            .to_psd()
        }

        pvals <- .mvt_p_values(res[[stat]][!zeros], vp, df)
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
                   index = seq_along(est),
                   ci.type = ci.type,
                   p.value = TRUE,
                   null = null,
                   simultaneous = simultaneous)

      res$p.value <- s[, fncol(s)]

      # Estimates with no variance have no p-value, as with Wald inference
      is.na(res$p.value[.zero_var(ss(.boot[["t"]], j = seq_along(est)))]) <- TRUE
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

        t_crit <- .mvt_crit(conf_level, vp, df)
      }
      else {
        t_crit <- switch(stat,
                         "z" = abs(qnorm((1 - conf_level) / 2)),
                         "t" = abs(qt((1 - conf_level) / 2, df)))
      }

      # Reverse transformation
      res$conf.low[] <- inv_transform(est - t_crit * se)
      res$conf.high[] <- inv_transform(est + t_crit * se)

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

  # Use original estimate in output
  object$estimate <- {
    if (do_transform) inv_transform(est)
    else est0
  }

  attr(object, "simultaneous") <- simultaneous
  attr(object, "conf_level") <- conf_level

  if (conf_level > 0 && !is_null(t_crit)) {
    attr(t_crit, "df") <- df
    attr(object, "crit") <- t_crit
  }

  attr(object, "ci.type") <- ci.type

  if (!anyNA(null)) {
    attr(object, "null") <- null
  }

  class(object) <- unique(c("summary.curve_est", class(object)))

  object
}

#' @exportS3Method coef curve_est
#' @rdname summary.curve_est
coef.curve_est <- function(object, ...) {
  object[["estimate"]] |>
    setNames(.get_curve_est_labels(object))
}

#' @exportS3Method vcov curve_est
#' @rdname summary.curve_est
vcov.curve_est <- function(object, ...) {
  v <- .attr(object, ".vcov")

  # No covariance matrix when vcov = "none" in the original call to adrf()
  if (is_null(v)) {
    return(NULL)
  }

  dimnames(v) <- rep.int(list(.get_curve_est_labels(object)), 2L)

  v
}

#' @exportS3Method print summary.curve_est
print.summary.curve_est <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  arg::arg_whole_number(digits)

  .reference <- .attr(x, ".reference")
  .contrast <- .attr(x, ".contrast")

  label <- {
    if (inherits(x, "summary.curve_est_contrast")) "Point Contrasts"
    else "Estimates"
  }

  main <- paste(c(get_curve_type(x),
                  if (!is_null(.reference)) "Reference",
                  if (!is_null(.contrast)) "Contrast",
                  label),
                collapse = " ")

  .print_estimate_table(x = x, digits = digits, topn = Inf,
                        main = main, rownames = FALSE,
                        ...)

  .print_inference(x)

  invisible(x)
}

.get_curve_est_labels <- function(x) {
  .est_names <- .attr(x, ".proj_coefficient_names") %or%
    .attr(x, ".contrast_names") %or%
    sprintf("%s = %s", .attr(x, ".treat"), .attr(x, "values"))

  .contrast <- .attr(x, ".contrast")
  .by_grid <- .attr(x, ".by_grid")

  if (is_null(.by_grid) && is_null(.contrast)) {
    return(.est_names)
  }

  .labels <- .contrast %or% get_by_grid_labels(.by_grid)

  sprintf("%s | %s",
          rep(.labels, each = length(.est_names)),
          rep.int(.est_names, length(.labels)))
}

.print_inference <- function(x, .print = TRUE) {
  .vcov_type <- .attr(x, ".vcov_type")
  null <- .attr(x, "null")
  conf_level <- .attr(x, "conf_level")

  out <- character(0L)

  do_inference <- (.vcov_type != "none") &&
    ((!is_null(conf_level) && conf_level > 0) ||
       (!is_null(null) && !anyNA(null))) &&
    !is_null(stats::vcov(x))

  if (do_inference) {
    simultaneous <- .attr(x, "simultaneous")

    out <- c(out,
             .it(sprintf("Inference: %s", .get_vcov_type_name(.vcov_type, simultaneous))))

    #Additional info
    if (conf_level > 0) {
      .crit <- .attr(x, "crit")

      if (is_null(.crit)) {
        stat_msg <- ""
      }
      else {
        df <- .attr(.crit, "df")

        stat_msg <- {
          if (is.finite(df)) {
            sprintf(" (t* = %s, df = %s)",
                    round(.crit, 3),
                    round(df, 1))
          }
          else {
            sprintf(" (z* = %s)",
                    round(.crit, 3))
          }
        }
      }

      out <- c(out,
               .it(sprintf("Confidence level: %s%%%s",
                           round(100 * conf_level, 2L),
                           stat_msg)))
    }
  }

  .reference <- .attr(x, ".reference")
  .treat <- .attr(x, ".treat")

  null_ref <- c(if (!is_null(.reference)) sprintf("Reference: %s = %s",
                                                     .treat, .reference),
                if (!is_null(null)) sprintf("Null value: %s",
                                               format(null, digits = 4L, drop0trailing = TRUE)))

  if (!is_null(null_ref)) {
    out <- c(out, paste(.it(null_ref), collapse = " | "))
  }

  if (.print && !is_null(out)) {
    cli::cat_line(out)
  }

  invisible(paste(out, collapse = "\n"))
}
