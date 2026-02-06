#' Plot an effect curve
#'
#' @description
#' `plot()` plots an effect curve and its confidence band.
#'
#' @param x an `effect_curve` object; the output of a call to [adrf()] or functions that modify it.
#' @param conf_level the desired confidence level. Default is .95 for 95% confidence bands.
#' @param simultaneous `logical`; whether the computed confidence bands should be simultaneous (`TRUE`) or pointwise (`FALSE`). Simultaneous (also known as uniform) bands cover the full line at the desired confidence level, whereas pointwise confidence bands only cover each point at the desired level. Default is `TRUE`. See Details.
#' @param null the value at which to plot a horizontal reference line. Default is to plot a line with a y-intercept of 0 when the effect curve is an AMEF, a curve contrast, or a reference effect curve, and to omit a line otherwise. Set to `NA` to manually omit the line.
#' @param subset an optional logical expression indicating the subset of the subgroups to plot. Can only be used when `by` was supplied to the original call to [adrf()], and only to refer to variables defining subgroups.
#' @param proj an optional `curve_projection` object, the output of a call to [curve_projection()]. Supplying this adds the projection curve to the effect curve plot.
#' @param transform whether to compute intervals on the transformed estimates. Allowable options include `TRUE`, `FALSE`, or a function specifying a transformation. Ignored unless `x` is an ADRF. See Details.
#' @param df the "denominator" degrees of freedom to use for the critical test statistics for confidence bands. Default is to use the residual degrees of freedom from the original model if it is a linear model and `Inf` otherwise.
#' @param ci.type string; when bootstrapping is used in the original effect curve, what type of confidence interval is to be computed. Allowable options include `"perc"` for percentile intervals, `"wald"` for Wald intervals, and other options allowed by \pkgfun{fwb}{summary.fwb}. When `simultaneous = TRUE`, only `"perc"` and `"wald"` are allowed. Default is `"perc"`. Ignored when bootstrapping is not used.
#' @param ... ignored.
#'
#' @returns
#' A `ggplot` object that can modified using functions in \pkg{ggplot2}. Below are some common tasks and instructions to perform them. Note all should be run after running `library("ggplot2")`.
#'
#' \itemize{
#' \item{Change the position of the legend:
#' \preformatted{
#' theme(legend.position = "{POSITION}")}
#' }
#' \item{Remove the legend:
#' \preformatted{
#' theme(legend.position = "none")}
#' }
#' \item{Change the color of the plotted line:
#' \preformatted{
#' theme(geom = element_geom(ink = "{COLOR}")}
#' }
#' \item{Change the color scheme of the plotted lines:
#' \preformatted{
#' scale_color_brewer(aesthetics = c("color", "fill"),
#'                    palette = "{PALETTE}")}
#' }
#' \item{Change the title, subtitle, or axis labels:
#' \preformatted{
#' labs(title = "{TITLE}", subtitle = "{SUBTITLE}",
#'      x = "{X-AXIS}", y = "{Y-AXIS}")}
#' }
#' \item{Change the y-axis range:
#' \preformatted{
#' coord_cartesian(ylim = c({LOWER}, {UPPER}),
#'                 expand = FALSE)}
#' }
#' }
#'
#' Strings in brackets are to be changed by the user. Refer to the \CRANpkg{ggplot2} documentation for other options.
#'
#' @details
#' `plot()` displays the effect curve in a plot. The solid line corresponds to the effect curve and the ribbon around it corresponds to its confidence band. When `null` is not `NA`, an additional flat line at `null` is displayed. When `proj` is supplied, a dashed line corresponding to the projection is added.
#'
#' When `by` is supplied to [adrf()], `plot()` produces an effect curve plot for each subgroup. When `x` is the output of a call to [curve_contrast()], `plot()` produces an effect curve plot for each treatment contrast.
#'
#' ## Transform
#'
#' The usual confidence bands assume the estimates along the effect curve are normally distributed (or t-distributed when `df` is not `Inf`). However, when the outcome is bounded (e.g., a probability bounded between 0 and 1), this assumption may not be valid for the ADRF in finite samples. `transform` transforms the estimates to ones that are unbounded and computes the corresponding distribution of transformed estimates using the delta method. By default, if a generalized linear model is used for the outcome with a non-identity link function, the estimates are transformed by the link function to be on an unbounded scale. Note this is not the same as using the linear predictor for the effect curve; this is simple a transformation of the estimated points along the curve already computed. Confidence bands are computed using the transformed estimates before being back-transformed to ensure they are within the bounds of the outcome.
#'
#' ## Simultaneous confidence bands
#'
#' Simultaneous confidence bands ensure the whole effect curve, not just a given individual point, is contained within the band at the given confidence level. These are wider than pointwise bands to reflect that they are covering multiple estimates, which otherwise would decrease the true coverage rate from that specified. `plot()` uses the "sup-t" simultaneous confidence band, which is the smallest one-parameter band that covers the whole effect curve at the desired rate.
#'
#' @seealso
#' * [adrf()] for computing the ADRF
#' * [summary.effect_curve()] for testing hypotheses about an effect curve
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
#' adrf1 <- adrf(fit, treat = "logBLL",
#'               n = 50) #using 50 to speed up examples
#'
#' # Plot the ADRF; simultaneous inference,
#' # CIs computed on transformed estimates
#' plot(adrf1)
#'
#' # Plot the ADRF; simultaneous inference,
#' # CIs computed on original estimates
#' plot(adrf1, transform = FALSE)
#'
#' # Plot the ADRF; pointwise inference
#' plot(adrf1, simultaneous = FALSE)
#'
#' # ADRF within subgroups
#' adrf2 <- adrf(fit, treat = "logBLL",
#'               by = ~Male, n = 50)
#'
#' # Plot subgroup ADRFs
#' plot(adrf2)
#'
#' # Plot ADRF in one subgroup
#' plot(adrf2, subset = Male == 1)
#'
#' # ADRF contrast
#' adrf_contrast <- curve_contrast(adrf2)
#'
#' plot(adrf_contrast)


#' @exportS3Method plot effect_curve
plot.effect_curve <- function(x, conf_level = 0.95, simultaneous = TRUE, null = NULL,
                              subset = NULL, proj = NULL, transform = TRUE,
                              ci.type = "perc", df = NULL, ...) {
  check_proj(x, proj)

  .est0 <- .attr(x, ".est")
  .vcov <- .attr(x, ".vcov")
  .vcov_type <- .attr(x, ".vcov_type")
  .treat <- .attr(x, ".treat")
  .values <- .attr(x, ".values")
  .boot <- .attr(x, ".boot")
  .draws <- .attr(x, ".draws")
  .by_grid <- .attr(x, ".by_grid")
  .contrast <- .attr(x, ".contrast")

  if (.vcov_type == "none") {
    conf_level <- 0
    transform <- FALSE
    simultaneous <- FALSE
  }
  else {
    conf_level <- process_conf_level(conf_level)
  }

  null <- process_null(null, x)

  transform_list <- process_transform(transform, x, .est0)

  # subset
  .subset <- process_subset_by_grid(substitute(subset),
                                    .by_grid = .by_grid,
                                    .contrast = .contrast)

  if (is_null(.by_grid) || is_not_null(.contrast)) {
    .s <- seq_along(.est0)
  }
  else {
    .s <- which(rep(.subset, each = length(.est0) / nrow(.by_grid)))
    .est0 <- .est0[.s]
  }

  # First, apply transforms if necessary
  transform <- transform_list$transform
  inv_transform <- transform_list$inv_transform

  .est <- transform(.est0)

  do_inference <- (.vcov_type != "none") &&
    is_not_null(.vcov) && (conf_level > 0)

  if (do_inference) {
    simultaneous <- process_simultaneous(simultaneous, .est)

    ci.type <- process_ci.type(ci.type, .vcov_type, simultaneous)

    res_names <- c(.treat, "estimate",
                   "conf.low", "conf.high")

    res <- make_df(res_names, length(.est))

    do_transform <- ci.type != "perc" &&
      !all(check_if_zero(.est - .est0)) &&
      all(check_if_zero(inv_transform(.est) - .est0))

    t_crit <- vp <- NULL

    # Compute CI
    if (conf_level > 0) {
      if (ci.type == "wald") {
        # Process df
        df <- df %or% .attr(x, ".df")

        arg_number(df)
        arg_gt(df, 0)

        stat <- if (is.finite(df)) "t" else "z"

        if (.vcov_type == "bootstrap") {
          .vcov <- .boot[["t"]] |>
            ss(j = .s) |>
            transform() |>
            cov()
        }
        else if (.vcov_type == "posterior") {
          .vcov <- .draws |>
            ss(j = .s) |>
            transform() |>
            cov()
        }
        else if (do_transform) {
          d_transform <- transform_list$d_transform
          .vcov <- quad_mult(diag(d_transform(.est0)), ss(.vcov, .s, .s))
        }
        else {
          .vcov <- ss(.vcov, .s, .s)
        }

        se <- sqrt(diag(.vcov))
        zeros <- se < 1e-10

        if (simultaneous) {
          vp <- ss(.vcov, !zeros, !zeros) |>
            cov2cor() |>
            .to_psd()

          t_crit <- try(mvtnorm::qmvt(conf_level,
                                      tail = "both.tails",
                                      df = df, corr = vp,
                                      keepAttr = FALSE,
                                      maxiter = 1e5,
                                      abseps = 1e-5)$quantile,
                        silent = TRUE)

          if (null_or_error(t_crit)) {
            error <- conditionMessage(.attr(t_crit, "condition"))
            .err("There was an error computing simultaneous confidence intervals:\n{error}")
          }
        }
        else {
          t_crit <- switch(stat,
                           "z" = abs(qnorm((1 - conf_level) / 2)),
                           "t" = abs(qt((1 - conf_level) / 2, df)))
        }
        # Reverse transformation
        res$conf.low[] <- inv_transform(.est - t_crit * se)
        res$conf.high[] <- inv_transform(.est + t_crit * se)
      }
      else if (.vcov_type == "bootstrap") {
        boot_ci <- confint(.boot, level = conf_level,
                           parm = .s,
                           ci.type = ci.type,
                           simultaneous = simultaneous)

        res$conf.low[] <- boot_ci[, 1L]
        res$conf.high[] <- boot_ci[, 2L]
      }
      else if (.vcov_type == "posterior") {
        draw_ci <- posterior_ci(.draws, level = conf_level,
                                parm = .s,
                                simultaneous = simultaneous)

        res$conf.low[] <- draw_ci[, 1L]
        res$conf.high[] <- draw_ci[, 2L]
      }
    }

    if (is_not_null(t_crit)) {
      attr(t_crit, "df") <- df
      attr(conf_level, "crit") <- t_crit
    }
  }
  else {
    res_names <- c(.treat, "estimate")
    res <- make_df(res_names, length(.est0))
  }

  # Use original estimate in output
  res$estimate <- inv_transform(.est)

  if (is_not_null(proj)) {
    res$proj_estimate <- .attr(proj, ".est")[.s]

    proj_lt <- "55"
  }

  # Apply by_grid and contrast labels if necessary
  if (is_not_null(.contrast)) {
    contrast <- factor(rep(.contrast, each = length(.values)),
                       levels = .contrast)

    res[[.treat]] <- rep.int(.values, length(.contrast))

    res <- cbind(contrast = contrast, res)

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$contrast)) +
      labs(color = "Contrast")

    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$contrast),
                    alpha = .2) +
        labs(fill = "Contrast")
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$contrast),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }

    p <- p +
      scale_color_brewer(palette = "Set1",
                         aesthetics = c("color", "fill"),
                         breaks = .contrast,
                         limits = .contrast)
  }
  else if (is_not_null(.by_grid)) {
    res_names <- names(res)
    val_id <- seq_along(.values)
    by_id <- seq_row(.by_grid)

    res$.vi <- rep.int(val_id, sum(.subset))
    res$.by_id <- rep(which(.subset), each = length(val_id))

    res <- merge(res,
                 cbind(.by_grid, .by_id = by_id),
                 by = ".by_id",
                 all.x = TRUE, all.y = FALSE,
                 sort = FALSE) |>
      roworderv(c(".by_id", ".vi"))

    res[[.treat]] <- .values[res$.vi]

    .labels <- get_by_grid_labels(.by_grid, sep = ", ")

    res$.by_id <- factor(res$.by_id, levels = seq_row(.by_grid),
                         labels = .labels)

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate,
                    color = .data$.by_id)) +
      labs(color = "Subgroup")


    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high,
                        fill = .data$.by_id),
                    alpha = .2) +
        labs(fill = "Subgroup")
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate,
                             color = .data$.by_id),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }

    p <- p + scale_color_brewer(palette = "Set1",
                                aesthetics = c("color", "fill"),
                                breaks = .labels[.subset],
                                limits = .labels)
  }
  else {
    res[[.treat]] <- .values

    p <- ggplot(res, aes(x = .data[[.treat]])) +
      geom_line(aes(y = .data$estimate))

    if (do_inference) {
      p <- p +
        geom_ribbon(aes(ymin = .data$conf.low,
                        ymax = .data$conf.high),
                    alpha = .2)
    }

    if (is_not_null(proj)) {
      p <- p + geom_line(aes(y = .data$proj_estimate),
                         linetype = proj_lt,
                         show.legend = FALSE)
    }
  }

  if (!anyNA(null)) {
    p <- p + geom_hline(yintercept = null, color = "black")
  }

  p + labs(x = .treat, y = .get_y_axis(x),
           title = .get_title(x),
           subtitle = .get_subtitle(x, conf_level, simultaneous)) +
    coord_cartesian(expand = c(left = FALSE, right = FALSE),
                    default = TRUE) +
    theme_bw() +
    theme(legend.position = "right",
          legend.direction = "vertical",
          geom = element_geom(ink = "#E41A1C"),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
}

.get_y_axis <- function(x) {
  .response <- .attr(x, ".response")
  .treat <- .attr(x, ".treat")
  .reference <- .attr(x, ".reference")
  .contrast <- .attr(x, ".contrast")

  curve_atom <- switch(get_curve_type(x),
                       ADRF = bquote(E * group("[", .(.response), "]")),
                       AMEF = bquote(over(d ~ E * group("[", .(.response), "]"),
                                          d ~ .(.treat))))

  label <- curve_atom

  if (is_not_null(.reference)) {
    label <- bquote(.(label) ~ .(bquote(~ - ~ .(curve_atom) * group("(", .(.reference), ")"))))
  }

  if (is_not_null(.contrast)) {
    label <- bquote(.(label) ~ .("Difference"))
  }

  label
}

.get_title <- function(x) {
  .reference <- .attr(x, ".reference")
  .contrast <- .attr(x, ".contrast")

  paste(c(get_curve_type(x),
          if (is_not_null(.reference)) "Reference",
          if (is_not_null(.contrast)) "Contrast",
          if (inherits(x, "curve_projection")) "Projection"),
        collapse = " ")
}

.get_subtitle <- function(x, conf_level, simultaneous) {
  .vcov_type <- .attr(x, ".vcov_type")

  do_inference <- (.vcov_type != "none") &&
    is_not_null(conf_level) && conf_level > 0

  if (!do_inference) {
    return(NULL)
  }

  sprintf("%s%% confidence bands [%s]",
          round(100 * conf_level, 2L),
          .get_vcov_type_name(.vcov_type, simultaneous))
}
