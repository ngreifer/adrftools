#' @exportS3Method print effect_curve
print.effect_curve <- function(x, ...) {
  cli::format_inline("An {.cls effect_curve} object\n") |>
    cli::cat_line()

  .print_curve_type(x)

  .print_response(x)

  .print_outcome(x)

  .print_treat(x)

  .print_range(x)

  .print_reference(x)

  .print_projection(x)

  .print_by(x)

  .print_contrast(x)

  .print_sample(x)

  .print_vcov_type(x)

  cli::cat_line()

  if (interactive()) {
    width <- min(60, getOption("width", 60))

    sprintf(
      "{.emph Use {.help [{.fun plot}](%1$s::plot.effect_curve)} to plot the curve, {.help [{.fun summary}](%1$s::summary.effect_curve)} to test the curve, or {.help [{'`{object}(values)`'}](%1$s::effect_curve-class)} to compute estimates.}",
      utils::packageName()
    ) |>
      cli::format_inline() |>
      cli::ansi_strwrap(width) |>
      cli::cat_line()
  }

  invisible(x)
}

.get_full_curve_type_name <- function(x) {
  .curve_type_name <- .curve_type_label(x)

  if (inherits(x, "reference_curve")) {
    .curve_type_name <- paste(.curve_type_name, "reference")
  }

  if (inherits(x, "contrast_curve")) {
    .curve_type_name <- paste(.curve_type_name, "contrast")
  }

  if (inherits(x, "curve_projection")) {
    .curve_type_name <- paste(.curve_type_name, "projection")
  }

  .curve_type_name
}

.print_curve_type <- function(x) {
  cli::format_inline(' - curve type: {(.get_full_curve_type_name(x))}') |>
    cli::cat_line()
}

.print_response <- function(x) {
  cli::format_inline(' - response: {(.attr(x, ".response"))}') |>
    cli::cat_line()
}

.print_treat <- function(x) {
  .shift <- .attr(x, ".shift")

  cli::format_inline(' - treatment: {(.shift[["treat"]] %or% .attr(x, ".treat"))}') |>
    cli::cat_line()

  if (!is_null(.shift)) {
    cli::format_inline(' - shift: {(.shift[["label"]])}') |>
      cli::cat_line()
  }
}

.print_range <- function(x) {
  .values <- .attr(x, ".values")
  sd_units <- .attr(x, ".shift")[["sd_units"]]

  sd_msg <- {
    if (is_null(sd_units)) ""
    else sprintf(" (%s to %s SD)", round(sd_units[1L], 4), round(sd_units[2L], 4))
  }

  cli::format_inline("   + range: {round(min(.values), 4)} to {round(max(.values), 4)}{sd_msg}") |>
    cli::cat_line()

  cli::format_inline("   + grid points: {length(.values)}") |>
    cli::cat_line()
}

.print_reference <- function(x) {
  if (inherits(x, "reference_curve")) {
    .reference <- .attr(x, ".reference")

    if (!is_null(.reference)) {
      cli::format_inline(" - reference level: {round(.reference, 4L)}") |>
        cli::cat_line()
    }
  }
}

.print_by <- function(x) {
  if (!inherits(x, "contrast_curve")) {
    .by_grid <- .attr(x, ".by_grid")
    by_vars <- setdiff(names(.by_grid), .attr(x, ".outcome"))

    if (!is_null(by_vars)) {
      by_vec <- cli::cli_vec(by_vars, style = list("vec-last" = ", "))

      # Subgroups, not counting the categories of the outcome
      n_subgroups <- ss(.by_grid, j = by_vars) |>
        funique() |>
        fnrow()

      cli::format_inline(" - by: {by_vec} ({n_subgroups} subgroup{?s})") |>
        cli::cat_line()
    }
  }
}

.print_outcome <- function(x) {
  .outcome <- .attr(x, ".outcome")

  if (!is_null(.outcome)) {
    categories <- .attr(x, ".by_grid")[[.outcome]] |>
      as.character() |>
      unique()

    cat_vec <- cli::cli_vec(categories, style = list("vec-sep2" = ", ", "vec-last" = ", "))
    cli::format_inline(" - outcome categor{?y/ies}: {cat_vec}") |>
      cli::cat_line()
  }
}

.print_contrast <- function(x) {
  if (inherits(x, "contrast_curve")) {
    .contrast <- .attr(x, ".contrast") |>
      unique()

    if (!is_null(.contrast)) {
      contrast_vec <- cli::cli_vec(.contrast, style = list("vec-last" = ", "))
      cli::format_inline(" - contrast{?s}: {contrast_vec}") |>
        cli::cat_line()
    }
  }
}

# With `details = TRUE` (in print()), the small-sample adjustment, the type of
# bootstrap or of the covariance matrix of the coefficients, and the number of
# bootstrap replicates (`R`) or posterior draws (`n_draws`) are given in
# parentheses along with whether the inference is clustered
.get_vcov_type_name <- function(.vcov_type, simultaneous = NULL, details = FALSE,
                                R = NULL, n_draws = NULL) {
  out <- as.character(.vcov_type)

  if (details) {
    fmt <- function(k) {
      format(k, big.mark = ",", scientific = FALSE)
    }

    info <- c(if (isTRUE(.attr(.vcov_type, "clustered"))) "clustered",
              .attr(.vcov_type, "adjustment"),
              .attr(.vcov_type, "method"),
              if (!is_null(R)) sprintf("%s replicates", fmt(R)),
              if (!is_null(n_draws)) sprintf("%s draws", fmt(n_draws)))

    if (!is_null(info)) {
      out <- sprintf("%s (%s)", out, toString(info))
    }
  }
  else if (isTRUE(.attr(.vcov_type, "clustered"))) {
    out <- paste(out, "(clustered)")
  }

  if (isTRUE(simultaneous)) {
    out <- toString(c(out, "simultaneous"))
  }
  else if (isFALSE(simultaneous)) {
    out <- toString(c(out, "pointwise"))
  }

  out
}

.print_vcov_type <- function(x) {
  .boot <- .attr(x, ".boot")
  .draws <- .attr(x, ".draws")

  vcov_name <- .get_vcov_type_name(.attr(x, ".vcov_type"), details = TRUE,
                                   R = if (!is_null(.boot)) fnrow(.boot[["t"]]),
                                   n_draws = if (!is_null(.draws)) fnrow(.draws))

  cli::format_inline(" - inference: {vcov_name}") |>
    cli::cat_line()
}

.print_sample <- function(x) {
  .sample <- .attr(x, ".sample")

  if (is_null(.sample)) {
    return(invisible(NULL))
  }

  n <- format(range(.sample[["n"]]), big.mark = ",", scientific = FALSE)
  m <- .sample[["m"]]

  n_msg <- {
    if (m == 1L) n[1L]
    else if (n[1L] == n[2L]) sprintf("%s in each imputed dataset", n[1L])
    else sprintf("%s to %s in each imputed dataset", n[1L], n[2L])
  }

  cli::format_inline(" - sample size: {n_msg}") |>
    cli::cat_line()

  if (m > 1L) {
    cli::format_inline(" - imputed datasets: {m}") |>
      cli::cat_line()
  }
}

.print_projection <- function(x) {
  if (inherits(x, "curve_projection")) {
    .formula <- .attr(x, ".proj_formula")

    if (!is_null(.formula)) {
      cli::format_inline(" - projection model: {(.formula)}") |>
        cli::cat_line()
    }

    if (identical(.attr(x, ".proj_integrate"), "empirical")) {
      cli::format_inline(" - weighted by the distribution of {(.attr(x, '.treat'))}") |>
        cli::cat_line()
    }
  }
}
