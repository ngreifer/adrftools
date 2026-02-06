#' @exportS3Method print effect_curve
print.effect_curve <- function(x, ...) {
  cli::format_inline("An {.cls effect_curve} object\n") |>
    cli::cat_line()

  .print_curve_type(x)

  .print_response(x)

  .print_treat(x)

  .print_range(x)

  .print_reference(x)

  .print_projection(x)

  .print_by(x)

  .print_contrast(x)

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
  .curve_type_name <- get_curve_type(x)

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
  cli::format_inline(' - treatment: {(.attr(x, ".treat"))}') |>
    cli::cat_line()
}

.print_range <- function(x) {
  .values <- .attr(x, ".values")
  cli::format_inline("   + range: {round(min(.values), 4)} to {round(max(.values), 4)}") |>
    cli::cat_line()
}

.print_reference <- function(x) {
  if (inherits(x, "reference_curve")) {
    .reference <- .attr(x, ".reference")

    if (is_not_null(.reference)) {
      cli::format_inline(" - reference level: {round(.reference, 4L)}") |>
        cli::cat_line()
    }
  }
}

.print_by <- function(x) {
  if (!inherits(x, "contrast_curve")) {
    .by_grid <- .attr(x, ".by_grid")

    if (is_not_null(.by_grid)) {
      by_vec <- cli::cli_vec(names(.by_grid), style = list("vec-last" = ", "))
      cli::format_inline(" - by: {by_vec}") |>
        cli::cat_line()
    }
  }
}

.print_contrast <- function(x) {
  if (inherits(x, "contrast_curve")) {
    .contrast <- .attr(x, ".contrast")

    if (is_not_null(.contrast)) {
      contrast_vec <- cli::cli_vec(.contrast, style = list("vec-last" = ", "))
      cli::format_inline(" - contrast{?s}: {contrast_vec}") |>
        cli::cat_line()
    }
  }
}

.get_vcov_type_name <- function(.vcov_type, simultaneous = NULL) {
  out <- as.character(.vcov_type)

  if (isTRUE(.attr(.vcov_type, "clustered"))) {
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
  cli::format_inline(' - inference: {(.get_vcov_type_name(.attr(x, ".vcov_type")))}') |>
    cli::cat_line()
}

.print_projection <- function(x) {
  if (inherits(x, "curve_projection")) {
    .formula <- .attr(x, ".proj_formula")

    if (is_not_null(.formula)) {
      cli::format_inline(" - projection model: {(.formula)}") |>
        cli::cat_line()
    }
  }
}
