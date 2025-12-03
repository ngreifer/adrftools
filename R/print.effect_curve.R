#' @exportS3Method print effect_curve
print.effect_curve <- function(x, ...) {
  cat(sprintf("An %s object\n\n",
              .it("effect_curve")))

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

  width <- min(60, getOption("width", 60))

  sprintf(
    "{.emph Use {.help [{.fun plot}](%1$s::plot.effect_curve)} to plot the curve, {.help [{.fun summary}](%1$s::summary.effect_curve)} to test the curve, or {.help [{'`{object}(values)`'}](%1$s::effect_curve-class)} to compute estimates.}",
    utils::packageName()
  ) |>
    cli::format_inline() |>
    cli::ansi_strwrap(width) |>
    cat()

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
  cli::cat_line(" - curve type: ", .get_full_curve_type_name(x))
}

.print_response <- function(x) {
  cli::cat_line(" - response: ", .attr(x, ".response"))
}

.print_treat <- function(x) {
  cli::cat_line(sprintf(" - treatment: %s", .attr(x, ".treat")))
}

.print_range <- function(x) {
  .values <- .attr(x, ".values")

  cli::cat_line(sprintf("   + range: %s to %s",
                        round(min(.values), 4),
                        round(max(.values), 4)))
}

.print_reference <- function(x) {
  if (inherits(x, "reference_curve")) {
    .reference <- .attr(x, ".reference")

    if (is_not_null(.reference)) {
      cli::cat_line(" - reference level: ", round(.reference, 4L))
    }
  }
}

.print_by <- function(x) {
  if (!inherits(x, "contrast_curve")) {
    .by_grid <- .attr(x, ".by_grid")

    if (is_not_null(.by_grid)) {
      cli::cat_line(" - by: ", toString(names(.by_grid)))
    }
  }
}

.print_contrast <- function(x) {
  if (inherits(x, "contrast_curve")) {
    .contrast <- .attr(x, ".contrast")

    if (is_not_null(.contrast)) {
      cli::cat_line(sprintf(" - %s: %s",
                            ngettext(length(.contrast), "contrast", "contrasts"),
                            toString(.contrast)))
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
  cli::cat_line(" - inference: ", .get_vcov_type_name(.attr(x, ".vcov_type")))
}

.print_projection <- function(x) {
  if (inherits(x, "curve_projection")) {
    .formula <- .attr(x, ".proj_formula")

    if (is_not_null(.formula)) {
      cli::cat_line(" - projection model: ", .formula)
    }
  }
}
