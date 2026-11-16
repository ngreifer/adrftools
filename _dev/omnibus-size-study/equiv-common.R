# Shared code for the comparison of omnibus tests whose null hypotheses are
# equivalent (run-equiv.R, analyze-equiv.R): a flat ADRF against a reference
# curve equal to 0 (at three reference values) and against a flat reference
# curve, a flat ADRF against an AMEF equal to 0, and a linear ADRF against a
# flat AMEF. Uses the data-generating process and outcome model of common.R
# with the package's default (HC2) unconditional variance.

source(file.path(Sys.getenv("ADRF_ROOT", "."), "_dev/omnibus-size-study/common.R"), local = TRUE)

# Effects of the treatment on the outcome; under "none" every null is true,
# under "linear" only the linear-ADRF and flat-AMEF nulls are
equiv_effects <- list(
  none = function(a) 0 * a,
  linear = function(a) .1 * a,
  quadratic = function(a) .1 * a^2
)

ref_values <- c(low = grid[1L], mid = grid[11L], high = grid[21L])

one_rep_equiv <- function(n, effect) {
  d <- gen_data(n)
  d$Y <- d$Y + equiv_effects[[effect]](d$A)

  fit <- lm(fmla, data = d)
  a <- adrf(fit, treat = "A", n = n_grid, range = a_range, data = d)
  m <- amef(a)

  p <- function(x, h) summary(x, hypothesis = h, method = "imhof")$p.value

  out <- c(adrf_flat = p(a, "flat"),
           vapply(ref_values, function(a0) p(reference_curve(a, a0), 0), numeric(1L)) |>
             setNames(paste0("ref0_", names(ref_values))),
           ref_flat_low = p(reference_curve(a, ref_values[["low"]]), "flat"),
           amef_0 = p(m, 0),
           adrf_linear = p(a, "linear"),
           amef_flat = p(m, "flat"))

  list(n = n, effect = effect, version = as.character(getNamespaceVersion("adrftools")),
       p = out)
}
