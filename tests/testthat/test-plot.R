test_that("plot() draws the curve with confidence bands", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  p <- plot(a)

  expect_s3_class(p, "ggplot")
  expect_equal(plot_data(p)$y, get_est(a), tolerance = 1e-8)
  expect_identical(p$labels$title, "ADRF")
  expect_match(p$labels$subtitle, "95% confidence bands \\[unconditional, simultaneous\\]")

  # Simultaneous bands are wider than pointwise bands
  rb_s <- ribbon_data(p)
  rb_p <- ribbon_data(plot(a, simultaneous = FALSE))
  expect_true(all(rb_s$ymax - rb_s$ymin > rb_p$ymax - rb_p$ymin))

  # Pointwise bands match summary() at the grid points
  s <- summary(a(logBLL = get_values(a)), simultaneous = FALSE)
  expect_equal(rb_p$ymin, s$conf.low, tolerance = 1e-8)
  expect_equal(rb_p$ymax, s$conf.high, tolerance = 1e-8)

  # No bands when conf_level = 0
  p0 <- plot(a, conf_level = 0)
  expect_length(p0$layers, 1L)
  expect_null(p0$labels$subtitle)
})

test_that("plot() adds a null line for curves tested against zero", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1L], character(1L))

  expect_false("GeomHline" %in% geoms(plot(a)))
  expect_true("GeomHline" %in% geoms(plot(a, null = 8)))
  expect_true("GeomHline" %in% geoms(plot(amef(a))))
  expect_false("GeomHline" %in% geoms(plot(amef(a), null = NA)))

  expect_identical(plot(amef(a))$labels$title, "AMEF")
  expect_identical(plot(reference_curve(a, 1))$labels$title, "ADRF Reference")
})

test_that("plot() works for subgroups, subsets, and contrasts", {
  fit <- lm(Math ~ poly(logBLL, 2) * Race + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Race)

  p <- plot(a)
  expect_identical(nrow(plot_data(p)), 44L)
  expect_identical(length(unique(plot_data(p)$colour)), 4L)

  p_sub <- plot(a, subset = Race %in% c("White", "Black"))
  expect_identical(nrow(plot_data(p_sub)), 22L)

  pc <- plot(curve_contrast(a), simultaneous = FALSE)
  expect_identical(nrow(plot_data(pc)), 66L)
  expect_identical(pc$labels$title, "ADRF Contrast")

  expect_error(plot(curve_contrast(a), subset = Race == "White"), "cannot be specified")
  expect_error(plot(adrf(fit, treat = "logBLL"), subset = Race == "White"), "by")
})

test_that("plot() overlays a projection", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)
  pr <- curve_projection(a, "linear")

  p <- plot(a, proj = pr)
  expect_identical(length(p$layers), 3L)

  pp <- plot(pr)
  expect_identical(pp$labels$title, "ADRF Projection")

  a_by <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)
  expect_s3_class(plot(a_by, proj = curve_projection(a_by, "linear")), "ggplot")

  cc <- curve_contrast(a_by)
  expect_s3_class(plot(cc, proj = curve_projection(cc, "flat")), "ggplot")

  expect_error(plot(a, proj = curve_projection(a_by, "linear")), "proj")
  expect_error(plot(a, proj = a), "curve_projection")
})

test_that("plot() transforms bands for a logistic ADRF", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)

  rb <- ribbon_data(plot(a, simultaneous = FALSE))
  expect_true(all(rb$ymin > 0 & rb$ymax < 1))

  rb_u <- ribbon_data(plot(a, simultaneous = FALSE, transform = FALSE))
  expect_false(isTRUE(all.equal(rb$ymin, rb_u$ymin)))

  rb_df <- ribbon_data(plot(a, simultaneous = FALSE, transform = FALSE, df = 5))
  expect_true(all(rb_df$ymin < rb_u$ymin))
})

test_that("plot() of an effect curve without uncertainty draws only the curve", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")

  p <- plot(a)
  expect_length(p$layers, 1L)
  expect_null(p$labels$subtitle)
})
