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

test_that("as.data.frame() returns the estimates and bands plot() displays", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 6, by = ~Male)

  set.seed(1)
  d <- as.data.frame(a)
  set.seed(1)
  rb <- ribbon_data(plot(a))

  expect_s3_class(d, "data.frame")
  expect_named(d, c("Male", "logBLL", "estimate", "conf.low", "conf.high"))
  expect_identical(nrow(d), 12L)
  expect_identical(d$Male, rep(c(0L, 1L), each = 6))
  expect_identical(d$logBLL, rep(get_values(a), 2))
  expect_equal(d$estimate, get_est(a), tolerance = 1e-12)

  # The bands match the plotted ribbon, which is drawn per subgroup
  expect_equal(sort(d$conf.low), sort(rb$ymin), tolerance = 1e-10)
  expect_equal(sort(d$conf.high), sort(rb$ymax), tolerance = 1e-10)

  # Pointwise bands match summary() at the grid points
  d_p <- as.data.frame(a, simultaneous = FALSE, subset = Male == 1)
  s <- summary(a(logBLL = get_values(a), subset = Male == 1), simultaneous = FALSE)
  expect_identical(unique(d_p$Male), 1L)
  expect_equal(d_p$conf.low, s$conf.low, tolerance = 1e-8)
  expect_equal(d_p$conf.high, s$conf.high, tolerance = 1e-8)

  # The projection is added as a column
  p <- curve_projection(a, "linear")
  d_proj <- as.data.frame(a, proj = p, simultaneous = FALSE)
  expect_equal(d_proj$proj_estimate, attr(p, ".est"), tolerance = 1e-12)

  # row.names are applied
  expect_identical(rownames(as.data.frame(a, row.names = letters[1:12], simultaneous = FALSE)),
                   letters[1:12])
})

test_that("as.data.frame() works for contrasts, curves without inference, and shift curves", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  cc <- curve_contrast(adrf(fit, treat = "logBLL", n = 4, by = ~Male))
  d_cc <- as.data.frame(cc, simultaneous = FALSE)
  expect_named(d_cc, c("contrast", "logBLL", "estimate", "conf.low", "conf.high"))
  expect_identical(unique(d_cc$contrast), "[Male = 1] - [Male = 0]")
  expect_equal(d_cc$estimate, get_est(cc), tolerance = 1e-12)

  a0 <- adrf(fit, treat = "logBLL", n = 4, vcov = "none")
  expect_named(as.data.frame(a0), c("logBLL", "estimate"))

  sc <- shift_curve(fit, treat = "logBLL", n = 5)
  d_sc <- as.data.frame(sc, simultaneous = FALSE)
  expect_named(d_sc, c("shift", "estimate", "conf.low", "conf.high"))
  expect_equal(d_sc$shift, get_values(sc), tolerance = 1e-12)
})
