test_that("amef() recovers the slope of a linear ADRF", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  m <- amef(a)

  expect_s3_class(m, "amef_curve")
  expect_false(inherits(m, "adrf_curve"))
  expect_equal(get_est(m), rep(coef(fit)[["logBLL"]], 11), tolerance = 1e-5)

  # The covariate-sampling term cancels in a difference of ADRF values, so the
  # unconditional variance of the slope is its HC0 sandwich variance
  hc0 <- sandwich::vcovHC(fit, type = "HC0")["logBLL", "logBLL"]
  expect_equal(get_se(m), rep(sqrt(hc0), 11), tolerance = 1e-4)

  s <- summary(m(logBLL = c(0, 1)), simultaneous = FALSE)
  expect_identical(attr(s, "null"), 0)
  expect_output(print(s), "AMEF Estimates")
})

test_that("amef() matches the analytic derivative of a quadratic ADRF", {
  fit <- lm(Math ~ logBLL + I(logBLL^2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  m <- amef(a)
  b <- coef(fit)
  vals <- get_values(m)

  # Interior points use central differences and endpoints use one-sided
  # differences
  inner <- 2:20
  expect_equal(get_est(m)[inner], unname(b[2] + 2 * b[3] * vals)[inner], tolerance = 1e-6)
  expect_equal(get_est(m)[-inner], unname(b[2] + 2 * b[3] * vals)[-inner], tolerance = 1e-5)
})

test_that("local polynomial interpolation reproduces a cubic ADRF exactly", {
  fit <- lm(Math ~ logBLL + I(logBLL^2) + I(logBLL^3) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 51)

  vals <- get_values(a)
  b <- coef(fit)
  adrf_true <- function(x) {
    unname(b[1] + b[2] * x + b[3] * x^2 + b[4] * x^3 +
             b[5] * mean(nhanes3lead$Male) + b[6] * mean(nhanes3lead$PIR))
  }

  # Points just inside the ends and between grid points, where the kernel
  # weights span many orders of magnitude
  x <- c(vals[1] + 1e-3, (vals[1] + vals[2]) / 2, vals[25] + 1e-3, vals[51] - 1e-3)
  L <- adrftools:::get_locpoly_w(x = x, v = vals)
  expect_equal(drop(L %*% get_est(a)), adrf_true(x), tolerance = 1e-10)
})

test_that("amef() works with by and validates its input", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  m <- amef(a)
  expect_length(get_est(m), 22L)

  expect_error(amef(m), "amef_curve")
  expect_error(amef(a, eps = -1))
  expect_error(amef(fit), "effect_curve")
  expect_error(amef())
})

test_that("reference_curve() subtracts the estimate at the reference value", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  r <- reference_curve(a, reference = 1)

  expect_s3_class(r, "reference_curve")
  expect_identical(attr(r, ".reference"), 1)
  expect_true(1 %in% get_values(r))
  expect_length(get_values(r), 12L)

  vals <- c(0, .5, 2)
  expect_equal(r(logBLL = vals)$estimate,
               a(logBLL = vals)$estimate - a(logBLL = 1)$estimate,
               tolerance = 1e-8)

  est_r <- r(logBLL = 2)
  est_a <- a(logBLL = c(1, 2))
  L <- matrix(c(-1, 1), nrow = 1)
  expect_equal(unname(vcov(est_r)), unname(L %*% vcov(est_a) %*% t(L)),
               tolerance = 1e-6)

  # A reference already on the grid does not add a point
  r_grid <- reference_curve(a, reference = get_values(a)[3])
  expect_length(get_values(r_grid), 11L)

  expect_output(print(r), "reference level: 1")
  expect_output(print(summary(r(logBLL = 2))), "Reference: logBLL = 1")
})

test_that("reference_curve() validates its input", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)
  r <- reference_curve(a, reference = 1)

  expect_error(reference_curve(a), "reference")
  expect_error(reference_curve(a, reference = 100), "within the range")
  expect_error(reference_curve(a, reference = c(0, 1)))
  expect_error(reference_curve(r, reference = 0), "reference_curve")
})

test_that("reference curves can be built from AMEFs and subgroup curves", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  r <- reference_curve(amef(a), reference = .5)
  expect_s3_class(r, c("reference_curve", "amef_curve"))
  expect_length(get_est(r), 24L)

  est <- r(logBLL = .5)
  expect_equal(est$estimate, c(0, 0), tolerance = 1e-10)
})

test_that("amef() of a reference curve is the AMEF of the original curve", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  m <- amef(reference_curve(a, reference = .5))
  expect_false(inherits(m, "reference_curve"))
  expect_null(attr(m, ".reference"))

  m0 <- amef(a)
  expect_equal(m(logBLL = c(0, .5, 1))$estimate, m0(logBLL = c(0, .5, 1))$estimate,
               tolerance = 1e-8)
  expect_equal(vcov(m(logBLL = c(0, 1))), vcov(m0(logBLL = c(0, 1))), tolerance = 1e-8)

  # A reference AMEF can be built from it
  r <- reference_curve(m, reference = .5)
  expect_equal(r(logBLL = .5)$estimate, 0, tolerance = 1e-10)
})

test_that("curve_contrast() computes differences between subgroup curves", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  cc <- curve_contrast(a)

  expect_s3_class(cc, "contrast_curve")
  expect_identical(attr(cc, ".contrast"), "[Male = 1] - [Male = 0]")

  est <- get_est(a)
  expect_equal(get_est(cc), est[12:22] - est[1:11], tolerance = 1e-10)

  L <- cbind(-diag(11), diag(11))
  expect_equal(get_vcov(cc), unname(L %*% get_vcov(a) %*% t(L)), tolerance = 1e-8)

  e <- cc(logBLL = c(0, 1))
  expect_named(e, c("contrast", "logBLL", "estimate"))
  expect_identical(as.character(e$contrast), rep("[Male = 1] - [Male = 0]", 2))

  s <- summary(e)
  expect_identical(attr(s, "null"), 0)
  expect_output(print(s), "ADRF Contrast Estimates")
  expect_output(print(cc), "contrast: \\[Male = 1\\] - \\[Male = 0\\]")
  expect_identical(names(formals(cc)), "logBLL")
})

test_that("curve_contrast() includes all pairwise comparisons", {
  fit <- lm(Math ~ poly(logBLL, 2) * Race + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Race)

  cc <- curve_contrast(a)
  expect_length(attr(cc, ".contrast"), choose(4, 2))
  expect_length(get_est(cc), 5 * choose(4, 2))
  expect_match(attr(cc, ".contrast")[1], '^\\[Race = "Hispanic"\\] - \\[Race = "Black"\\]$')

  # AMEF contrasts
  cm <- curve_contrast(amef(a))
  expect_s3_class(cm, c("contrast_curve", "amef_curve"))

  # Reference curves of contrasts
  rc <- reference_curve(cc, reference = 1)
  expect_s3_class(rc, c("reference_curve", "contrast_curve"))
  expect_length(get_est(rc), 6 * 6)
})

test_that("curve_contrast() validates its input", {
  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)

  expect_error(curve_contrast(adrf(fit, treat = "logBLL")), "by")

  cc <- curve_contrast(adrf(fit, treat = "logBLL", by = ~Male))
  expect_error(curve_contrast(cc), "contrast_curve")
  expect_error(curve_contrast())
})

test_that("point_contrast() computes all pairwise differences", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  est <- a(logBLL = c(0, 1, 2))
  pc <- point_contrast(est)

  expect_s3_class(pc, "curve_est_contrast")
  expect_named(pc, c("term", "estimate"))
  expect_identical(pc$term, c("[logBLL = 1] - [logBLL = 0]",
                              "[logBLL = 2] - [logBLL = 0]",
                              "[logBLL = 2] - [logBLL = 1]"))

  b <- coef(est)
  expect_equal(pc$estimate, unname(c(b[2] - b[1], b[3] - b[1], b[3] - b[2])),
               tolerance = 1e-10)

  L <- rbind(c(-1, 1, 0), c(-1, 0, 1), c(0, -1, 1))
  expect_equal(unname(vcov(pc)), unname(L %*% vcov(est) %*% t(L)), tolerance = 1e-8)

  s <- summary(pc, simultaneous = FALSE)
  se <- sqrt(diag(vcov(pc)))
  expect_equal(s$std.error, unname(se), tolerance = 1e-10)
  expect_equal(s$p.value, unname(2 * pt(-abs(pc$estimate / se), df.residual(fit))),
               tolerance = 1e-10)
  expect_output(print(s), "ADRF Point Contrasts")

  s_sim <- summary(pc)
  expect_true(all(s_sim$p.value >= s$p.value - 1e-8))

  s_z <- summary(pc, df = Inf, null = NA)
  expect_false("p.value" %in% names(s_z))

  s_noinf <- summary(pc, conf_level = 0, null = NA)
  expect_named(s_noinf, c("term", "estimate"))
})

test_that("point_contrast() of a summary gives the same contrasts as of the estimates", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  est <- a(logBLL = c(0, 1, 2))
  pc <- point_contrast(est)
  pc_s <- point_contrast(summary(est))

  expect_false(inherits(pc_s, "summary.curve_est"))
  expect_null(attr(pc_s, "crit"))
  expect_null(attr(pc_s, "conf_level"))
  expect_equal(as.data.frame(pc_s), as.data.frame(pc), tolerance = 1e-10)
  expect_equal(vcov(pc_s), vcov(pc), tolerance = 1e-10)

  expect_output(print(pc_s), "Point Contrasts")
  expect_equal(summary(pc_s, simultaneous = FALSE)$conf.low,
               summary(pc, simultaneous = FALSE)$conf.low, tolerance = 1e-10)
})

test_that("point_contrast() validates its input", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  expect_error(point_contrast(a(logBLL = 1)), "two or more")
  expect_error(point_contrast(data.frame(x = 1)), "curve_est")

  a_none <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")
  pc_none <- point_contrast(a_none(logBLL = c(0, 1)))
  expect_named(summary(pc_none), c("term", "estimate"))
})

test_that("point_contrast() works within subgroups", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  est <- a(logBLL = c(0, 1, 2))
  pc <- point_contrast(est)

  expect_s3_class(pc, "curve_est_contrast")
  expect_named(pc, c("Male", "term", "estimate"))
  expect_identical(pc$Male, rep(c(0L, 1L), each = 3))

  b <- coef(est)
  L1 <- rbind(c(-1, 1, 0), c(-1, 0, 1), c(0, -1, 1))
  L <- rbind(cbind(L1, matrix(0, 3, 3)), cbind(matrix(0, 3, 3), L1))

  expect_equal(pc$estimate, unname(drop(L %*% b)), tolerance = 1e-10)
  expect_equal(unname(vcov(pc)), unname(L %*% vcov(est) %*% t(L)), tolerance = 1e-8)
  expect_identical(names(coef(pc))[4L], "Male = 1 | [logBLL = 1] - [logBLL = 0]")

  s <- summary(pc, simultaneous = FALSE)
  se <- sqrt(diag(L %*% vcov(est) %*% t(L)))
  expect_equal(s$std.error, unname(se), tolerance = 1e-8)
  expect_output(print(s), "ADRF Point Contrasts")
})
