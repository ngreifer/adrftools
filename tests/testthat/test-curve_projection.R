test_that("projecting a linear ADRF onto a line recovers the slope", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p <- curve_projection(a, "linear")

  expect_s3_class(p, c("curve_projection", "adrf_curve", "effect_curve"))
  expect_named(coef(p), c("(Intercept)", "logBLL"))
  expect_equal(coef(p)[["logBLL"]], coef(fit)[["logBLL"]], tolerance = 1e-8)
  expect_equal(get_est(p), get_est(a), tolerance = 1e-8)

  # The slope's unconditional variance is its HC2 variance
  hc2 <- sandwich::vcovHC(fit, type = "HC2")["logBLL", "logBLL"]
  expect_equal(vcov(p)["logBLL", "logBLL"], hc2, tolerance = 1e-6)

  expect_identical(deparse(formula(p)), "~logBLL")
  mm <- model.matrix(p)
  expect_identical(dim(mm), c(21L, 2L))
  expect_equal(unname(mm[, 2]), get_values(a))

  expect_output(print(p), "projection model: ~logBLL")
})

test_that("projection models can be given by name or formula", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p_q <- curve_projection(a, "quadratic")
  p_f <- curve_projection(a, ~ logBLL + I(logBLL^2))
  expect_equal(unname(coef(p_q)), unname(coef(p_f)), tolerance = 1e-8)
  expect_equal(get_est(p_q), get_est(p_f), tolerance = 1e-8)

  p_poly <- curve_projection(a, ~ poly(logBLL, 2))
  expect_equal(get_est(p_poly), get_est(p_q), tolerance = 1e-8)

  # The flat projection is the trapezoidal average of the curve
  p_flat <- curve_projection(a, "flat")
  vals <- get_values(a)
  w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
  expect_equal(unname(coef(p_flat)), sum(w * get_est(a)) / sum(w), tolerance = 1e-8)
})

test_that("named projection models have the expected number of coefficients", {
  skip_if_not_installed("patrick")

  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  patrick::with_parameters_test_that(
    "model = \"{model}\"",
    {
      expect_length(coef(curve_projection(a, model)), n_coef)
    },
    model = c("flat", "quadratic", "cubic"),
    n_coef = c(1L, 3L, 4L)
  )
})

test_that("curve_projection() validates its input", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  expect_error(curve_projection(a), "model")
  expect_error(curve_projection(a, "wiggly"))
  expect_error(curve_projection(a, Math ~ logBLL), "one-sided")
  expect_error(curve_projection(a, ~ logBLL + PIR), "is allowed to appear")
  expect_error(curve_projection(curve_projection(a, "linear"), "flat"),
               "curve_projection")
})

test_that("summary() of a projection gives coefficient tests", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)
  p <- curve_projection(a, "linear")

  s <- summary(p)

  expect_s3_class(s, "summary.curve_projection")
  expect_named(s, c("term", "estimate", "std.error", "t", "p.value",
                    "conf.low", "conf.high"))
  expect_identical(s$term, c("(Intercept)", "logBLL"))

  se <- sqrt(diag(vcov(p)))
  expect_equal(s$std.error, unname(se), tolerance = 1e-10)
  expect_equal(s$p.value, unname(2 * pt(-abs(coef(p) / se), df.residual(fit))),
               tolerance = 1e-10)
  expect_equal(s$conf.low, unname(coef(p) - qt(.975, df.residual(fit)) * se),
               tolerance = 1e-10)

  expect_equal(coef(s), coef(p))
  expect_equal(vcov(s), vcov(p))

  expect_output(print(s), "ADRF Projection Coefficients")

  s_z <- summary(p, df = Inf, conf_level = 0)
  expect_named(s_z, c("term", "estimate", "std.error", "z", "p.value"))

  s_null <- summary(p, null = NA)
  expect_false("p.value" %in% names(s_null))
})

test_that("projections of subgroup curves can be subset in summary()", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21, by = ~Male)
  p <- curve_projection(a, "linear")

  expect_length(coef(p), 4L)
  expect_identical(names(coef(p))[3L], "Male = 1 | (Intercept)")

  s <- summary(p)
  s1 <- summary(p, subset = Male == 1)

  expect_identical(nrow(s1), 2L)
  expect_equal(s1$estimate, s$estimate[3:4], tolerance = 1e-10)
  expect_equal(s1$std.error, s$std.error[3:4], tolerance = 1e-10)

  pc <- curve_projection(curve_contrast(a), "linear")
  expect_length(coef(pc), 2L)
  expect_output(print(summary(pc)), "Contrast Projection")
})

test_that("projections of a logistic ADRF use the link transformation", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p <- curve_projection(a, "linear")
  expect_true(all(get_est(p) > 0 & get_est(p) < 1))

  # The fitted curve is linear on the logit scale
  expect_equal(diff(diff(qlogis(get_est(p)))), rep(0, 19), tolerance = 1e-8)

  p_id <- curve_projection(a, "linear", transform = FALSE)
  expect_equal(diff(diff(get_est(p_id))), rep(0, 19), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(coef(p), coef(p_id))))
})

test_that("anova() keeps long formulas in a single heading", {
  fit <- lm(Math ~ poly(logBLL, 6) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p5 <- curve_projection(a, ~ logBLL + I(logBLL^2) + I(logBLL^3) + I(logBLL^4) + I(logBLL^5))
  p1 <- curve_projection(a, "linear")

  av <- anova(p5, p1)
  expect_length(attr(av, "heading"), 3L)

  out <- capture.output(print(av))
  expect_true(any(grepl("I(logBLL^5)", out, fixed = TRUE)))
  expect_true(any(grepl("Model 2: ~logBLL", out, fixed = TRUE)))
  expect_true(any(grepl("Extra terms", out, fixed = TRUE)))
})

test_that("anova() gives an informative error when the extra terms cannot be identified", {
  fit <- lm(Math ~ poly(logBLL, 4) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  # The ADRF is a quartic polynomial, so the fifth power cannot be identified
  p5 <- curve_projection(a, ~ logBLL + I(logBLL^2) + I(logBLL^3) + I(logBLL^4) + I(logBLL^5))
  p1 <- curve_projection(a, "linear")

  expect_error(anova(p5, p1), "singular")

  p4 <- curve_projection(a, ~ logBLL + I(logBLL^2) + I(logBLL^3) + I(logBLL^4))
  expect_no_error(anova(p4, p1))
})

test_that("anova() compares nested projections", {
  fit <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p1 <- curve_projection(a, "linear")
  p0 <- curve_projection(a, "flat")
  p3 <- curve_projection(a, "cubic")

  av <- anova(p1, p0)
  expect_s3_class(av, "anova.curve_projection")
  expect_identical(av$df, 1L)

  t_slope <- coef(p1)[["logBLL"]] / sqrt(vcov(p1)["logBLL", "logBLL"])
  expect_equal(av$F, t_slope^2, tolerance = 1e-8)
  expect_equal(av$p.value, pf(t_slope^2, 1, df.residual(fit), lower.tail = FALSE),
               tolerance = 1e-8)

  av3 <- anova(p3, p0, df = Inf)
  expect_identical(av3$df, 3L)
  expect_true("Chisq" %in% names(av3))

  expect_output(print(av), "Wald Test")
  expect_output(print(av3), "Chisq|χ")

  expect_error(anova(p0, p1), "nested")
  expect_error(anova(p1, curve_projection(a, ~ I(logBLL^2))), "nested")
  expect_error(anova(p1), "object2")
  expect_error(anova(p1, a), "curve_projection")

  a2 <- adrf(fit, treat = "logBLL", n = 11)
  expect_error(anova(p1, curve_projection(a2, "flat")), "same")
})

test_that("projections of curves without uncertainty can be summarized", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")

  p <- curve_projection(a, "linear")

  expect_null(vcov(p))
  expect_equal(coef(p)[["logBLL"]], coef(fit)[["logBLL"]], tolerance = 1e-8)

  s <- summary(p)
  expect_named(s, c("term", "estimate"))
  expect_equal(s$estimate, unname(coef(p)), tolerance = 1e-10)
  expect_output(print(s), "Projection Coefficients")
})

test_that("projection models producing non-finite values give an informative error", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  expect_error(suppressWarnings(curve_projection(a, ~ log(logBLL))), "non-finite")
  expect_error(suppressWarnings(summary(a, ~ log(logBLL))), "non-finite")
})

test_that("the covariance of a transformed projection matches the numerical Jacobian", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead,
             family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 21)

  p <- curve_projection(a, "linear", transform = TRUE)

  # Derivative of the projection coefficients with respect to the curve estimates
  coef_at <- function(y) {
    attr(a, ".est") <- y
    coef(curve_projection(a, "linear", transform = TRUE))
  }

  est <- attr(a, ".est")
  J <- vapply(seq_along(est), function(i) {
    h <- 1e-6
    e1 <- e0 <- est
    e1[i] <- e1[i] + h
    e0[i] <- e0[i] - h
    (coef_at(e1) - coef_at(e0)) / (2 * h)
  }, numeric(2L))

  expect_equal(unname(vcov(p)), unname(J %*% attr(a, ".vcov") %*% t(J)), tolerance = 1e-5)
})

test_that("polynomial projections work for treatments on a large scale", {
  d <- nhanes3lead
  d$bw <- 1000 * (d$logBLL + 1.5)
  fit <- lm(Math ~ poly(bw, 3) + Male + PIR, data = d)
  a <- adrf(fit, treat = "bw", data = d)

  p_big <- curve_projection(a, "quadratic")

  fit1 <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)
  p_unit <- curve_projection(adrf(fit1, treat = "logBLL"), "quadratic")

  # The fitted curves and their covariance do not depend on the scale
  expect_equal(get_est(p_big), get_est(p_unit), tolerance = 1e-8)
  expect_equal(get_vcov(p_big), get_vcov(p_unit), tolerance = 1e-6)
})

test_that("anova() requires the same transformation in both projections", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead,
             family = binomial)
  a <- adrf(fit, treat = "logBLL")

  p3 <- curve_projection(a, "cubic", transform = TRUE)
  p1 <- curve_projection(a, "linear", transform = TRUE)
  p1_id <- curve_projection(a, "linear", transform = FALSE)

  expect_true(is.finite(anova(p3, p1)$p.value))
  expect_error(anova(p3, p1_id), "same `transform`")
})

test_that("projections and anova() give informative errors for unusable models", {
  fit <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)

  a3 <- adrf(fit, treat = "logBLL", n = 3)
  expect_error(curve_projection(a3, "cubic"), "only 3 grid points")
  expect_error(summary(a3, "cubic"), "only 3 grid points")
  expect_error(summary(a3, ~ logBLL + I(logBLL^2) + I(logBLL^3)), "only 3 grid points")

  a <- adrf(fit, treat = "logBLL", n = 11)
  expect_error(curve_projection(a, ~ logBLL + I(2 * logBLL)), "linearly dependent")

  # Raw powers of a treatment on a large scale are not linearly dependent
  d <- transform(nhanes3lead, big = logBLL * 1e9)
  a_big <- adrf(lm(Math ~ poly(big, 3) + Male, data = d), treat = "big", n = 11)
  expect_no_error(curve_projection(a_big, "cubic"))

  a_none <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")
  expect_error(anova(curve_projection(a_none, "quadratic"), curve_projection(a_none, "linear")),
               "none")
})
