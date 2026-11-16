test_that("pointwise Wald intervals and tests match hand computations", {
  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)
  a <- reference_curve(adrf(fit, treat = "logBLL", n = 11), reference = 0)

  est <- a(logBLL = c(1, 2))
  s <- summary(est, simultaneous = FALSE)

  expect_s3_class(s, "summary.curve_est")
  expect_named(s, c("logBLL", "estimate", "std.error", "t", "p.value",
                    "conf.low", "conf.high"))

  se <- sqrt(diag(vcov(est)))
  df <- df.residual(fit)
  crit <- qt(.975, df)

  expect_equal(s$std.error, unname(se), tolerance = 1e-10)
  expect_equal(s$t, unname(coef(est) / se), tolerance = 1e-10)
  expect_equal(s$p.value, unname(2 * pt(-abs(coef(est) / se), df)), tolerance = 1e-10)
  expect_equal(s$conf.low, unname(coef(est) - crit * se), tolerance = 1e-10)
  expect_equal(s$conf.high, unname(coef(est) + crit * se), tolerance = 1e-10)
  expect_equal(c(attr(s, "crit")), crit, tolerance = 1e-10)
  expect_false(attr(s, "simultaneous"))
  expect_identical(attr(s, "null"), 0)

  # df = Inf gives z statistics
  s_z <- summary(est, simultaneous = FALSE, df = Inf)
  expect_true("z" %in% names(s_z))
  expect_equal(s_z$p.value, unname(2 * pnorm(-abs(coef(est) / se))), tolerance = 1e-10)

  expect_output(print(s), "ADRF Reference Estimates")
  expect_output(print(s), "pointwise")
  expect_output(print(s_z), "z\\* = 1.96")
})

test_that("simultaneous intervals are wider than pointwise intervals", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  est <- a(logBLL = c(0, 1, 2))

  s_p <- summary(est, simultaneous = FALSE)
  s_s <- summary(est, simultaneous = TRUE)

  expect_true(all(s_s$conf.high - s_s$conf.low > s_p$conf.high - s_p$conf.low))
  expect_gt(c(attr(s_s, "crit")), c(attr(s_p, "crit")))

  # Simultaneous p-values are at least as large as pointwise p-values
  s_p0 <- summary(est, simultaneous = FALSE, null = 8)
  s_s0 <- summary(est, simultaneous = TRUE, null = 8)
  expect_true(all(s_s0$p.value >= s_p0$p.value - 1e-8))

  # A single estimate cannot have a simultaneous interval
  s1 <- summary(a(logBLL = 1), simultaneous = TRUE)
  expect_false(attr(s1, "simultaneous"))
})

test_that("ADRFs omit tests by default and use null when supplied", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  s <- summary(a(logBLL = c(0, 1)))
  expect_false("p.value" %in% names(s))
  expect_null(attr(s, "null"))

  s_null <- summary(a(logBLL = c(0, 1)), null = 8)
  expect_true("p.value" %in% names(s_null))
  expect_output(print(s_null), "Null value: 8")

  s_noci <- summary(a(logBLL = c(0, 1)), conf_level = 0, null = 8)
  expect_false("conf.low" %in% names(s_noci))

  # No inference at all just returns the estimates
  s_none <- summary(a(logBLL = c(0, 1)), conf_level = 0, null = NA)
  expect_named(s_none, c("logBLL", "estimate"))

  expect_error(summary(a(logBLL = 1), conf_level = 1))
  expect_error(summary(a(logBLL = 1), conf_level = -.1))
  expect_error(summary(a(logBLL = 1), null = c(1, 2)))
  expect_error(summary(a(logBLL = 1), df = 0))
})

test_that("transformed intervals for a logistic ADRF stay within [0, 1]", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)

  est <- a(logBLL = c(0, 1, 2))

  s_t <- summary(est, simultaneous = FALSE)
  s_u <- summary(est, simultaneous = FALSE, transform = FALSE)

  # Transformed intervals are asymmetric and do not report standard errors
  expect_false("std.error" %in% names(s_t))
  expect_true("std.error" %in% names(s_u))
  expect_true(all(s_t$conf.low > 0 & s_t$conf.high < 1))
  expect_equal(s_t$estimate, s_u$estimate, tolerance = 1e-8)

  lo <- qlogis(s_t$estimate) - qlogis(s_t$conf.low)
  hi <- qlogis(s_t$conf.high) - qlogis(s_t$estimate)
  expect_equal(lo, hi, tolerance = 1e-6)

  # Delta method on the logit scale
  p <- coef(est)
  se_logit <- sqrt(diag(vcov(est))) / (p * (1 - p))
  expect_equal(s_t$conf.low, unname(plogis(qlogis(p) - qnorm(.975) * se_logit)),
               tolerance = 1e-8)

  # Null values are transformed too
  s_null <- summary(est, simultaneous = FALSE, null = .1)
  z <- (qlogis(p) - qlogis(.1)) / se_logit
  expect_equal(s_null$p.value, unname(2 * pnorm(-abs(z))), tolerance = 1e-8)

  expect_error(summary(est, transform = "logit"), "transform")
})

test_that("other specifications of the logit transformation are equivalent", {
  skip_if_not_installed("patrick")

  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 1, 2))

  s_t <- summary(est, simultaneous = FALSE)

  patrick::with_parameters_test_that(
    "transform = {label}",
    {
      s <- summary(est, simultaneous = FALSE, transform = transform)
      expect_equal(s$conf.low, s_t$conf.low, tolerance = 1e-8)
    },
    patrick::cases(
      list(label = "qlogis", transform = qlogis),
      list(label = "binomial()", transform = binomial()),
      list(label = "make.link(\"logit\")", transform = make.link("logit"))
    )
  )
})

test_that("transform works with other built-in functions", {
  skip_if_not_installed("patrick")

  fit <- glm(Block ~ logBLL + Male + PIR, data = nhanes3lead, family = poisson)
  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 2))

  s_def <- summary(est, simultaneous = FALSE)

  # Each of these is the log link of the Poisson model, specified another way
  patrick::with_parameters_test_that(
    "transform = {label}",
    {
      s <- summary(est, simultaneous = FALSE, transform = transform)
      expect_equal(s$conf.low, s_def$conf.low, tolerance = 1e-8)
    },
    patrick::cases(
      list(label = "log", transform = log),
      list(label = "poisson()", transform = poisson()),
      list(label = "make.link(\"log\")", transform = make.link("log"))
    )
  )

  patrick::with_parameters_test_that(
    "transform = {label}",
    {
      s <- summary(est, simultaneous = FALSE, transform = transform)
      expect_true(all(s$conf.low < s$estimate))
    },
    patrick::cases(
      list(label = "sqrt", transform = sqrt),
      list(label = "exp", transform = exp)
    )
  )
})

test_that("transform is ignored for curves that are not pure ADRFs", {
  fit <- glm(I(Block >= 14) ~ logBLL + Male + PIR, data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)

  r <- reference_curve(a, reference = 0)
  s <- summary(r(logBLL = c(1, 2)), simultaneous = FALSE)

  expect_true("std.error" %in% names(s))
  se <- sqrt(diag(vcov(r(logBLL = c(1, 2)))))
  expect_equal(s$conf.low, unname(coef(r(logBLL = c(1, 2))) - qnorm(.975) * se),
               tolerance = 1e-8)
})

test_that("estimates with zero variance are handled", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  r <- reference_curve(adrf(fit, treat = "logBLL", n = 11), reference = 0)

  s <- summary(r(logBLL = c(0, 1, 2)))

  expect_equal(s$estimate[1], 0, tolerance = 1e-10)
  expect_true(is.na(s$p.value[1]))
  expect_false(anyNA(s$p.value[-1]))
  expect_equal(s$conf.low[1], 0, tolerance = 1e-10)
})

test_that("summary() of a curve_est with vcov = 'none' only gives estimates", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")

  s <- summary(a(logBLL = c(0, 1)))
  expect_named(s, c("logBLL", "estimate"))
  expect_output(print(s), "Estimate")
})

test_that("subsetting subgroups changes simultaneous but not pointwise inference", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", by = ~Male, n = 11)

  vals <- c(.5, 1.5, 2.2)
  est_all <- a(logBLL = vals)
  est_1 <- a(logBLL = vals, subset = Male == 1)
  rows_1 <- est_all$Male == 1

  # Pointwise intervals and tests depend only on each estimate
  s_all <- summary(est_all, simultaneous = FALSE, null = 8)
  s_1 <- summary(est_1, simultaneous = FALSE, null = 8)
  expect_equal(s_1$p.value, s_all$p.value[rows_1], tolerance = 1e-10)
  expect_equal(s_1$conf.low, s_all$conf.low[rows_1], tolerance = 1e-10)

  # Simultaneous inference covers fewer estimates in the subset, so its
  # intervals are narrower and its p-values smaller
  s_all_sim <- summary(est_all, simultaneous = TRUE, null = 8)
  s_1_sim <- summary(est_1, simultaneous = TRUE, null = 8)
  expect_true(all((s_1_sim$conf.high - s_1_sim$conf.low) <
                    (s_all_sim$conf.high - s_all_sim$conf.low)[rows_1]))
  expect_true(all(s_1_sim$p.value <= s_all_sim$p.value[rows_1] + 1e-8))
})

test_that("a single estimate on a transformed scale can be summarized", {
  fit <- glm(I(Block >= 14) ~ logBLL + Male, data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 5)

  s1 <- summary(a(logBLL = 1))
  s3 <- summary(a(logBLL = c(0, 1, 2)), simultaneous = FALSE)

  # A single estimate has a pointwise interval, the same as in a larger summary
  expect_equal(s1$conf.low, s3$conf.low[2], tolerance = 1e-10)
  expect_equal(s1$conf.high, s3$conf.high[2], tolerance = 1e-10)
  expect_true(s1$conf.low > 0 && s1$conf.high < 1)

  s1_null <- summary(a(logBLL = 1), null = .1)
  expect_true(is.finite(s1_null$p.value))
})

test_that("vcov() of estimates without uncertainty is NULL", {
  fit <- lm(Math ~ logBLL + Male, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5, vcov = "none")

  expect_null(vcov(a(logBLL = c(0, 1))))
  expect_named(coef(a(logBLL = c(0, 1))), c("logBLL = 0", "logBLL = 1"))
})

test_that("a transformation that cannot be applied to every estimate gives an error", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  est <- a(logBLL = c(0, 1, 2))
  expect_true(all(coef(est) > 1))

  # The estimates are outside the domain of the logit and probit
  for (tr in list(qlogis, binomial("probit"))) {
    expect_error(summary(est, simultaneous = FALSE, transform = tr), "cannot be applied")
    expect_error(plot(a, simultaneous = FALSE, transform = tr), "cannot be applied")
    expect_error(summary(a, hypothesis = "flat", method = "satterthwaite", transform = tr),
                 "cannot be applied")
    expect_error(curve_projection(a, "linear", transform = tr), "cannot be applied")
  }
})

test_that("transform accepts only transformations with a known inverse and derivative", {
  fit <- glm(I(Block >= 14) ~ logBLL + Male + PIR, data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 1, 2))

  for (tr in list(function(p) qlogis(p), "logit",
                  list(transform = qlogis, inv_transform = plogis))) {
    expect_error(summary(est, transform = tr), "transform")
    expect_error(plot(a, transform = tr), "transform")
    expect_error(summary(a, transform = tr), "transform")
    expect_error(curve_projection(a, "linear", transform = tr), "transform")
  }
})

test_that("non-integer degrees of freedom work with simultaneous inference", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 1, 2))

  # mvtnorm requires integer degrees of freedom, so they are rounded down
  set.seed(1)
  s1 <- summary(est, df = 10.5, null = 0)
  set.seed(1)
  s2 <- summary(est, df = 10, null = 0)
  expect_equal(s1$conf.low, s2$conf.low, tolerance = 1e-10)
  expect_equal(s1$p.value, s2$p.value, tolerance = 1e-10)
})

test_that("estimates on a very small scale are printed with enough digits", {
  fit <- lm(I(Math * 1e-6) ~ logBLL + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 3)
  s <- summary(a(logBLL = get_values(a)), simultaneous = FALSE)

  out <- capture.output(print(s))
  expect_true(any(grepl("0.000008", out, fixed = TRUE)))
  expect_false(any(grepl(" 0 ", out, fixed = TRUE)))
})

test_that("a null value outside the domain of the transformation gives an error", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 1, 2))

  for (nl in c(0, 1, 1.5, -.2)) {
    expect_error(summary(est, null = nl, simultaneous = FALSE), "outside the domain")
    expect_error(summary(a, hypothesis = nl, method = "satterthwaite"), "hypothesis")

    # Without the transformation, the test is performed on the original scale
    s_u <- summary(est, null = nl, simultaneous = FALSE, transform = FALSE)
    expect_true(all(is.finite(s_u$p.value)))
    expect_no_error(summary(a, hypothesis = nl, method = "satterthwaite", transform = FALSE))
  }

  # A null value inside the domain is transformed
  s_in <- summary(est, null = .1, simultaneous = FALSE)
  expect_false("std.error" %in% names(s_in))
})
