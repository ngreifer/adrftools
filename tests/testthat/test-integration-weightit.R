skip_if_not_installed("WeightIt")

test_that("glm_weightit() with estimated weights gives correct point estimates", {
  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age + Race, data = nhanes3lead,
                          method = "glm")
  fit <- WeightIt::glm_weightit(Math ~ poly(logBLL, 2) * Male + PIR,
                                data = nhanes3lead, weightit = W)

  a <- adrf(fit, treat = "logBLL", n = 5)

  # The balancing weights are not target weights
  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))
})

test_that("glm_weightit() conditional variance uses the model's M-estimation covariance", {
  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = nhanes3lead, method = "glm")
  fit <- WeightIt::glm_weightit(Math ~ logBLL + Male + PIR, data = nhanes3lead,
                                weightit = W)

  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional")

  vals <- get_values(a)
  G <- cbind(1, vals, mean(nhanes3lead$Male), mean(nhanes3lead$PIR))
  expect_equal(get_vcov(a), unname(G %*% vcov(fit) %*% t(G)), tolerance = 1e-8)
})

test_that("glm_weightit() with vcov = 'HC0' treats the weights as fixed in the unconditional variance", {
  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = nhanes3lead, method = "glm")
  f <- Math ~ poly(logBLL, 2) * Male + PIR

  fit_hc0 <- WeightIt::glm_weightit(f, data = nhanes3lead, weightit = W, vcov = "HC0")
  fit_asympt <- WeightIt::glm_weightit(f, data = nhanes3lead, weightit = W, vcov = "asympt")

  d <- nhanes3lead
  d$w <- W$weights
  fit_lm <- lm(f, data = d, weights = w)

  a_hc0 <- adrf(fit_hc0, treat = "logBLL", n = 4)
  a_asympt <- adrf(fit_asympt, treat = "logBLL", n = 4)
  a_lm <- adrf(fit_lm, treat = "logBLL", n = 4, data = d)

  # With fixed weights, the model is a weighted regression
  expect_equal(get_vcov(a_hc0), get_vcov(a_lm), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(get_vcov(a_hc0), get_vcov(a_asympt))))
})

test_that("glm_weightit() uses sampling weights as target weights", {
  set.seed(1)
  d <- transform(nhanes3lead, sw = rexp(nrow(nhanes3lead)))

  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = d, method = "glm",
                          s.weights = "sw")
  fit <- WeightIt::glm_weightit(I(Block >= 14) ~ logBLL * Male + PIR,
                                data = d, weightit = W, family = binomial)

  a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2))

  expect_equal(get_est(a), manual_adrf(fit, d, "logBLL", get_values(a), w = d$sw),
               tolerance = 1e-8)
  expect_s3_class(summary(a(logBLL = c(0, 1))), "summary.curve_est")
})

test_that("lm_weightit() is supported", {
  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = nhanes3lead, method = "glm")
  fit <- WeightIt::lm_weightit(Math ~ logBLL * Male + PIR, data = nhanes3lead,
                               weightit = W)

  a <- adrf(fit, treat = "logBLL", n = 3)
  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
})

test_that("ordinal_weightit() models are supported on the link scale", {
  fit <- WeightIt::ordinal_weightit(ordered(cut(Math, 4)) ~ logBLL + PIR + Male,
                                    data = nhanes3lead)

  # WeightIt (> 2.1.0) has a hatvalues() method for these models, so HC2 is
  # used, scaling each unit's influence by 1 / sqrt(1 - h); without it, HC1 is
  # used with a message
  if (is_null(utils::getS3method("hatvalues", "ordinal_weightit", optional = TRUE))) {
    expect_message(a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), type = "link"),
                   "HC2")
  }
  else {
    expect_no_message(a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), type = "link"))

    h <- stats::hatvalues(fit)
    expect_equal(sum(h), ncol(sandwich::estfun(fit)), tolerance = 1e-8)
    expect_equal(adrftools:::.coef_influence(fit, hc2 = TRUE),
                 adrftools:::.coef_influence(fit, hc2 = FALSE) / sqrt(1 - h),
                 tolerance = 1e-10)
  }

  ref <- vapply(get_values(a), function(v) {
    nd <- transform(nhanes3lead, logBLL = v)
    mean(marginaleffects::get_predict(fit, newdata = nd, type = "link")$estimate)
  }, numeric(1L))

  expect_equal(get_est(a), ref, tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))

  # On the response scale, there is one prediction per category
  expect_error(suppressMessages(adrf(fit, treat = "logBLL", n = 3)),
               "each category of the outcome")
})

test_that("multinom_weightit() models give an error", {
  fit <- WeightIt::multinom_weightit(cut(Math, 3) ~ logBLL + PIR + Male,
                                     data = nhanes3lead)

  expect_error(suppressMessages(adrf(fit, treat = "logBLL", n = 3)),
               "each category of the outcome")
})
