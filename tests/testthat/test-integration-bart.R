skip_on_cran()

bart_data <- function() nhanes3lead[1:600, ]

test_that("bartisan ADRFs average posterior predictions draw by draw", {
  skip_if_not_installed("bartisan")

  d <- bart_data()
  set.seed(1)
  fit <- bartisan::bartisan(Math ~ logBLL + Male + PIR + Age, data = d,
                            family = gaussian(),
                            control = bartisan::bartisan_control(num_burn = 50,
                                                                 num_draws = 100,
                                                                 num_trees = 20))

  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_identical(as.character(attr(a, ".vcov_type")), "posterior")

  draws <- vapply(get_values(a), function(v) {
    p <- marginaleffects::get_predict(fit, newdata = transform(d, logBLL = v))
    colMeans(attr(p, "posterior_draws"))
  }, numeric(100L))

  expect_equal(unname(attr(a, ".draws")), unname(draws), tolerance = 1e-8)
  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)

  # Untransformed inference does not depend on the family
  s <- summary(a(logBLL = c(0, 1)), transform = FALSE, simultaneous = FALSE)
  expect_true(all(s$conf.low < s$estimate & s$estimate < s$conf.high))

  expect_s3_class(plot(a, transform = FALSE, simultaneous = FALSE), "ggplot")
  expect_s3_class(plot(amef(a)), "ggplot")
})

test_that("bartisan binomial ADRFs are probabilities", {
  skip_if_not_installed("bartisan")

  d <- bart_data()
  set.seed(1)
  fit <- bartisan::bartisan(I(Block >= 12) ~ logBLL + Male + PIR + Age, data = d,
                            family = binomial(),
                            control = bartisan::bartisan_control(num_burn = 50,
                                                                 num_draws = 100,
                                                                 num_trees = 20))

  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_true(all(attr(a, ".draws") > 0 & attr(a, ".draws") < 1))
  expect_identical(attr(a, ".family")$link, "logit")
})

test_that("dbarts ADRFs are supported", {
  skip_if_not_installed("dbarts")

  d <- bart_data()
  set.seed(1)
  fit <- dbarts::bart2(Math ~ logBLL + Male + PIR + Age, data = d,
                       keepTrees = TRUE, verbose = FALSE,
                       n.chains = 1, n.samples = 100, n.burn = 50)

  # dbarts models are not supported by insight, so data must be supplied
  a <- adrf(fit, treat = "logBLL", n = 5, data = d)

  expect_identical(as.character(attr(a, ".vcov_type")), "posterior")

  # predict.bart() returns one row per draw
  draws <- vapply(get_values(a), function(v) {
    rowMeans(predict(fit, newdata = transform(d, logBLL = v), type = "ev"))
  }, numeric(100L))

  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-6)

  s <- summary(a(logBLL = c(0, 1)), simultaneous = FALSE)
  expect_s3_class(s, "summary.curve_est")
  expect_s3_class(plot(a, simultaneous = FALSE), "ggplot")
})

test_that("stan4bart ADRFs are supported", {
  skip_if_not_installed("stan4bart")

  d <- bart_data()
  fit <- suppressWarnings(
    stan4bart::stan4bart(Math ~ bart(logBLL + PIR + Age) + Male, data = d,
                         chains = 1, iter = 200, warmup = 100, verbose = -1,
                         seed = 1, bart_args = list(keepTrees = TRUE))
  )

  a <- adrf(fit, treat = "logBLL", n = 3, data = d)

  draws <- vapply(get_values(a), function(v) {
    colMeans(predict(fit, newdata = transform(d, logBLL = v), type = "ev"))
  }, numeric(100L))

  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-6)
  expect_s3_class(summary(a(logBLL = c(0, 1)), simultaneous = FALSE), "summary.curve_est")
})

test_that("bartisan families support transformed inference", {
  skip_if_not_installed("bartisan")

  d <- bart_data()
  set.seed(1)
  fit <- bartisan::bartisan(Math ~ logBLL + Male + PIR, data = d, family = gaussian(),
                            control = bartisan::bartisan_control(num_burn = 50,
                                                                 num_draws = 100,
                                                                 num_trees = 20))
  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_identical(attr(a, ".family")$link, "identity")
  expect_s3_class(summary(a(logBLL = c(0, 1))), "summary.curve_est")
  expect_s3_class(plot(a), "ggplot")
})
