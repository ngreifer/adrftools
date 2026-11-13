skip_on_cran()
skip_if_not_installed("rstanarm")

fit_stan <- function(formula, family = gaussian(), data = nhanes3lead) {
  suppressWarnings(
    rstanarm::stan_glm(formula, data = data, family = family,
                       chains = 1, iter = 500, refresh = 0, seed = 1)
  )
}

# Posterior draws of the average counterfactual prediction, computed directly
manual_adrf_draws <- function(fit, data, treat, values) {
  vapply(values, function(v) {
    data[[treat]] <- v
    rowMeans(rstanarm::posterior_epred(fit, newdata = data))
  }, numeric(nrow(as.matrix(fit))))
}

test_that("Bayesian ADRFs average the posterior predictions draw by draw", {
  fit <- fit_stan(Math ~ logBLL + I(logBLL^2) + Male + PIR)

  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_identical(as.character(attr(a, ".vcov_type")), "posterior")
  expect_identical(attr(a, ".df"), Inf)

  draws <- manual_adrf_draws(fit, nhanes3lead, "logBLL", get_values(a))

  expect_equal(unname(attr(a, ".draws")), unname(draws), tolerance = 1e-8)
  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)
  expect_equal(get_vcov(a), unname(cov(draws)), tolerance = 1e-8)

  a_none <- adrf(fit, treat = "logBLL", n = 5, vcov = "none")
  expect_null(attr(a_none, ".draws"))
  expect_equal(get_est(a_none), get_est(a), tolerance = 1e-8)

  expect_identical(as.character(attr(adrf(fit, treat = "logBLL", n = 3, vcov = "posterior"),
                                      ".vcov_type")), "posterior")
  expect_error(adrf(fit, treat = "logBLL", vcov = "HC0"), "Bayesian")
  expect_warning(adrf(fit, treat = "logBLL", n = 3, cluster = ~Race), "ignored")
})

test_that("Bayesian ADRFs support percentile and Wald inference", {
  fit <- fit_stan(I(Block >= 14) ~ logBLL + Male + PIR, family = binomial())

  a <- adrf(fit, treat = "logBLL", n = 11)
  est <- a(logBLL = c(0, 1, 2))

  s <- summary(est, simultaneous = FALSE)
  expect_identical(attr(s, "ci.type"), "perc")

  dr <- attr(est, ".draws")
  expect_equal(s$conf.low, unname(apply(dr, 2, quantile, .025)), tolerance = 1e-6)
  expect_equal(s$conf.high, unname(apply(dr, 2, quantile, .975)), tolerance = 1e-6)

  s_sim <- summary(est)
  expect_true(all(s_sim$conf.high - s_sim$conf.low > s$conf.high - s$conf.low))

  s_w <- summary(est, ci.type = "wald", simultaneous = FALSE)
  expect_true(all(s_w$conf.low > 0 & s_w$conf.high < 1))

  s_null <- summary(est, simultaneous = FALSE, null = .1)
  # The smallest alpha at which the type 7 percentile interval excludes .1
  p_manual <- vapply(seq_len(ncol(dr)), function(j) {
    alphas <- seq(0, 1, by = 1e-4)
    excl <- vapply(alphas, function(a) {
      q <- quantile(dr[, j], c(a / 2, 1 - a / 2), type = 7)
      .1 < q[1] || .1 > q[2]
    }, logical(1L))
    min(c(alphas[excl], 1))
  }, numeric(1L))
  # Equal up to the step of the grid of alphas
  expect_lte(max(abs(s_null$p.value - p_manual)), 1e-4 + 1e-10)

  expect_error(summary(est, ci.type = "bc"))
})

test_that("Bayesian ADRFs work with curve modifiers, tests, and plots", {
  skip_if_not_installed("CompQuadForm")

  fit <- fit_stan(Math ~ logBLL + I(logBLL^2) + Male + PIR)
  a <- adrf(fit, treat = "logBLL", n = 11)

  m <- amef(a)
  expect_identical(dim(attr(m, ".draws")), c(250L, 11L))

  r <- reference_curve(a, reference = 1)
  expect_equal(r(logBLL = 1)$estimate, 0, tolerance = 1e-10)

  s_flat <- summary(a)
  expect_lt(s_flat$p.value, .05)
  expect_true(is.finite(summary(a, "linear")$p.value))
  expect_true(is.finite(summary(a, hypothesis = 8)$p.value))

  pc <- summary(point_contrast(a(logBLL = c(0, 1, 2))), simultaneous = FALSE)
  expect_identical(nrow(pc), 3L)

  expect_s3_class(plot(a), "ggplot")
  expect_s3_class(plot(a, simultaneous = FALSE), "ggplot")
  expect_s3_class(plot(a, ci.type = "wald"), "ggplot")
  expect_s3_class(plot(amef(a), simultaneous = FALSE), "ggplot")
})

test_that("brms models are supported", {
  skip_if_not_installed("brms")
  skip_if_not(identical(Sys.getenv("ADRFTOOLS_TEST_BRMS"), "true"),
              "Set ADRFTOOLS_TEST_BRMS=true to compile a brms model")

  fit <- suppressMessages(suppressWarnings(
    brms::brm(Math ~ logBLL + Male + PIR, data = nhanes3lead,
              chains = 1, iter = 500, refresh = 0, seed = 1, silent = 2)
  ))

  a <- adrf(fit, treat = "logBLL", n = 3)

  draws <- vapply(get_values(a), function(v) {
    rowMeans(brms::posterior_epred(fit, newdata = transform(nhanes3lead, logBLL = v)))
  }, numeric(250L))

  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)
})

test_that("Bayesian subgroup ADRFs are ordered by subgroup", {
  fit <- fit_stan(Math ~ logBLL * Male + PIR)
  a <- adrf(fit, treat = "logBLL", n = 3, by = ~Male)

  vals <- get_values(a)
  draws <- cbind(manual_adrf_draws(fit, nhanes3lead[nhanes3lead$Male == 0, ], "logBLL", vals),
                 manual_adrf_draws(fit, nhanes3lead[nhanes3lead$Male == 1, ], "logBLL", vals))

  expect_equal(unname(attr(a, ".draws")), unname(draws), tolerance = 1e-8)
  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)
  expect_identical(a(logBLL = vals)$Male, rep(0:1, each = 3))
})

test_that("Bayesian ADRFs can be projected", {
  fit <- fit_stan(Math ~ logBLL + I(logBLL^2) + Male + PIR)
  a <- adrf(fit, treat = "logBLL", n = 11)

  p <- curve_projection(a, "linear")

  expect_identical(dim(attr(p, ".proj_draws")), c(250L, 2L))
  expect_equal(unname(vcov(p)), unname(cov(attr(p, ".proj_draws"))), tolerance = 1e-10)

  s <- summary(p)
  expect_identical(attr(s, "ci.type"), "perc")
  expect_false(anyNA(s$p.value))
  expect_s3_class(plot(a, proj = p), "ggplot")
})

test_that("curve contrasts and flat projections carry the posterior draws", {
  fit <- fit_stan(Math ~ logBLL * Male + PIR)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Male)
  d_a <- attr(a, ".draws")

  cc <- curve_contrast(a)
  expect_equal(unname(attr(cc, ".draws")), unname(d_a[, 6:10] - d_a[, 1:5]),
               tolerance = 1e-10)

  a1 <- adrf(fit, treat = "logBLL", n = 5)
  p <- curve_projection(a1, "flat")
  vals <- get_values(a1)
  w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
  w <- w / sum(w)

  expect_identical(dim(attr(p, ".proj_draws")), c(250L, 1L))
  expect_equal(drop(attr(p, ".proj_draws")), unname(drop(attr(a1, ".draws") %*% w)),
               tolerance = 1e-8)

  pc <- summary(point_contrast(a(logBLL = c(0, 2))), simultaneous = FALSE)
  expect_identical(nrow(pc), 2L)
})
