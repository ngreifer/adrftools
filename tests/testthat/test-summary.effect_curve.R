test_that("omnibus test of a flat ADRF rejects a strong dose-response", {
  skip_if_not_installed("CompQuadForm")

  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  s <- summary(a)

  expect_s3_class(s, "summary.effect_curve")
  expect_identical(attr(s, "hypothesis"), "flat")
  expect_identical(attr(s, "method"), "imhof")
  expect_lt(s$p.value, 1e-6)

  expect_output(print(s), "Omnibus Curve Test")
  expect_output(print(s), "ADRF is flat")
  expect_output(print(s), "Imhof")
})

test_that("an exactly linear ADRF is not rejected by a test of linearity", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  expect_equal(summary(a, "linear", method = "satterthwaite")$p.value, 1)
  expect_equal(summary(a, ~logBLL, method = "satterthwaite")$p.value, 1)
})

test_that("p-value approximation methods agree", {
  skip_if_not_installed("CompQuadForm")

  set.seed(1)
  d <- nhanes3lead
  d$y <- d$Math + rnorm(nrow(d), sd = 15)
  fit <- lm(y ~ poly(logBLL, 2) * Male + PIR, data = d)
  a <- adrf(fit, treat = "logBLL", n = 21)

  r <- reference_curve(a, reference = 1)

  p_imhof <- summary(r, method = "imhof")$p.value

  # The p-value should not be trivial for the comparison to be informative
  expect_gt(p_imhof, 1e-4)
  expect_lt(p_imhof, .9)

  skip_if_not_installed("patrick")

  # `nsim` is used only by the simulation method
  patrick::with_parameters_test_that(
    "method = {method}",
    {
      set.seed(2)
      s <- summary(r, method = method, nsim = 1e5)

      expect_equal(s$p.value, p_imhof, tolerance = tol)
      expect_output(print(s), printed)
    },
    .cases = data.frame(
      method = c("davies", "saddlepoint", "liu", "satterthwaite", "sim"),
      tol = c(1e-3, .05, .2, .2, .05),
      printed = c("Davies", "saddlepoint", "Liu", "Satterthwaite", "100,000")
    )
  )
})

test_that("the imhof p-value matches a direct computation", {
  skip_if_not_installed("CompQuadForm")

  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)
  r <- reference_curve(a, reference = get_values(a)[1])

  s <- summary(r, df = Inf)

  est <- get_est(r)
  V <- get_vcov(r)
  vals <- get_values(r)
  w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
  w <- w / sum(w)

  stat <- sum(w * est^2)
  ev <- eigen(sqrt(w) * t(sqrt(w) * V), symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > 1e-10 * max(ev)]

  p <- CompQuadForm::imhof(stat, ev, epsabs = 1e-8, epsrel = 1e-8)$Qq

  expect_equal(s$p.value, max(p, 1e-12), tolerance = 1e-4)
})

test_that("numeric and default hypotheses depend on the curve type", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  s_num <- summary(a, hypothesis = 8, method = "satterthwaite")
  expect_identical(attr(s_num, "hypothesis"), 8)
  expect_output(print(s_num), "ADRF is 8")

  s_amef <- summary(amef(a), method = "satterthwaite")
  expect_identical(attr(s_amef, "hypothesis"), 0)
  expect_output(print(s_amef), "AMEF is 0")

  s_ref <- summary(reference_curve(a, 1), method = "satterthwaite")
  expect_output(print(s_ref), "difference from reference")
})

test_that("the default method depends on installed packages", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  local_mocked_bindings(is_installed = function(pkg, ...) FALSE, .package = "rlang")
  expect_identical(attr(summary(a), "method"), "saddlepoint")
})

test_that("subgroup curves are tested separately and can be subset", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  s <- summary(a, "linear", method = "satterthwaite")
  expect_identical(nrow(s), 2L)
  expect_named(s, c("Male", "p.value"))

  s1 <- summary(a, "linear", method = "satterthwaite", subset = Male == 1)
  expect_identical(nrow(s1), 1L)
  expect_equal(s1$p.value, s$p.value[2], tolerance = 1e-10)

  sc <- summary(curve_contrast(a), method = "satterthwaite")
  expect_identical(nrow(sc), 1L)
  expect_output(print(sc), "ADRF contrast")
})

test_that("subset gives an error without by or for a contrast", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  expect_error(summary(a, "linear", subset = Male == 1), "by")
  expect_error(summary(curve_projection(a, "linear"), subset = Male == 1), "by")

  a_by <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)
  cc <- curve_contrast(a_by)

  expect_error(summary(cc, "linear", subset = Male == 1), "curve_contrast")
  expect_error(summary(curve_projection(cc, "linear"), subset = Male == 1), "curve_contrast")
})

test_that("transformed tests are used for logistic ADRFs", {
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)
  a <- adrf(fit, treat = "logBLL", n = 11)

  h <- mean(get_est(a))
  s_t <- summary(a, hypothesis = h, method = "satterthwaite")
  s_u <- summary(a, hypothesis = h, method = "satterthwaite", transform = FALSE)

  expect_true(attr(s_t, "transformed"))
  expect_false(attr(s_u, "transformed"))
  expect_false(isTRUE(all.equal(s_t$p.value, s_u$p.value)))
  expect_output(print(s_t), "transformed ADRF")

  s_lin <- summary(a, hypothesis = "linear", method = "satterthwaite")
  expect_output(print(s_lin), "\\(transformed\\)")
})

test_that("summary.effect_curve() validates its input", {
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  expect_error(summary(a, hypothesis = "wiggly"))
  expect_error(summary(a, hypothesis = Math ~ logBLL), "one-sided")
  expect_error(summary(a, hypothesis = ~ PIR), "Only the treatment")
  expect_error(summary(a, hypothesis = list(1)), "string, formula, or a number")
  expect_error(summary(a, method = "exact"))
  expect_error(summary(a, method = "sim", nsim = 10))
  expect_error(summary(a, df = -1))

  a_none <- adrf(fit, treat = "logBLL", n = 11, vcov = "none")
  expect_error(summary(a_none), "vcov = \"none\"")
})

test_that("omnibus tests with subset use the matching bootstrap columns", {
  skip_if_not_installed("fwb")

  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  set.seed(1)
  a <- adrf(fit, treat = "logBLL", by = ~Male, vcov = "fwb", n = 11,
            fwb.args = list(R = 99, verbose = FALSE))

  s <- summary(a, "linear", method = "satterthwaite")
  s1 <- summary(a, "linear", method = "satterthwaite", subset = Male == 1)
  expect_equal(s1$p.value, s$p.value[2], tolerance = 1e-10)

  s_num <- summary(a, 8, method = "satterthwaite")
  s1_num <- summary(a, 8, method = "satterthwaite", subset = Male == 1)
  expect_equal(s1_num$p.value, s_num$p.value[2], tolerance = 1e-10)
})

test_that("omnibus p-values do not depend on the scale of the outcome", {
  skip_if_not_installed("CompQuadForm")

  set.seed(1)
  d <- transform(nhanes3lead, y = Math + rnorm(nrow(nhanes3lead), sd = 15))

  p_at <- function(sc, ...) {
    d$ys <- d$y * sc
    a <- adrf(lm(ys ~ poly(logBLL, 2) * Male + PIR, data = d), treat = "logBLL", n = 21)
    c(flat = summary(a, "flat", ...)$p.value,
      linear = summary(a, "linear", ...)$p.value,
      num = summary(a, hypothesis = 8 * sc, ...)$p.value,
      ref = summary(reference_curve(a, 1), ...)$p.value)
  }

  for (mt in c("imhof", "davies", "saddlepoint", "satterthwaite")) {
    p1 <- p_at(1, method = mt)
    expect_equal(p_at(1e-8, method = mt), p1, tolerance = 1e-6, label = mt)
    expect_equal(p_at(1e8, method = mt), p1, tolerance = 1e-6, label = mt)
  }
})

test_that("omnibus tests of hypotheses that hold exactly give p = 1 at any scale", {
  for (sc in c(1e-8, 1, 1e8)) {
    d <- transform(nhanes3lead, ys = Math * sc)
    a <- adrf(lm(ys ~ logBLL + Male + PIR, data = d), treat = "logBLL", n = 21)

    expect_equal(summary(a, "linear")$p.value, 1, label = sprintf("scale %g", sc))
    expect_equal(summary(a, "quadratic")$p.value, 1, label = sprintf("scale %g", sc))
  }
})

test_that("omnibus p-values are valid probabilities for strong effects", {
  skip_if_not_installed("CompQuadForm")

  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR + Race, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 31)

  for (mt in c("imhof", "davies", "saddlepoint", "satterthwaite", "liu")) {
    p <- summary(a, "flat", method = mt)$p.value
    expect_true(p >= 0 && p <= 1, label = mt)
    expect_lt(p, 1e-4)
  }

  # Imhof's and Davies's methods agree with the saddlepoint approximation where
  # they hand off to it
  expect_equal(summary(a, "flat", method = "imhof")$p.value,
               summary(a, "flat", method = "saddlepoint")$p.value, tolerance = 1e-10)
})

test_that(".quad_form_tail() is accurate in the tails", {
  skip_if_not_installed("CompQuadForm")

  # With equal eigenvalues, the distribution is a scaled chi-squared or F
  # distribution, so exact tail probabilities are available
  ev <- rep(2, 4)

  for (target in c(.3, 1e-3, 1e-6, 1e-10)) {
    stat <- 2 * qchisq(target, 4, lower.tail = FALSE)
    stat_f <- 2 * 4 * qf(target, 4, 30, lower.tail = FALSE)

    for (mt in c("imhof", "davies", "saddlepoint")) {
      p <- adrftools:::.quad_form_tail(stat, ev, method = mt)
      p_f <- adrftools:::.quad_form_tail(stat_f, ev, method = mt, df = 30)

      expect_equal(p, target, tolerance = .15, label = sprintf("%s, df = Inf, p = %g", mt, target))
      expect_equal(p_f, target, tolerance = .15, label = sprintf("%s, df = 30, p = %g", mt, target))
    }
  }

  # The saddlepoint approximation needs no other packages and matches survey's
  skip_if_not_installed("survey")
  ev2 <- c(1, .5, .1)
  expect_equal(adrftools:::.quad_form_tail(8, ev2, method = "saddlepoint"),
               survey::pchisqsum(8, rep(1, 3), ev2, lower.tail = FALSE, method = "saddlepoint"),
               tolerance = 1e-6)
})

test_that(".omnibus_p_value() handles degenerate and out-of-range cases", {
  # No variance in the tested quantity
  expect_identical(adrftools:::.omnibus_p_value(1e-30, c(1e-25, 1e-26), scale = 1,
                                                method = "satterthwaite"), 1)
  expect_identical(adrftools:::.omnibus_p_value(2, c(1e-25, 1e-26), scale = 1,
                                                method = "satterthwaite"), 0)

  # A statistic of 0
  expect_identical(adrftools:::.omnibus_p_value(0, c(1, .5), scale = 1,
                                                method = "satterthwaite"), 1)

  # Results are restricted to [0, 1]
  local_mocked_bindings(.quad_form_tail = function(...) -.01)
  expect_identical(adrftools:::.omnibus_p_value(5, c(1, .5), scale = 1, method = "liu"), 0)
})
