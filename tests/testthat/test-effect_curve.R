test_that("calling an effect_curve at grid points returns the grid estimates", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11)

  vals <- get_values(a)
  est <- a(logBLL = vals[c(2, 5, 9)])

  expect_s3_class(est, "curve_est")
  expect_s3_class(est, "data.frame")
  expect_named(est, c("logBLL", "estimate"))
  expect_equal(est$estimate, get_est(a)[c(2, 5, 9)], tolerance = 1e-10)
  expect_equal(unname(vcov(est)), get_vcov(a)[c(2, 5, 9), c(2, 5, 9)],
               tolerance = 1e-10)

  # The treatment argument is named after the treatment and can be positional
  expect_identical(names(formals(a))[1L], "logBLL")
  expect_equal(a(vals[2])$estimate, get_est(a)[2], tolerance = 1e-10)
})

test_that("interpolation between grid points is exact for a cubic ADRF", {
  # Without treatment-covariate interactions, the ADRF of a cubic model is a
  # cubic polynomial, which local cubic regression reproduces exactly
  fit <- lm(Math ~ poly(logBLL, 3, raw = TRUE) + Male + PIR, data = nhanes3lead)

  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "n = {n}",
    {
      a <- adrf(fit, treat = "logBLL", n = n, range = c(0, 2))
      new_vals <- c(.123, .777, 1.5, 1.91)

      expect_equal(a(logBLL = new_vals)$estimate,
                   manual_adrf(fit, nhanes3lead, "logBLL", new_vals),
                   tolerance = 1e-8)
    },
    n = c(4, 11, 51)
  )
})

test_that("interpolated estimates are linear combinations of grid estimates", {
  fit <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  est <- a(logBLL = c(.1, .2, .3))

  V <- vcov(est)
  expect_identical(dim(V), c(3L, 3L))
  expect_true(isSymmetric(unname(V)))
  expect_true(all(eigen(V, only.values = TRUE)$values > 0))
})

test_that("coef() and vcov() methods label the estimates", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5)

  est <- a(logBLL = c(0, 1))

  expect_identical(names(coef(est)), c("logBLL = 0", "logBLL = 1"))
  expect_identical(dimnames(vcov(est)), list(names(coef(est)), names(coef(est))))
})

test_that("effect curves estimated within subgroups can be subset when called", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Male)

  vals <- get_values(a)[c(1, 3)]

  est <- a(logBLL = vals)
  expect_named(est, c("Male", "logBLL", "estimate"))
  expect_identical(est$Male, c(0L, 0L, 1L, 1L))
  expect_equal(est$estimate, get_est(a)[c(1, 3, 6, 8)], tolerance = 1e-10)
  expect_identical(names(coef(est))[3L], sprintf("Male = 1 | logBLL = %s", vals[1]))

  est1 <- a(logBLL = vals, subset = Male == 1)
  expect_identical(est1$Male, c(1L, 1L))
  expect_equal(est1$estimate, get_est(a)[c(6, 8)], tolerance = 1e-10)
  expect_equal(unname(vcov(est1)), get_vcov(a)[c(6, 8), c(6, 8)], tolerance = 1e-10)

  expect_error(a(logBLL = 1, subset = Male == 3), "real subset")
  expect_error(a(logBLL = 1, subset = PIR > 1), "present in the original")
  expect_error(a(logBLL = 1, subset = 1), "logical")
})

test_that("subset is not an argument when by was not supplied", {
  fit <- lm(Math ~ logBLL + Male, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_identical(names(formals(a)), "logBLL")
  expect_error(a(logBLL = 1, subset = Male == 1), "unused argument")
})

test_that("calling an effect_curve validates the treatment values", {
  fit <- lm(Math ~ logBLL + Male, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_error(a(logBLL = "a"), "logBLL")
  expect_error(a(), "logBLL")
})

test_that("print methods produce output", {
  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Male)

  expect_output(print(a), "<effect_curve>")
  expect_output(print(a), "curve type: ADRF")
  expect_output(print(a), "response: Math")
  expect_output(print(a), "treatment: logBLL")
  expect_output(print(a), "by: Male")
  expect_output(print(a), "inference: unconditional")

  a_cl <- adrf(fit, treat = "logBLL", n = 5, cluster = ~Race)
  expect_output(print(a_cl), "unconditional \\(clustered, CR2\\)")

  expect_output(print(a(logBLL = c(0, 1))), "ADRF Estimates")

  # Printing a curve_est returns it invisibly
  est <- a(logBLL = 1)
  expect_output(expect_invisible(print(est)))
})

test_that("interpolation is exact for a linear ADRF with very few grid points", {
  skip_if_not_installed("patrick")

  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  # With fewer than four grid points, the local polynomial's degree is reduced
  # to one less than the number of points used
  patrick::with_parameters_test_that(
    "n = {n}",
    {
      a <- adrf(fit, treat = "logBLL", n = n, range = c(0, 2))
      expect_equal(a(logBLL = c(.3, 1.7))$estimate,
                   manual_adrf(fit, nhanes3lead, "logBLL", c(.3, 1.7)),
                   tolerance = 1e-8)
    },
    n = 2:3
  )
})

test_that("interpolation and the AMEF do not depend on the scale of the treatment", {
  est_at <- function(scale) {
    d <- nhanes3lead
    d$tr <- d$logBLL * scale
    fit <- lm(Math ~ poly(tr, 5) + Male + PIR, data = d)
    a <- adrf(fit, treat = "tr", data = d)
    v <- get_values(a)

    list(mid = coef(a(tr = (v[5] + v[6]) / 2)),
         amef = get_est(amef(a))[c(1, 10, 25, 51)] * scale)
  }

  ref <- est_at(1)

  for (scale in c(1e-9, 1e9)) {
    res <- est_at(scale)
    expect_equal(unname(res$mid), unname(ref$mid), tolerance = 1e-8)
    # The AMEF uses numerical differences, so allow for rounding error
    expect_equal(res$amef, ref$amef, tolerance = 1e-5)
  }
})

test_that("values just outside the grid are moved to its ends and others give an error", {
  fit <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 51)
  v <- get_values(a)
  w <- diff(range(v))

  # Within 0.1% of the width, the estimate at the end point, labeled with the
  # requested value
  e <- a(logBLL = c(min(v) - .0009 * w, max(v) + .0009 * w))
  e0 <- a(logBLL = range(v))
  expect_equal(e$estimate, e0$estimate, tolerance = 1e-12)
  expect_equal(vcov(e), vcov(e0), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(e$logBLL, c(min(v) - .0009 * w, max(v) + .0009 * w))

  # Farther out, which used to give NA or unstable extrapolation
  expect_error(a(logBLL = max(v) + .005 * w), "within the range")
  expect_error(a(logBLL = c(1, min(v) - .1 * w)), "within the range")
  expect_error(amef(a)(logBLL = max(v) + .005 * w), "within the range")

  # reference_curve() follows the same rule
  r <- reference_curve(a, reference = max(v) + .0009 * w)
  expect_equal(attr(r, ".reference"), max(v))
  expect_false(anyNA(get_est(r)))
  expect_error(reference_curve(a, reference = max(v) + .005 * w), "within the range")
})

test_that("print() reports the grid, subgroups, and sample size", {
  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 7, by = ~Race, subset = Age > 8)
  n <- sum(nhanes3lead$Age > 8)

  expect_output(print(a), "grid points: 7")
  expect_output(print(a), sprintf("by: Race \\(%s subgroups\\)",
                                  length(unique(nhanes3lead$Race))))
  expect_output(print(a), sprintf("sample size: %s",
                                  format(n, big.mark = ",")))
  expect_no_match(capture.output(print(a)), "imputed datasets")

  # Curves derived from the ADRF keep the sample size and inference details
  expect_output(print(amef(a)), sprintf("sample size: %s", format(n, big.mark = ",")))
  expect_output(print(reference_curve(a, 1)), "inference: unconditional \\(HC2\\)")
  expect_output(print(curve_contrast(a)), "inference: unconditional \\(HC2\\)")

  # A shift curve reports the number of shifts
  expect_output(print(shift_curve(fit, treat = "logBLL", n = 9)), "grid points: 9")
})

test_that("print() reports the small-sample adjustment or covariance matrix used", {
  skip_if_not_installed("patrick")

  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)

  patrick::with_parameters_test_that(
    "{label}",
    {
      a <- do.call("adrf", c(list(fit, treat = "logBLL", n = 3), args))
      expect_output(print(a), pattern)
    },
    patrick::cases(
      list(label = "unconditional", args = list(),
           pattern = "inference: unconditional \\(HC2\\)"),
      list(label = "unconditional with clustering", args = list(cluster = ~Race),
           pattern = "inference: unconditional \\(clustered, CR2\\)"),
      list(label = "conditional with the model's covariance", args = list(vcov = "conditional"),
           pattern = "inference: conditional\n"),
      list(label = "conditional with a named type", args = list(vcov = "HC3"),
           pattern = "inference: conditional \\(HC3\\)"),
      list(label = "conditional with a matrix", args = list(vcov = vcov(fit)),
           pattern = "inference: conditional \\(user-supplied\\)"),
      list(label = "conditional with clustering", args = list(vcov = "HC1", cluster = ~Race),
           pattern = "inference: conditional \\(clustered, HC1\\)"),
      list(label = "no inference", args = list(vcov = "none"),
           pattern = "inference: none\n")
    )
  )
})

test_that("print() reports HC1 and CR1 when HC2 and CR2 are not available", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("betareg")

  # bam() fits have no leverages
  fit_bam <- mgcv::bam(Math ~ s(logBLL) + Male, data = nhanes3lead)
  a_bam <- suppressMessages(adrf(fit_bam, treat = "logBLL", n = 3))
  expect_output(print(a_bam), "inference: unconditional \\(HC1\\)")

  # CR2 is available only for lm() and glm() models
  d <- nhanes3lead
  d$y <- (d$Math + 1) / (max(d$Math) + 2)
  fit_beta <- betareg::betareg(y ~ logBLL + Male, data = d)
  a_beta <- suppressMessages(adrf(fit_beta, treat = "logBLL", n = 3, cluster = ~Race))
  expect_output(print(a_beta), "inference: unconditional \\(clustered, CR1\\)")
})
