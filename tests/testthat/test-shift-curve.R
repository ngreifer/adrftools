# Weighted mean predictions with every unit's treatment shifted, computed
# directly with predict(), as an independent check on shift_curve()
manual_shift <- function(fit, data, treat, shifts, shift_fun, w = NULL, type = "response") {
  if (is.null(w)) {
    w <- rep(1, nrow(data))
  }

  vapply(shifts, function(v) {
    nd <- data
    nd[[treat]] <- shift_fun(data[[treat]], v)
    stats::weighted.mean(stats::predict(fit, newdata = nd, type = type), w)
  }, numeric(1L))
}

d_pos <- nhanes3lead
d_pos$BLL <- exp(d_pos$logBLL)

test_that("shift curve estimates are mean predictions at the shifted treatment", {
  skip_if_not_installed("patrick")

  fit <- glm(Block >= 12 ~ poly(BLL, 2) * Male + PIR, data = d_pos, family = binomial)

  patrick::with_parameters_test_that(
    "shift_fun = {label}",
    {
      sc <- shift_curve(fit, treat = "BLL", shift = shift, shift_fun = shift_fun,
                        n = 5, wts = d_pos$Age, vcov = "none")
      shifts <- get_values(sc)

      expect_equal(get_est(sc),
                   manual_shift(fit, d_pos, "BLL", shifts, f, w = d_pos$Age),
                   tolerance = 1e-10)
    },
    label = c("additive", "multiplicative", "function"),
    shift = list(c(-1, 1), .2, c(1, 3)),
    shift_fun = list("additive", "multiplicative", function(a, delta) pmin(a, delta)),
    f = list(function(a, v) a + v, function(a, v) a * (1 + v), function(a, v) pmin(a, v))
  )
})

test_that("shift curves within subgroups use the units in each subgroup", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)

  sc <- shift_curve(fit, treat = "logBLL", shift = c(-.5, .5), n = 3, by = ~Male)

  for (m in 0:1) {
    d_m <- nhanes3lead[nhanes3lead$Male == m, ]
    expect_equal(get_est(sc)[3 * m + 1:3],
                 manual_shift(fit, d_m, "logBLL", c(-.5, 0, .5), function(a, v) a + v),
                 tolerance = 1e-10)
  }
})

test_that("a model with only the treatment is not collapsed to one unit", {
  fit <- lm(Math ~ poly(logBLL, 2), data = nhanes3lead)

  sc <- shift_curve(fit, treat = "logBLL", shift = c(-.5, .5), n = 3)

  expect_equal(get_est(sc),
               manual_shift(fit, nhanes3lead, "logBLL", c(-.5, 0, .5), function(a, v) a + v),
               tolerance = 1e-10)

  # The covariate-sampling term is not 0, since predictions differ by unit
  V <- manual_uncond_shift_vcov_lm(fit, nhanes3lead, "logBLL", c(-.5, 0, .5))
  expect_equal(get_vcov(sc), V, tolerance = 1e-8)
})

test_that("the unconditional variance matches a hand-computed influence function", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR + Age, data = nhanes3lead)

  sc <- shift_curve(fit, treat = "logBLL", shift = c(-.4, .4), n = 3)

  expect_equal(get_vcov(sc),
               manual_uncond_shift_vcov_lm(fit, nhanes3lead, "logBLL", c(-.4, 0, .4)),
               tolerance = 1e-8)
})

test_that("the conditional variance matches marginaleffects", {
  fit <- glm(Block >= 12 ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead,
             family = binomial)

  sc <- shift_curve(fit, treat = "logBLL", shift = c(-.4, .4), n = 3, vcov = "HC0")

  me_se <- vapply(c(-.4, 0, .4), function(v) {
    nd <- nhanes3lead
    nd$logBLL <- nd$logBLL + v
    marginaleffects::avg_predictions(fit, newdata = nd, vcov = "HC0")$std.error
  }, numeric(1L))

  expect_equal(get_se(sc), me_se, tolerance = 1e-5)
})

test_that("the default shift is half a weighted standard deviation on either side of 0", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  w <- nhanes3lead$Age

  sc <- shift_curve(fit, treat = "logBLL", wts = w, by = ~Male)

  # The weighted SD of cobalt::col_w_sd(), which does not depend on the scale
  # of the weights
  wn <- w / sum(w)
  A <- nhanes3lead$logBLL
  sd_w <- sqrt(sum(wn * (A - sum(wn * A))^2) / (1 - sum(wn^2)))
  v <- get_values(sc)

  expect_length(v, 51L)
  expect_equal(range(v), c(-.5, .5) * sd_w, tolerance = 1e-10)
  expect_identical(v[26L], 0)
  expect_equal(attr(sc, ".shift")[["sd_units"]], c(-.5, .5))

  # Rescaling the weights, including to sum to 1, does not change the grid
  for (k in c(1 / sum(w), 1e3)) {
    expect_equal(get_values(shift_curve(fit, treat = "logBLL", wts = w * k, by = ~Male)), v,
                 tolerance = 1e-10)
  }

  # A single number is in standard deviations for an additive shift and in the
  # units of the shift otherwise
  expect_equal(range(get_values(shift_curve(fit, treat = "logBLL", shift = 1, n = 3))),
               c(-1, 1) * sd(nhanes3lead$logBLL), tolerance = 1e-10)

  fit_pos <- lm(Math ~ poly(BLL, 2) * Male + PIR, data = d_pos)
  sc_m <- shift_curve(fit_pos, treat = "BLL", shift = .1, shift_fun = "multiplicative", n = 3)
  expect_equal(get_values(sc_m), c(-.1, 0, .1))
  expect_null(attr(sc_m, ".shift")[["sd_units"]])
})

test_that("the ADE matches marginaleffects::avg_slopes()", {
  fit <- glm(Block >= 12 ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead,
             family = binomial)

  a <- ade(shift_curve(fit, treat = "logBLL", vcov = "HC0"))
  s <- marginaleffects::avg_slopes(fit, variables = "logBLL", vcov = "HC0")

  expect_s3_class(a, "ade_est")
  expect_equal(a$estimate, s$estimate, tolerance = 1e-5)
  expect_equal(get_se(a), s$std.error, tolerance = 1e-4)
})

test_that("the ADE is the derivative of the shift curve at 0", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  sc <- shift_curve(fit, treat = "logBLL")

  a <- ade(sc)
  d0 <- amef(sc)(shift = 0)

  expect_equal(a$estimate, d0$estimate)
  expect_equal(get_vcov(a), get_vcov(d0))
  expect_named(coef(a), "ADE")

  # The exact average derivative, by central differences of the predictions
  eps <- 1e-6
  dp <- dm <- nhanes3lead
  dp$logBLL <- dp$logBLL + eps
  dm$logBLL <- dm$logBLL - eps
  exact <- mean(predict(fit, dp) - predict(fit, dm)) / (2 * eps)

  expect_equal(a$estimate, exact, tolerance = 1e-6)

  # The ADE does not depend on the range of the shift curve
  expect_equal(ade(shift_curve(fit, treat = "logBLL", shift = c(-.1, .3), n = 41))$estimate,
               exact, tolerance = 1e-6)
})

test_that("the ADE of a multiplicative shift weights the slopes by the treatment", {
  fit <- lm(Math ~ poly(BLL, 3) * Male + PIR, data = d_pos)

  a <- ade(shift_curve(fit, treat = "BLL", shift = .1, shift_fun = "multiplicative"))

  eps <- 1e-6
  dp <- dm <- d_pos
  dp$BLL <- dp$BLL + eps
  dm$BLL <- dm$BLL - eps
  slopes <- (predict(fit, dp) - predict(fit, dm)) / (2 * eps)

  expect_equal(a$estimate, mean(d_pos$BLL * slopes), tolerance = 1e-6)
})

test_that("ade() works within subgroups and on curve contrasts", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  sc <- shift_curve(fit, treat = "logBLL", by = ~Male)

  a <- ade(sc)
  expect_equal(nrow(a), 2L)

  a1 <- ade(sc, subset = Male == 1)
  expect_equal(a1$estimate, a$estimate[2L])

  ac <- ade(curve_contrast(sc))
  L <- matrix(c(-1, 1), nrow = 1L)
  expect_equal(abs(ac$estimate), abs(diff(a$estimate)), tolerance = 1e-10)
  expect_equal(get_vcov(ac), L %*% get_vcov(a) %*% t(L), tolerance = 1e-10)

  # A reference curve has the same derivative
  expect_equal(ade(reference_curve(sc, 0))$estimate, a$estimate, tolerance = 1e-8)
})

test_that("shift curves work with the other effect curve functions", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  sc <- shift_curve(fit, treat = "logBLL")

  expect_s3_class(sc, c("shift_curve", "adrf_curve", "effect_curve"))
  expect_named(sc(shift = 0), c("shift", "estimate"))

  r <- reference_curve(sc, 0)
  expect_equal(r(shift = .2)$estimate, sc(shift = .2)$estimate - sc(shift = 0)$estimate,
               tolerance = 1e-10)

  expect_s3_class(summary(sc, hypothesis = "flat"), "summary.effect_curve")
  expect_length(coef(curve_projection(sc, "linear")), 2L)

  expect_output(print(sc), "shift curve")
  expect_output(print(sc), "SD")
  expect_output(print(summary(ade(sc))), "Average Derivative Effect")

  p <- plot(sc)
  expect_identical(p$labels$x, "Shift in logBLL")
})

test_that("shift_curve() and ade() validate their input", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  sc <- shift_curve(fit, treat = "logBLL", n = 11)

  expect_error(shift_curve(fit, treat = "logBLL", shift_fun = "multiplicative"),
               "nonnegative")
  expect_error(shift_curve(fit, treat = "logBLL", shift_fun = function(a, d) a + d),
               "must be supplied")
  expect_error(shift_curve(fit, treat = "logBLL", shift_fun = "rotate"))
  expect_error(shift_curve(fit, treat = "logBLL", shift_fun = 3), "must be")
  expect_error(shift_curve(fit, treat = "logBLL", shift_fun = function(a) a),
               "two arguments")
  expect_error(shift_curve(fit, treat = "logBLL", shift = c(-1, 1),
                           shift_fun = function(a, d) a[-1L] + d),
               "one value for each")
  expect_error(suppressWarnings(shift_curve(fit, treat = "logBLL", shift = c(-1, 1),
                                            shift_fun = function(a, d) log(a + d))),
               "non-finite")
  expect_error(shift_curve(fit, treat = "logBLL", shift = c(1, 1)), "must differ")
  expect_error(shift_curve(fit, treat = "logBLL", shift = -1))

  fit_log <- lm(Math ~ log(BLL) * Male + PIR, data = d_pos)
  expect_error(suppressWarnings(shift_curve(fit_log, treat = "BLL", shift = c(-5, 0))),
               "narrower")

  expect_error(summary(sc, integrate = "empirical"), "shift curve")

  expect_error(ade(adrf(fit, treat = "logBLL", n = 11)), "shift curve")
  expect_error(ade(shift_curve(fit, treat = "logBLL", shift = c(.1, .5), n = 11)),
               "includes 0")
  expect_error(ade(shift_curve(fit, treat = "logBLL", shift = c(-1, 1), n = 11,
                               shift_fun = function(a, d) pmin(a, d))),
               "unchanged")
  expect_error(ade(amef(sc)), "amef_curve")
  expect_error(ade(sc, subset = Male == 1), "by")

  # Multiplicative shifts below -1 would make the treatment negative
  fit_pos <- lm(Math ~ poly(BLL, 2) * Male + PIR, data = d_pos)
  expect_error(shift_curve(fit_pos, treat = "BLL", shift = 2, shift_fun = "multiplicative"),
               "below -1")
  expect_no_error(shift_curve(fit_pos, treat = "BLL", shift = c(-1, 0),
                              shift_fun = "multiplicative", n = 3))

  # Arguments of the other function
  expect_error(shift_curve(fit, treat = "logBLL", range = c(-1, 1)), "shift")
  expect_error(adrf(fit, treat = "logBLL", shift_fun = "additive"), "shift_curve")
  expect_error(adrf(fit, treat = "logBLL", shift = 1), "shift_curve")

  # Messages refer to the shift, not the treatment range
  expect_error(sc(shift = 2), "shift_curve")
  expect_error(reference_curve(sc, 2), "shift_curve")
  expect_error(curve_contrast(sc), "shift_curve")

  sc_by <- shift_curve(fit, treat = "logBLL", by = ~Male, n = 11)
  expect_error(point_contrast(ade(sc_by)), "curve_contrast")
})

test_that("a primitive function can be supplied to shift_fun", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  expect_equal(shift_curve(fit, treat = "logBLL", shift = c(-.5, .5), n = 5,
                           shift_fun = `+`)(shift = c(-.5, .5))$estimate,
               shift_curve(fit, treat = "logBLL", shift = c(-.5, .5), n = 5)(shift = c(-.5, .5))$estimate)
})

test_that("a projection must come from a shift curve with the same shift", {
  fit_pos <- lm(Math ~ poly(BLL, 2) * Male + PIR, data = d_pos)
  sa <- shift_curve(fit_pos, treat = "BLL", shift = c(-.2, .2), n = 11)
  sm <- shift_curve(fit_pos, treat = "BLL", shift = c(-.2, .2), n = 11,
                    shift_fun = "multiplicative")

  expect_error(plot(sa, proj = curve_projection(sm, "linear")), "proj")
  expect_error(anova(curve_projection(sa, "flat"), curve_projection(sm, "linear")))
})

test_that("simultaneous bands work for highly correlated estimates", {
  # The AMEF of a cubic polynomial over a narrow range has a correlation matrix
  # of rank 3
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  d <- amef(shift_curve(fit, treat = "logBLL"))
  v <- get_values(d)

  s <- summary(d(shift = v[c(5, 25, 45)]))
  expect_true(all(is.finite(s$conf.low)))
  expect_true(all(is.finite(s$p.value)))
  expect_s3_class(plot(d), "ggplot")

  # A correlation matrix of rank 1 gives the univariate critical value
  expect_equal(adrftools:::.mvt_crit(.95, adrftools:::.to_psd(matrix(1, 5, 5)), df = 30),
               qt(.975, 30), tolerance = 1e-3)
})

test_that("svyglm models with data-dependent bases match the equivalent glm", {
  skip_if_not_installed("survey")

  des <- survey::svydesign(~1, weights = ~MEC_wt, data = nhanes3lead)
  f <- Math ~ splines::ns(logBLL, 4) * (Male + Age + PIR)
  fit_s <- survey::svyglm(f, design = des)
  fit_g <- glm(f, data = nhanes3lead, weights = MEC_wt)

  expect_equal(adrf(fit_s, treat = "logBLL", n = 11)(logBLL = c(0, 1, 2))$estimate,
               adrf(fit_g, treat = "logBLL", n = 11, wts = TRUE)(logBLL = c(0, 1, 2))$estimate,
               tolerance = 1e-8)
  expect_equal(shift_curve(fit_s, treat = "logBLL", n = 11)(shift = c(-.3, 0, .3))$estimate,
               shift_curve(fit_g, treat = "logBLL", n = 11, wts = TRUE)(shift = c(-.3, 0, .3))$estimate,
               tolerance = 1e-8)
})

test_that("shift curves are supported with bootstrapping", {
  skip_if_not_installed("fwb")

  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  set.seed(123)
  sc <- shift_curve(fit, treat = "logBLL", n = 11, vcov = "fwb",
                    fwb.args = list(R = 20, verbose = FALSE))

  expect_equal(get_est(sc),
               manual_shift(fit, nhanes3lead, "logBLL", get_values(sc),
                            function(a, v) a + v),
               tolerance = 1e-10)
  expect_identical(dim(attr(sc, ".boot")[["t"]]), c(20L, 11L))

  a <- ade(sc)
  expect_true(is.finite(get_se(a)))
})

test_that("shift curves are supported with multiply imputed data", {
  skip_if_not_installed("mice")

  d <- nhanes3lead[, c("Math", "logBLL", "PIR", "Male")]
  set.seed(11)
  d$PIR[sample(nrow(d), 300)] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 2, printFlag = FALSE, seed = 12))
  fits <- with(imp, lm(Math ~ poly(logBLL, 2) * Male + PIR))

  sc <- shift_curve(fits, treat = "logBLL", shift = c(-.5, .5), n = 3)

  per_imp <- lapply(1:2, function(i) {
    di <- mice::complete(imp, i)
    shift_curve(lm(Math ~ poly(logBLL, 2) * Male + PIR, data = di), treat = "logBLL",
                shift = c(-.5, .5), n = 3, data = di)
  })

  expect_equal(get_est(sc), colMeans(do.call("rbind", lapply(per_imp, get_est))),
               tolerance = 1e-10)
})
