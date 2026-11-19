skip_if_not_installed("survey")

svy_data <- function() {
  transform(nhanes3lead, sw = MEC_wt / mean(MEC_wt))
}

test_that("svyglm ADRFs use the survey weights as target weights", {
  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)
  fit <- survey::svyglm(Math ~ (logBLL + I(logBLL^2)) * Male + PIR, design = des)

  a <- adrf(fit, treat = "logBLL", n = 5, range = c(0, 2))

  expect_equal(get_est(a),
               manual_adrf(fit, d, "logBLL", get_values(a), w = d$sw),
               tolerance = 1e-8)

  # The design-based variance assumes the target population is the one the
  # design represents, so other target weights cannot be supplied
  expect_identical(get_est(adrf(fit, treat = "logBLL", n = 5, range = c(0, 2), wts = TRUE)),
                   get_est(a))
  expect_error(adrf(fit, treat = "logBLL", n = 5, wts = FALSE),
               "cannot be supplied with <svyglm> models")
  expect_error(adrf(fit, treat = "logBLL", n = 5, wts = rep(1, nrow(d))),
               "cannot be supplied with <svyglm> models")
  expect_error(shift_curve(fit, treat = "logBLL", n = 3, wts = "PIR"),
               "cannot be supplied with <svyglm> models")
})

# For a model with only the treatment, every unit has the same prediction, so
# the variance of the ADRF is that of the coefficients, which survey computes
# from the design
svy_treat_only_se <- function(fit, a) {
  G <- cbind(1, get_values(a))
  sqrt(diag(G %*% stats::vcov(fit) %*% t(G)))
}

test_that("svyglm unconditional variance handles missing values", {
  d <- svy_data()
  d$logBLL[c(3, 50, 120)] <- NA

  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)

  # Calibrated designs cannot drop rows, so the units with missing values stay
  # in the design with weights of 0
  cal <- survey::calibrate(des, ~Male, c(`(Intercept)` = 2600, Male = 1300))

  for (design in list(des, cal)) {
    for (influence in c(FALSE, TRUE)) {
      fit <- survey::svyglm(Math ~ logBLL, design = design, influence = influence)
      a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2))

      expect_equal(get_se(a), svy_treat_only_se(fit, a), tolerance = 1e-8)
    }
  }

  # With covariates, the calibrated design gives the same variance as the
  # design restricted to the complete cases, which keeps the other units with
  # weights of 0
  f <- Math ~ logBLL * Male + PIR
  a_na <- adrf(survey::svyglm(f, design = cal), treat = "logBLL", n = 3, range = c(0, 2))
  a_cc <- adrf(survey::svyglm(f, design = subset(cal, !is.na(logBLL))),
               treat = "logBLL", n = 3, range = c(0, 2))

  expect_equal(get_est(a_na), get_est(a_cc), tolerance = 1e-8)
  expect_equal(get_vcov(a_na), get_vcov(a_cc), tolerance = 1e-8)
})

test_that("svyglm unconditional variance works with two-phase designs", {
  d <- svy_data()
  d$id <- seq_len(nrow(d))
  d$phase2 <- d$Male == 1 | d$id %% 3L == 0L

  des2 <- survey::twophase(id = list(~id, ~id), strata = list(NULL, ~Male),
                           subset = ~phase2, data = d)
  fit <- survey::svyglm(Math ~ logBLL, design = des2)
  a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2))

  expect_equal(get_se(a), svy_treat_only_se(fit, a), tolerance = 1e-8)
})

test_that("svyglm unconditional variance does not use variables in the design", {
  d <- svy_data()
  f <- Math ~ logBLL * Male + PIR

  a <- adrf(survey::svyglm(f, design = survey::svydesign(ids = ~1, weights = ~sw, data = d)),
            treat = "logBLL", n = 3, range = c(0, 2))

  d$.adrftools_z1 <- 0
  a_z <- adrf(survey::svyglm(f, design = survey::svydesign(ids = ~1, weights = ~sw, data = d)),
              treat = "logBLL", n = 3, range = c(0, 2))

  expect_identical(get_vcov(a_z), get_vcov(a))
})

test_that("svyglm unconditional variance matches the equivalent weighted lm", {
  # The design-based variance of svyglm has no HC2 adjustment, so it is
  # compared with the unadjusted variance of the weighted lm
  local_mocked_bindings(.hc2_adjust = function(infl, model) infl)

  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)
  fit <- survey::svyglm(Math ~ (logBLL + I(logBLL^2)) * Male + PIR, design = des)
  fit_lm <- lm(Math ~ (logBLL + I(logBLL^2)) * Male + PIR, data = d, weights = sw)

  a <- adrf(fit, treat = "logBLL", n = 5, range = c(0, 2))
  a_lm <- adrf(fit_lm, treat = "logBLL", n = 5, range = c(0, 2), wts = d$sw)

  expect_equal(get_est(a), get_est(a_lm), tolerance = 1e-8)
  # The design-based variance of a with-replacement design includes n / (n - 1)
  n <- nrow(d)
  expect_equal(get_vcov(a), get_vcov(a_lm) * n / (n - 1), tolerance = 1e-6)

  # svyglm's stored influence functions give the same answer
  fit_i <- survey::svyglm(Math ~ (logBLL + I(logBLL^2)) * Male + PIR, design = des,
                          influence = TRUE)
  a_i <- adrf(fit_i, treat = "logBLL", n = 5, range = c(0, 2))
  expect_equal(get_vcov(a_i), get_vcov(a), tolerance = 1e-6)
})

test_that("svyglm conditional variance uses the design-based coefficient covariance", {
  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)
  fit <- survey::svyglm(Math ~ logBLL + Male + PIR, design = des)

  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional")

  vals <- get_values(a)
  G <- cbind(1, vals, weighted.mean(d$Male, d$sw), weighted.mean(d$PIR, d$sw))
  expect_equal(get_vcov(a), unname(G %*% vcov(fit) %*% t(G)), tolerance = 1e-8)
})

test_that("svyglm logistic ADRFs work with by and downstream functions", {
  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)
  fit <- survey::svyglm(I(Block >= 14) ~ (logBLL + I(logBLL^2)) * Male + PIR, design = des,
                        family = quasibinomial)

  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male)

  sub <- d$Male == 1
  expect_equal(get_est(a)[12:22],
               manual_adrf(fit, d[sub, ], "logBLL", get_values(a), w = d$sw[sub]),
               tolerance = 1e-8)

  s <- summary(a(logBLL = c(0, 1)), simultaneous = FALSE)
  expect_true(all(s$conf.low > 0 & s$conf.high < 1))

  expect_s3_class(plot(curve_contrast(a)), "ggplot")
  expect_s3_class(summary(curve_contrast(a), method = "satterthwaite"),
                  "summary.effect_curve")
})

test_that("svyglm unconditional variance does not depend on the subset size", {
  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)

  # Treatment-only model, which uses the single-row shortcut
  fit <- survey::svyglm(Math ~ logBLL, design = des)
  a <- adrf(fit, treat = "logBLL", n = 3)
  G <- cbind(1, get_values(a))
  expect_equal(get_se(a), sqrt(diag(G %*% vcov(fit) %*% t(G))), tolerance = .01)

  # Subset; the weighted lm is compared without its HC2 adjustment, which
  # svyglm does not have
  local_mocked_bindings(.hc2_adjust = function(infl, model) infl)

  fit2 <- survey::svyglm(Math ~ logBLL + PIR + Male, design = des)
  fit2_lm <- lm(Math ~ logBLL + PIR + Male, data = d, weights = sw)
  a2 <- adrf(fit2, treat = "logBLL", n = 3, subset = Male == 1)
  a2_lm <- adrf(fit2_lm, treat = "logBLL", n = 3, subset = Male == 1, wts = d$sw)
  n <- nrow(d)
  expect_equal(get_vcov(a2), get_vcov(a2_lm) * n / (n - 1), tolerance = 1e-6)
})

test_that("svyglm unconditional variance uses the survey design", {
  d <- svy_data()
  set.seed(5)
  d$stratum <- as.integer(cut(rank(d$PIR, ties.method = "first"), 12))
  d$psu <- interaction(d$stratum, sample(8, nrow(d), replace = TRUE), drop = TRUE)

  des <- survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~sw, data = d,
                           nest = TRUE)
  f <- Math ~ logBLL * Male + PIR

  a <- adrf(survey::svyglm(f, design = des), treat = "logBLL", n = 3, range = c(0, 2))

  # Design-based jackknife: refit the model and recompute the ADRF in each
  # replicate
  rdes <- survey::as.svrepdesign(des, type = "JKn")
  jk <- survey::withReplicates(rdes, function(w, data) {
    data$.w <- w
    fit <- lm(f, data = data, weights = .w)
    vapply(c(0, 1, 2), function(v) {
      nd <- data
      nd$logBLL <- v
      weighted.mean(predict(fit, newdata = nd), w)
    }, numeric(1L))
  })

  # Linearization and the jackknife differ slightly in finite samples
  expect_equal(get_se(a), sqrt(diag(attr(jk, "var"))), tolerance = .02)

  # The same model fit to the replicate-weight design gives the same estimates
  # and, for a total of linear influence values, the same variance
  a_rep <- adrf(survey::svyglm(f, design = rdes), treat = "logBLL", n = 3,
                range = c(0, 2))
  expect_equal(get_est(a_rep), get_est(a), tolerance = 1e-8)
  # survey centers replicate variances at the mean of the replicates
  expect_equal(get_vcov(a_rep), get_vcov(a), tolerance = 1e-4)

  expect_warning(adrf(survey::svyglm(f, design = des), treat = "logBLL", n = 3,
                      cluster = ~psu),
                 "survey design is used instead")
})

test_that("svyglm models cannot be bootstrapped", {
  skip_if_not_installed("fwb")

  d <- svy_data()
  des <- survey::svydesign(ids = ~1, weights = ~sw, data = d)
  fit <- survey::svyglm(Math ~ logBLL + PIR, design = des)

  expect_error(adrf(fit, treat = "logBLL", vcov = "fwb"), "survey design")
})
