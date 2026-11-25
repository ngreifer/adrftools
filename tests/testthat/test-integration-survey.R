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

  expect_output(print(a), "inference: unconditional \\(survey design\\)")

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

# The design-based standard errors of the ADRF of `treat` at `vals` for a
# linear svyglm model, by Binder's linearization computed directly from the
# influence matrix survey saves with `influence = TRUE` (one row per data row,
# with rows of 0 for units dropped for missing values). `var_fun` computes the
# design variance of a vector of per-unit contributions. Units the fit did not
# use have weights of 0 in the average.
svy_linearized_se <- function(fit_inf, treat, vals, var_fun) {
  des <- fit_inf$survey.design
  w <- stats::weights(des)
  infl <- attr(fit_inf, "influence")
  mf <- stats::model.frame(fit_inf)
  tt <- stats::delete.response(stats::terms(fit_inf))
  na_act <- as.integer(fit_inf$na.action)

  # Designs that keep the incomplete units (e.g., calibrated and PPS designs)
  # have a row for every data row; others drop them, but the influence matrix
  # is still padded
  if (nrow(des) == nrow(mf)) {
    rows <- seq_len(nrow(des))

    if (length(na_act) > 0L && nrow(infl) == nrow(mf) + length(na_act)) {
      infl <- infl[-na_act, , drop = FALSE]
    }
  }
  else {
    rows <- setdiff(seq_len(nrow(des)), na_act)
  }

  vapply(vals, function(v) {
    nd <- mf
    nd[[treat]] <- v

    Xv <- matrix(0, nrow(des), ncol(infl))
    Xv[rows, ] <- stats::model.matrix(tt, nd)

    p <- drop(Xv %*% stats::coef(fit_inf))

    wr <- w
    wr[-rows] <- 0

    theta <- sum(wr * p) / sum(wr)
    z <- wr * (p - theta) / sum(wr) + drop(infl %*% (colSums(wr * Xv) / sum(wr)))

    sqrt(var_fun(z, des))
  }, numeric(1L))
}

test_that("svyglm unconditional variance matches a hand linearization on calibrated, domain, and PPS designs", {
  skip_if_not_installed("patrick")

  data("api", package = "survey", envir = environment())

  # Variance of a total of per-unit contributions for linearization designs,
  # including calibration
  rec_var <- function(z, des) {
    as.numeric(survey::svyrecvar(z, des$cluster, des$strata, des$fpc,
                                 postStrata = des$postStrata))
  }

  # svyrecvar() does not handle PPS designs, so svytotal() is used, with the
  # per-unit values divided by the weights it multiplies them by
  total_var <- function(z, des) {
    w <- stats::weights(des)
    zz <- ifelse(w == 0, 0, z / w)
    as.numeric(stats::vcov(survey::svytotal(matrix(zz), des)))
  }

  dat_na <- apistrat
  dat_na$ell[c(3, 50, 120)] <- NA

  strat <- function(dat) {
    survey::svydesign(id = ~1, strata = ~stype, weights = ~pw, data = dat, fpc = ~fpc)
  }

  # Calibrating to the population total of api99 moves every weight
  cal <- function(des) {
    survey::calibrate(des, ~stype + api99,
                      c(`(Intercept)` = 6194, stypeH = 755, stypeM = 1018,
                        api99 = 3914069))
  }

  expect_true(all(abs(stats::weights(cal(strat(dat_na))) /
                        stats::weights(strat(dat_na)) - 1) > 1e-8))

  # survey warns that units with weights of 0 are not used for the dispersion
  quiet_svyglm <- function(...) {
    withCallingHandlers(survey::svyglm(...), warning = function(w) {
      if (grepl("zero weight", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }

  patrick::with_parameters_test_that(
    "{label}",
    {
      f <- api00 ~ ell + meals + stype
      vals <- c(10, 30, 50)

      fit <- quiet_svyglm(f, design = design)
      fit_inf <- quiet_svyglm(f, design = design, influence = TRUE)

      ref <- svy_linearized_se(fit_inf, "ell", vals, var_fun)

      a <- adrf(fit, treat = "ell", n = 3, range = c(10, 50))
      a_inf <- adrf(fit_inf, treat = "ell", n = 3, range = c(10, 50))

      expect_equal(get_values(a), vals)
      expect_equal(get_se(a), ref, tolerance = 1e-8)
      expect_equal(get_se(a_inf), ref, tolerance = 1e-8)
    },
    patrick::cases(
      list(label = "stratified design with missing values",
           design = strat(dat_na), var_fun = rec_var),
      list(label = "calibrated design with missing values",
           design = cal(strat(dat_na)), var_fun = rec_var),
      list(label = "domain of a calibrated design with missing values",
           design = subset(cal(strat(dat_na)), sch.wide == "Yes"), var_fun = rec_var),
      list(label = "domain of a calibrated design",
           design = subset(cal(strat(apistrat)), sch.wide == "Yes"), var_fun = rec_var),
      list(label = "PPS design with missing values",
           design = survey::svydesign(id = ~1, fpc = ~I(pw^-1), data = dat_na, pps = "brewer"),
           var_fun = total_var),
      list(label = "cluster sample with a finite population correction",
           design = survey::svydesign(id = ~dnum, weights = ~pw, data = apiclus1, fpc = ~fpc),
           var_fun = rec_var)
    )
  )
})
