skip_if_not_installed("fwb", "1.0.0")

boot_fit <- function() {
  lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
}

test_that("bootstrap point estimates equal the original estimates", {
  fit <- boot_fit()

  set.seed(123)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "fwb",
            fwb.args = list(R = 199, verbose = FALSE))

  expect_identical(as.character(attr(a, ".vcov_type")), "bootstrap")
  expect_s3_class(attr(a, ".boot"), "fwb")
  expect_identical(dim(attr(a, ".boot")$t), c(199L, 11L))

  a0 <- adrf(fit, treat = "logBLL", n = 11)
  expect_equal(get_est(a), get_est(a0), tolerance = 1e-8)
  expect_equal(get_vcov(a), unname(cov(attr(a, ".boot")$t)), tolerance = 1e-8)

  # The traditional bootstrap uses multinomial weights
  set.seed(123)
  a_b <- adrf(fit, treat = "logBLL", n = 5, vcov = "boot",
              fwb.args = list(R = 20, verbose = FALSE))
  expect_identical(attr(a_b, ".boot")$wtype, "multinom")
  expect_true(all(attr(a_b, ".boot")$t >= 0))


  expect_error(adrf(fit, treat = "logBLL", vcov = "fwb", fwb.args = 1))
  expect_error(adrf(fit, treat = "logBLL", vcov = "fwb", cluster = ~Race + Male),
               "Only one clustering variable")
})

test_that("bootstrap inference works throughout", {
  fit <- boot_fit()

  set.seed(123)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "fwb", by = ~Male,
            fwb.args = list(R = 199, verbose = FALSE))

  est <- a(logBLL = c(0, 1, 2), subset = Male == 1)

  # Percentile intervals by default
  s <- summary(est)
  expect_identical(attr(s, "ci.type"), "perc")
  expect_true(all(s$conf.low < s$estimate & s$estimate < s$conf.high))
  expect_output(print(s), "bootstrap, simultaneous")

  s_p <- summary(est, simultaneous = FALSE, null = 8)
  expect_true(all(s_p$conf.high - s_p$conf.low < s$conf.high - s$conf.low))
  expect_false(anyNA(s_p$p.value))

  s_w <- summary(est, ci.type = "wald", simultaneous = FALSE)
  se <- sqrt(diag(cov(attr(est, ".boot")$t)))
  expect_equal(s_w$std.error, unname(se), tolerance = 1e-8)

  s_bc <- summary(est, ci.type = "bc", simultaneous = FALSE)
  expect_identical(attr(s_bc, "ci.type"), "bc")
  expect_error(summary(est, ci.type = "bc", simultaneous = TRUE))

  # Wald intervals for a transformed bootstrap ADRF
  gfit <- glm(I(Block >= 14) ~ logBLL + Male, data = nhanes3lead, family = quasibinomial)
  set.seed(1)
  ag <- adrf(gfit, treat = "logBLL", n = 5, vcov = "fwb",
             fwb.args = list(R = 199, verbose = FALSE))
  sg <- summary(ag(logBLL = c(0, 1)), ci.type = "wald")
  expect_true(all(sg$conf.low > 0))

  # Curve modifiers propagate the bootstrap draws
  m <- amef(a)
  expect_identical(dim(attr(m, ".boot")$t), c(199L, 22L))
  expect_s3_class(summary(m(logBLL = 1)), "summary.curve_est")

  r <- reference_curve(a, reference = 1)
  expect_identical(dim(attr(r, ".boot")$t), c(199L, 24L))
  expect_equal(r(logBLL = 1)$estimate, c(0, 0), tolerance = 1e-10)

  # fwb warns that simultaneous percentile bands with few replicates use
  # extreme order statistics
  expect_s3_class(suppressWarnings(plot(a)), "ggplot")
  expect_s3_class(plot(a, ci.type = "wald", subset = Male == 0), "ggplot")
  expect_s3_class(plot(a, simultaneous = FALSE, ci.type = "bc"), "ggplot")
})

test_that("bootstrap omnibus tests and projections", {
  skip_if_not_installed("CompQuadForm")

  fit <- boot_fit()

  set.seed(123)
  a <- adrf(fit, treat = "logBLL", n = 11, vcov = "fwb",
            fwb.args = list(R = 199, verbose = FALSE))

  s_flat <- summary(a)
  expect_lt(s_flat$p.value, .05)

  pc <- point_contrast(a(logBLL = c(0, 1, 2)))
  spc <- summary(pc)
  expect_identical(nrow(spc), 3L)
  spcw <- summary(pc, ci.type = "wald", simultaneous = FALSE)
  expect_true("std.error" %in% names(spcw))

  s_num <- summary(a, hypothesis = 8)
  expect_true(is.finite(s_num$p.value))

  p <- curve_projection(a, "linear")
  expect_s3_class(attr(p, ".proj_boot"), "fwb")
  expect_identical(dim(attr(p, ".proj_boot")$t), c(199L, 2L))

  sp <- summary(p)
  expect_identical(attr(sp, "ci.type"), "perc")
  expect_false(anyNA(sp$p.value))

  sp_w <- summary(p, ci.type = "wald")
  expect_equal(sp_w$std.error, unname(sqrt(diag(vcov(p)))), tolerance = 1e-8)

  expect_s3_class(plot(a, proj = p, simultaneous = FALSE), "ggplot")
  expect_s3_class(anova(curve_projection(a, "quadratic"), p), "anova.curve_projection")
})

test_that("bootstrapping refits the model with its own weights", {
  set.seed(1)
  d <- transform(nhanes3lead, w = rexp(nrow(nhanes3lead)))

  # Model weights are kept
  fit_w <- lm(Math ~ logBLL * Male + PIR, data = d, weights = w)
  set.seed(2)
  a_w <- adrf(fit_w, treat = "logBLL", n = 3, vcov = "fwb",
              fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_est(a_w), get_est(adrf(fit_w, treat = "logBLL", n = 3)),
               tolerance = 1e-8)

  # Target weights are not used as model weights
  fit_u <- lm(Math ~ logBLL * Male + PIR, data = d)
  set.seed(2)
  a_u <- adrf(fit_u, treat = "logBLL", n = 3, wts = d$w, vcov = "fwb",
              fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_est(a_u), get_est(adrf(fit_u, treat = "logBLL", n = 3, wts = d$w)),
               tolerance = 1e-8)

  # Every replicate refits with the original model weights times the
  # bootstrap weights
  set.seed(3)
  fb <- fwb::fwb(d, function(data, w) {
    data$.w <- data$w * w
    f <- lm(Math ~ logBLL * Male + PIR, data = data, weights = .w)
    vapply(get_values(a_w), function(v) {
      nd <- data
      nd$logBLL <- v
      weighted.mean(predict(f, newdata = nd), w)
    }, numeric(1L))
  }, R = 20, verbose = FALSE)
  set.seed(3)
  a_w3 <- adrf(fit_w, treat = "logBLL", n = 3, vcov = "fwb",
               fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(unname(attr(a_w3, ".boot")$t), unname(fb$t), tolerance = 1e-8)
})

test_that("bootstrapping works for models fit inside functions and with missing data", {
  fitter <- function(dat) lm(Math ~ logBLL + PIR, data = dat)
  fit <- fitter(nhanes3lead)

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "fwb",
            fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_est(a), get_est(adrf(fit, treat = "logBLL", n = 3)), tolerance = 1e-8)

  d <- nhanes3lead
  d$PIR[1:30] <- NA
  fit_na <- lm(Math ~ logBLL + PIR, data = d)

  set.seed(1)
  a_na <- adrf(fit_na, treat = "logBLL", n = 3, vcov = "fwb",
               fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_est(a_na), get_est(adrf(fit_na, treat = "logBLL", n = 3)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a_na))))
})

test_that("bootstrapping glm_weightit models re-estimates the weights", {
  skip_if_not_installed("WeightIt")

  W <- WeightIt::weightit(logBLL ~ PIR + Male, data = nhanes3lead, method = "glm")
  fit <- WeightIt::glm_weightit(Math ~ logBLL + I(logBLL^2), data = nhanes3lead,
                                weightit = W)

  set.seed(2)
  a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), vcov = "fwb",
            fwb.args = list(R = 50, verbose = FALSE))

  # WeightIt's own bootstrap re-estimates the weights in each replicate; with
  # the same seed it gives the same coefficients, so the same ADRF
  set.seed(2)
  fit_fwb <- WeightIt::glm_weightit(Math ~ logBLL + I(logBLL^2), data = nhanes3lead,
                                    weightit = W, vcov = "FWB", R = 50,
                                    fwb.args = list(verbose = FALSE))
  G <- cbind(1, get_values(a), get_values(a)^2)

  expect_equal(get_vcov(a), unname(G %*% vcov(fit_fwb) %*% t(G)), tolerance = 1e-6)
})

test_that("cluster bootstrap resamples whole clusters", {
  d <- transform(nhanes3lead, cl = rep(seq_len(150), length.out = nrow(nhanes3lead)))
  fit <- lm(Math ~ logBLL * Male + PIR, data = d)

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), vcov = "fwb", cluster = ~cl,
            fwb.args = list(R = 199, verbose = FALSE))

  expect_identical(length(attr(a, ".boot")$cluster), nrow(d))

  # Close to the cluster-robust unconditional variance
  a_u <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), cluster = ~cl)
  expect_equal(get_se(a), get_se(a_u), tolerance = .15)

  # A cluster formula can contain a transformation of variables, such as their
  # interaction, and gives the same result as the variable computed beforehand
  # (this is one-way clustering on the combinations, not multiway clustering)
  d$RM <- interaction(d$Race, d$Male)
  fit2 <- lm(Math ~ logBLL * Male + PIR, data = d)
  set.seed(2)
  a_int <- adrf(fit, treat = "logBLL", n = 3, vcov = "fwb", cluster = ~interaction(Race, Male),
                fwb.args = list(R = 20, verbose = FALSE))
  set.seed(2)
  a_rm <- adrf(fit2, treat = "logBLL", n = 3, vcov = "fwb", cluster = ~RM,
               fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_vcov(a_int), get_vcov(a_rm), tolerance = 1e-10)

  # More than one clustering variable is an error
  expect_error(adrf(fit, treat = "logBLL", vcov = "boot", cluster = d[c("Race", "Male")]),
               "Only one clustering variable")
})

test_that("clusters are aligned with the rows used in the fit", {
  d <- nhanes3lead
  d$PIR[1:20] <- NA
  fit <- lm(Math ~ logBLL + PIR, data = d)

  a <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Race)
  expect_true(all(is.finite(get_se(a))))

  set.seed(1)
  a_b <- adrf(fit, treat = "logBLL", n = 3, vcov = "fwb", cluster = ~Race,
              fwb.args = list(R = 20, verbose = FALSE))
  expect_identical(length(attr(a_b, ".boot")$cluster), nobs(fit))

  d$Race[25] <- NA
  fit_na <- lm(Math ~ logBLL + PIR, data = d)
  expect_error(adrf(fit_na, treat = "logBLL", n = 3, cluster = ~Race), "not allowed")
})

test_that("curve contrasts and flat projections carry the bootstrap draws", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Male, vcov = "fwb",
            fwb.args = list(R = 50, verbose = FALSE))
  t_a <- attr(a, ".boot")$t

  # Each replicate of the contrast is the difference of the subgroup curves
  cc <- curve_contrast(a)
  expect_identical(dim(attr(cc, ".boot")$t), c(50L, 5L))
  expect_equal(unname(attr(cc, ".boot")$t), unname(t_a[, 6:10] - t_a[, 1:5]),
               tolerance = 1e-10)
  expect_equal(get_vcov(cc), unname(cov(t_a[, 6:10] - t_a[, 1:5])), tolerance = 1e-10)
  expect_s3_class(summary(cc(logBLL = 1)), "summary.curve_est")

  # Each replicate of a flat projection is the trapezoidal average of the curve
  set.seed(1)
  a1 <- adrf(fit, treat = "logBLL", n = 5, vcov = "fwb",
             fwb.args = list(R = 50, verbose = FALSE))
  p <- curve_projection(a1, "flat")
  vals <- get_values(a1)
  w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
  w <- w / sum(w)

  expect_identical(dim(attr(p, ".proj_boot")$t), c(50L, 1L))
  expect_equal(drop(attr(p, ".proj_boot")$t), unname(drop(attr(a1, ".boot")$t %*% w)),
               tolerance = 1e-8)
  expect_s3_class(summary(p), "summary.curve_projection")
})

test_that("point contrasts within subgroups carry the bootstrap draws", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", n = 5, by = ~Male, vcov = "fwb",
            fwb.args = list(R = 50, verbose = FALSE))

  est <- a(logBLL = c(0, 2))
  pc <- point_contrast(est)

  t_e <- attr(est, ".boot")$t
  expect_equal(unname(attr(pc, ".boot")$t), unname(t_e[, c(2, 4)] - t_e[, c(1, 3)]),
               tolerance = 1e-10)

  s <- summary(pc, simultaneous = FALSE)
  expect_identical(nrow(s), 2L)
  expect_false(anyNA(s$p.value))
})

test_that("bootstrapping works with weights and na.action = na.exclude", {
  d <- nhanes3lead
  d$PIR[1:20] <- NA
  set.seed(1)
  d$w <- runif(nrow(d))

  f_excl <- lm(Math ~ logBLL + PIR, data = d, weights = w, na.action = na.exclude)
  f_omit <- lm(Math ~ logBLL + PIR, data = d, weights = w)

  set.seed(9)
  a_excl <- adrf(f_excl, treat = "logBLL", n = 3, vcov = "boot",
                 fwb.args = list(R = 50, verbose = FALSE))
  set.seed(9)
  a_omit <- adrf(f_omit, treat = "logBLL", n = 3, vcov = "boot",
                 fwb.args = list(R = 50, verbose = FALSE))

  expect_equal(get_vcov(a_excl), get_vcov(a_omit), tolerance = 1e-10)
})

test_that("covariances of transformed bootstrap estimates equal those of the transformed replicates", {
  fit <- boot_fit()

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", by = ~Male, vcov = "fwb",
            fwb.args = list(R = 100, verbose = FALSE))

  cov_t <- function(x) unname(cov(attr(x, ".boot")$t))

  e <- a(logBLL = c(.3, 1.1, 2.2))
  expect_equal(unname(get_vcov(e)), cov_t(e), tolerance = 1e-8)

  m <- amef(a)(logBLL = c(0, 1))
  expect_equal(unname(get_vcov(m)), cov_t(m), tolerance = 1e-8)

  r <- reference_curve(a, 1)(logBLL = c(0, 2))
  expect_equal(unname(get_vcov(r)), cov_t(r), tolerance = 1e-8)

  pc <- point_contrast(a(logBLL = c(0, 1, 2)))
  expect_equal(unname(vcov(pc)), cov_t(pc), tolerance = 1e-8)

  # Estimates with no variance (the reference point) have no p-value
  s <- summary(reference_curve(a, 1)(logBLL = c(.5, 1, 2)), simultaneous = FALSE)
  expect_true(all(is.na(s$p.value[c(2, 5)])))
  expect_true(all(!is.na(s$p.value[-c(2, 5)])))
})

test_that("bootstrap replicates that alias a coefficient match refitting the model", {
  d <- nhanes3lead
  set.seed(2)
  d$g <- factor(ifelse(seq_len(nrow(d)) <= 2L, "rare", sample(c("a", "b"), nrow(d), TRUE)))
  fit <- lm(Math ~ poly(logBLL, 2) + Male + g, data = d)

  set.seed(1)
  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "boot",
            fwb.args = list(R = 100, verbose = FALSE))
  vals <- get_values(a)

  # The same resampling weights, with the model refit and predictions from
  # predict(); the last column flags replicates in which "rare" gets no weight
  stat <- function(data, w) {
    refit <- lm(Math ~ poly(logBLL, 2) + Male + g, data = data, weights = w)
    c(vapply(vals, function(v) {
      nd <- data
      nd$logBLL <- v
      weighted.mean(suppressWarnings(predict(refit, nd)), w)
    }, numeric(1L)),
    all(w[1:2] == 0))
  }

  set.seed(1)
  f <- fwb::fwb(d, stat, R = 100, wtype = "multinom", drop0 = FALSE, verbose = FALSE)

  expect_gt(sum(f$t[, 4L]), 0)
  expect_equal(unname(attr(a, ".boot")$t), unname(f$t[, 1:3]), tolerance = 1e-10)
})
