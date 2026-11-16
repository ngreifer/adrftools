test_that("adrf() point estimates equal average counterfactual predictions for lm", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR + Race, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 7)

  expect_s3_class(a, c("adrf_curve", "effect_curve", "function"), exact = TRUE)
  expect_length(get_values(a), 7L)
  expect_equal(get_est(a),
               manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_identical(attr(a, ".treat"), "logBLL")
  expect_identical(attr(a, ".response"), "Math")
  expect_identical(as.character(attr(a, ".vcov_type")), "unconditional")
  expect_equal(attr(a, ".df"), df.residual(fit))
})

test_that("adrf() point estimates equal average counterfactual predictions for glm", {
  fit <- glm(I(Block >= 12) ~ poly(logBLL, 2) + Male + PIR,
             data = nhanes3lead, family = binomial)

  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_equal(get_est(a),
               manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_identical(attr(a, ".family")$link, "logit")

  # Poisson with log link goes through the same glm path
  fit_p <- glm(Block ~ logBLL + Male + PIR, data = nhanes3lead, family = poisson)
  a_p <- adrf(fit_p, treat = "logBLL", n = 5)

  expect_equal(get_est(a_p),
               manual_adrf(fit_p, nhanes3lead, "logBLL", get_values(a_p)),
               tolerance = 1e-8)
})

test_that("adrf() with type = 'link' averages the linear predictor", {
  fit <- glm(I(Block >= 12) ~ logBLL + Male + PIR,
             data = nhanes3lead, family = binomial)

  a <- adrf(fit, treat = "logBLL", n = 5, type = "link")

  expect_equal(get_est(a),
               manual_adrf(fit, nhanes3lead, "logBLL", get_values(a), type = "link"),
               tolerance = 1e-6)
  expect_identical(attr(a, ".family")$link, "identity")
})

test_that("conditional variance matches marginaleffects::avg_predictions()", {
  skip_if_not_installed("marginaleffects")

  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 4, vcov = "conditional")

  expect_identical(as.character(attr(a, ".vcov_type")), "conditional")

  me <- marginaleffects::avg_predictions(fit,
                                         variables = list(logBLL = get_values(a)),
                                         by = "logBLL")

  expect_equal(get_est(a), me$estimate, tolerance = 1e-8)
  expect_equal(get_se(a), me$std.error, tolerance = 1e-6)

  # String vcov types are passed to marginaleffects::get_vcov()
  a_hc3 <- adrf(fit, treat = "logBLL", n = 4, vcov = "HC3")
  me_hc3 <- marginaleffects::avg_predictions(fit, vcov = "HC3",
                                             variables = list(logBLL = get_values(a)),
                                             by = "logBLL")

  expect_equal(get_se(a_hc3), me_hc3$std.error, tolerance = 1e-6)

  # A numeric covariance matrix
  a_mat <- adrf(fit, treat = "logBLL", n = 4, vcov = sandwich::vcovHC(fit, type = "HC3"))
  expect_equal(get_vcov(a_mat), get_vcov(a_hc3), tolerance = 1e-8)
})

test_that("unconditional variance matches a hand-computed influence function", {
  fit <- lm(Math ~ logBLL * Male + PIR + Race, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_equal(get_vcov(a),
               manual_uncond_vcov_lm(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)

  # Unconditional variance exceeds conditional variance because it adds the
  # sampling variability of the covariates
  a_c <- adrf(fit, treat = "logBLL", n = 5, vcov = "conditional")
  expect_true(all(get_se(a) > get_se(a_c)))
})

test_that("unconditional variance is the same for glm-path and numerical-gradient models", {
  fit <- glm(I(Block >= 12) ~ logBLL + Male + PIR,
             data = nhanes3lead, family = binomial)

  # The glm path uses the analytic gradient; type = "link" uses a numerical
  # gradient through marginaleffects
  a <- adrf(fit, treat = "logBLL", n = 4)

  fit_q <- glm(I(Block >= 12) ~ logBLL + Male + PIR,
               data = nhanes3lead, family = quasibinomial)
  a_q <- adrf(fit_q, treat = "logBLL", n = 4)

  expect_equal(get_est(a), get_est(a_q), tolerance = 1e-8)
  expect_equal(get_vcov(a), get_vcov(a_q), tolerance = 1e-6)
})

test_that("cluster-robust unconditional variance matches a hand computation", {
  fit <- lm(Math ~ logBLL + PIR + Male, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Race)

  expect_true(attr(attr(a, ".vcov_type"), "clustered"))

  # Recompute from the unclustered influence function contributions
  N <- nrow(nhanes3lead)
  IFb <- sandwich::estfun(fit) %*% sandwich::bread(fit)
  vals <- get_values(a)
  G <- cbind(1, vals, mean(nhanes3lead$PIR), mean(nhanes3lead$Male))
  P <- vapply(vals, function(v) predict(fit, transform(nhanes3lead, logBLL = v)),
              numeric(N))
  IF <- tcrossprod(IFb, G) + sweep(P, 2L, colMeans(P))
  cl <- nhanes3lead$Race
  g <- nlevels(cl)
  S <- rowsum(IF, cl)
  V <- crossprod(S / N) * g / (g - 1)

  expect_equal(get_vcov(a), V, tolerance = 1e-8)

  # Cluster can also be given as a data frame
  a_df <- adrf(fit, treat = "logBLL", n = 3, cluster = nhanes3lead["Race"])
  expect_equal(get_vcov(a_df), get_vcov(a), tolerance = 1e-10)
})

test_that("two-way clustering combines the cluster variances", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  a1 <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Race)
  a2 <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Male)
  a12 <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Race + Male)

  expect_false(isTRUE(all.equal(get_vcov(a12), get_vcov(a1))))
  expect_false(isTRUE(all.equal(get_vcov(a12), get_vcov(a2))))
  expect_true(attr(attr(a12, ".vcov_type"), "clustered"))
})

test_that("cluster-robust conditional variance uses sandwich::vcovCL()", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "HC1", cluster = ~Race)

  vals <- get_values(a)
  G <- cbind(1, vals, mean(nhanes3lead$PIR))
  V <- G %*% sandwich::vcovCL(fit, cluster = ~Race, type = "HC1") %*% t(G)

  expect_equal(get_vcov(a), V, tolerance = 1e-8)

  a_c <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional", cluster = ~Race)
  expect_equal(get_se(a_c), get_se(a), tolerance = 1e-6)

  # HC3 is passed to vcovCL() rather than ignored
  a3 <- adrf(fit, treat = "logBLL", n = 3, vcov = "HC3", cluster = ~Race)
  V3 <- G %*% sandwich::vcovCL(fit, cluster = ~Race, type = "HC3") %*% t(G)
  expect_equal(get_vcov(a3), V3, tolerance = 1e-8)

  expect_error(adrf(fit, treat = "logBLL", vcov = "HC4", cluster = ~Race),
               "cannot be \"HC4\" when")
  expect_error(adrf(fit, treat = "logBLL", vcov = sandwich::vcovHC, cluster = ~Race),
               "function or formula")
})

test_that("vcov = 'none' skips uncertainty", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 5, vcov = "none")

  expect_identical(as.character(attr(a, ".vcov_type")), "none")
  expect_null(get_vcov(a))
  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)

  expect_warning(adrf(fit, treat = "logBLL", vcov = "none", cluster = ~Race),
                 "ignored")
  expect_identical(as.character(attr(adrf(fit, treat = "logBLL", vcov = FALSE), ".vcov_type")),
                   "none")
  expect_identical(as.character(attr(adrf(fit, treat = "logBLL", vcov = NULL), ".vcov_type")),
                   "none")
})

test_that("vcov strings are partially matched", {
  skip_if_not_installed("patrick")

  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  patrick::with_parameters_test_that(
    "vcov = \"{v}\"",
    {
      a <- adrf(fit, treat = "logBLL", n = 3, vcov = v)
      expect_identical(as.character(attr(a, ".vcov_type")), expected)
    },
    v = c("uncond", "cond"),
    expected = c("unconditional", "conditional")
  )
})

test_that("unconditional variance falls back to conditional when sandwich methods are missing", {
  skip_if_not_installed("estimatr")

  fit <- estimatr::lm_robust(Math ~ logBLL + PIR, data = nhanes3lead)

  expect_warning(a <- adrf(fit, treat = "logBLL", n = 3),
                 "unconditional")
  expect_identical(as.character(attr(a, ".vcov_type")), "conditional")
})

test_that("range and n control the grid", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)
  tr <- nhanes3lead$logBLL

  a <- adrf(fit, treat = "logBLL", n = 11)
  expect_equal(range(get_values(a)), unname(quantile(tr, c(.025, .975))),
               tolerance = 1e-3)

  a1 <- adrf(fit, treat = "logBLL", range = 1, n = 6)
  expect_equal(get_values(a1), seq(min(tr), max(tr), length.out = 6))

  a2 <- adrf(fit, treat = "logBLL", range = c(2, 0), n = 5)
  expect_equal(get_values(a2), seq(0, 2, length.out = 5))

  expect_error(adrf(fit, treat = "logBLL", range = c(-10, 2)), "outside")
  expect_error(adrf(fit, treat = "logBLL", range = 1.5), "quantile")
  expect_error(adrf(fit, treat = "logBLL", range = c(0, 1, 2)))
  expect_error(adrf(fit, treat = "logBLL", range = "a"))
  expect_error(adrf(fit, treat = "logBLL", n = 1))
  expect_error(adrf(fit, treat = "logBLL", n = 2.5))
})

test_that("weighted quantiles define the default range when wts are supplied", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)

  w <- ifelse(nhanes3lead$logBLL > 1, 10, 1)

  a <- adrf(fit, treat = "logBLL", wts = w, n = 3)
  a0 <- adrf(fit, treat = "logBLL", n = 3)

  expect_gt(min(get_values(a)), min(get_values(a0)))
})

test_that("subset restricts the target population", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 4, subset = Male == 1)

  sub <- nhanes3lead[nhanes3lead$Male == 1, ]
  expect_equal(get_est(a), manual_adrf(fit, sub, "logBLL", get_values(a)),
               tolerance = 1e-8)

  expect_error(adrf(fit, treat = "logBLL", subset = c(TRUE, FALSE)),
               "one value for each unit")
  expect_error(adrf(fit, treat = "logBLL", subset = list(1)),
               "atomic")
  expect_error(adrf(fit, treat = "logBLL", subset = factor(nhanes3lead$Male)),
               "cannot be a factor")
  expect_error(adrf(fit, treat = "logBLL", subset = c("1", "not_a_row")),
               "not a row name")
  expect_error(adrf(fit, treat = "logBLL", subset = c(1, -2)),
               "both positive and negative")
  expect_error(adrf(fit, treat = "logBLL", subset = c(1, nrow(nhanes3lead) + 1)),
               "larger than the number of units")
  expect_error(adrf(fit, treat = "logBLL", subset = c(1.5, 2)),
               "whole numbers")

  skip_if_not_installed("patrick")

  # Other ways of selecting the same units give the same result; character
  # subsets refer to row names
  patrick::with_parameters_test_that(
    "subset as {label}",
    {
      a_s <- adrf(fit, treat = "logBLL", n = 4, subset = s)
      expect_equal(get_est(a_s), get_est(a), tolerance = 1e-10)
      expect_equal(get_vcov(a_s), get_vcov(a), tolerance = 1e-10)
    },
    patrick::cases(
      list(label = "positive indices", s = which(nhanes3lead$Male == 1)),
      list(label = "negative indices", s = -which(nhanes3lead$Male == 0)),
      list(label = "row names", s = rownames(nhanes3lead)[nhanes3lead$Male == 1])
    )
  )
})

test_that("subset can refer to variables in the calling environment", {
  fit <- lm(Math ~ logBLL + PIR + Male, data = nhanes3lead)

  sel <- nhanes3lead$PIR > 2
  a <- adrf(fit, treat = "logBLL", n = 3, subset = sel)

  expect_equal(get_est(a),
               manual_adrf(fit, nhanes3lead[sel, ], "logBLL", get_values(a)),
               tolerance = 1e-8)
})

test_that("by computes subgroup ADRFs in the correct order", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 4, by = ~Male)

  by_grid <- attr(a, ".by_grid")
  expect_identical(by_grid$Male, c(0L, 1L))

  vals <- get_values(a)
  manual <- c(manual_adrf(fit, nhanes3lead[nhanes3lead$Male == 0, ], "logBLL", vals),
              manual_adrf(fit, nhanes3lead[nhanes3lead$Male == 1, ], "logBLL", vals))

  expect_equal(get_est(a), manual, tolerance = 1e-8)
  expect_identical(dim(get_vcov(a)), c(8L, 8L))

  # Character specification is equivalent
  a_chr <- adrf(fit, treat = "logBLL", n = 4, by = "Male")
  expect_equal(get_est(a_chr), get_est(a), tolerance = 1e-10)

  expect_error(adrf(fit, treat = "logBLL", by = 1), "formula")
})

test_that("by with several variables and subset", {
  fit <- lm(Math ~ poly(logBLL, 2) * Male * Race + PIR, data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 3, by = ~Male + Race,
            subset = Race != "Other")

  by_grid <- attr(a, ".by_grid")
  expect_identical(nrow(by_grid), 6L)
  expect_false("Other" %in% as.character(by_grid$Race))

  sub <- nhanes3lead[nhanes3lead$Male == 1 & nhanes3lead$Race == "White", ]
  k <- which(by_grid$Male == 1 & by_grid$Race == "White")

  expect_equal(get_est(a)[(k - 1) * 3 + 1:3],
               manual_adrf(fit, sub, "logBLL", get_values(a)),
               tolerance = 1e-8)
})

test_that("wts generalizes the ADRF to a weighted population", {
  fit <- lm(Math ~ logBLL * Male + PIR, data = nhanes3lead)

  set.seed(123)
  w <- rexp(nrow(nhanes3lead))

  a <- adrf(fit, treat = "logBLL", n = 4, wts = w, range = c(0, 2))

  expect_equal(get_est(a),
               manual_adrf(fit, nhanes3lead, "logBLL", get_values(a), w = w),
               tolerance = 1e-8)

  # FALSE means unit weights
  a_f <- adrf(fit, treat = "logBLL", n = 4, wts = FALSE, range = c(0, 2))
  a_0 <- adrf(fit, treat = "logBLL", n = 4, range = c(0, 2))
  expect_equal(get_est(a_f), get_est(a_0), tolerance = 1e-10)

  expect_error(adrf(fit, treat = "logBLL", wts = w[-1]), "as many values")
  expect_error(adrf(fit, treat = "logBLL", wts = list(1)), "must be")
})

test_that("wts = TRUE uses the model weights", {
  set.seed(123)
  d <- transform(nhanes3lead, w = rexp(nrow(nhanes3lead)))
  fit <- lm(Math ~ logBLL * Male + PIR, data = d, weights = w)

  a_t <- adrf(fit, treat = "logBLL", n = 3, wts = TRUE, range = c(0, 2))
  a_w <- adrf(fit, treat = "logBLL", n = 3, wts = d$w, range = c(0, 2))

  expect_equal(get_est(a_t), get_est(a_w), tolerance = 1e-10)

  # String naming the model weights column
  a_s <- adrf(fit, treat = "logBLL", n = 3, wts = "(weights)", range = c(0, 2))
  expect_equal(get_est(a_s), get_est(a_w), tolerance = 1e-10)

  # By default, model weights are not used as target weights
  a_d <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2))
  expect_false(isTRUE(all.equal(get_est(a_d), get_est(a_w))))
})

test_that("a model with only the treatment uses the single-row shortcut", {
  fit <- lm(Math ~ poly(logBLL, 2), data = nhanes3lead)

  a <- adrf(fit, treat = "logBLL", n = 5, vcov = "conditional")
  a_slow <- adrf(fit, treat = "logBLL", n = 5, vcov = "conditional", .quick = FALSE)

  expect_equal(get_est(a), get_est(a_slow), tolerance = 1e-8)
  expect_equal(get_vcov(a), get_vcov(a_slow), tolerance = 1e-8)

  vals <- get_values(a)
  expect_equal(get_est(a), unname(predict(fit, data.frame(logBLL = vals))),
               tolerance = 1e-8)

  # Unconditional variance equals the HC0 conditional variance when there are
  # no covariates
  a_u <- adrf(fit, treat = "logBLL", n = 5)
  a_hc0 <- adrf(fit, treat = "logBLL", n = 5, vcov = "HC0")
  expect_equal(get_vcov(a_u), get_vcov(a_hc0), tolerance = 1e-6)
})

test_that("aliased coefficients are dropped", {
  d <- transform(nhanes3lead, PIR2 = 2 * PIR)
  fit <- lm(Math ~ logBLL + PIR + PIR2, data = d)
  fit0 <- lm(Math ~ logBLL + PIR, data = d)

  expect_true(anyNA(coef(fit)))

  a <- adrf(fit, treat = "logBLL", n = 3)
  a0 <- adrf(fit0, treat = "logBLL", n = 3)

  expect_equal(get_est(a), get_est(a0), tolerance = 1e-8)
  expect_equal(get_vcov(a), get_vcov(a0), tolerance = 1e-8)
})

test_that("models with missing data use the complete cases", {
  d <- nhanes3lead
  d$PIR[1:50] <- NA

  fit <- lm(Math ~ logBLL + PIR + Male, data = d)

  a <- adrf(fit, treat = "logBLL", n = 3)

  expect_equal(get_est(a),
               manual_adrf(fit, d[-(1:50), ], "logBLL", get_values(a)),
               tolerance = 1e-8)
})

test_that("splines and transformations of the treatment are handled", {
  skip_if_not_installed("splines")

  fit <- lm(Math ~ splines::ns(logBLL, 4) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 5)

  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)

  fit2 <- lm(Math ~ logBLL + I(logBLL^2) + log(PIR + 1), data = nhanes3lead)
  a2 <- adrf(fit2, treat = "logBLL", n = 5)

  expect_equal(get_est(a2), manual_adrf(fit2, nhanes3lead, "logBLL", get_values(a2)),
               tolerance = 1e-8)
})

test_that("negative binomial models are supported", {
  skip_if_not_installed("MASS")

  fit <- MASS::glm.nb(Block ~ logBLL + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 3)

  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))
})

test_that("models without a dedicated estimator use marginaleffects predictions", {
  skip_if_not_installed("betareg")

  d <- transform(nhanes3lead, y = (Math + .5) / (max(Math) + 1))
  fit <- betareg::betareg(y ~ logBLL + Male + PIR, data = d)

  a <- adrf(fit, treat = "logBLL", n = 3)

  expect_equal(get_est(a), manual_adrf(fit, d, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))

  a_c <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional")
  me <- marginaleffects::avg_predictions(fit, variables = list(logBLL = get_values(a)),
                                         by = "logBLL")
  expect_equal(get_se(a_c), me$std.error, tolerance = 1e-4)
})

test_that("adrf() validates its inputs", {
  fit <- lm(Math ~ logBLL + PIR + Race, data = nhanes3lead)

  expect_error(adrf(fit), "treat")
  expect_error(adrf(fit, treat = "Reading"), "must be the name of a variable")
  expect_error(adrf(fit, treat = "Race"), "numeric")
  expect_error(adrf(fit, treat = c("logBLL", "PIR")))
  expect_error(adrf(fit, treat = "logBLL", type = 1))
  expect_error(adrf(fit, treat = "logBLL", vcov = "nonsense"), "vcov")
  expect_error(suppressWarnings(adrf(fit, treat = "logBLL", vcov = diag(2))), "dimension")
  expect_error(adrf(fit, treat = "logBLL", cluster = 1:3), "number of observations")
  expect_error(adrf(fit, treat = "logBLL", cluster = ~nonexistent_var), "could not be extracted")

  d <- nhanes3lead
  d$cl <- NA
  fit2 <- lm(Math ~ logBLL + PIR, data = d)
  expect_error(adrf(fit2, treat = "logBLL", cluster = d["cl"]), "NA")
})

test_that("adrf() asks for data when it cannot be extracted", {
  fit <- lm(Math ~ logBLL + PIR, data = nhanes3lead)
  bad <- structure(list(), class = "not_a_model")

  expect_error(adrf(bad, treat = "logBLL"), "data")
})

test_that("models whose original data cannot be recovered ask for data", {
  # Fit inside a function, so insight finds only the model frame, which stores
  # the poly() basis rather than the treatment
  fitter <- function(formula, dat) lm(formula, data = dat)
  f <- local(Math ~ poly(logBLL, 2) + Male + PIR)
  fit <- fitter(f, nhanes3lead)

  expect_error(adrf(fit, treat = "logBLL", n = 3), "could not be recovered")

  a <- adrf(fit, treat = "logBLL", n = 3, data = nhanes3lead)

  expect_equal(range(get_values(a)), unname(quantile(nhanes3lead$logBLL, c(.025, .975))),
               tolerance = 1e-3)
  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)

  # Without transformations, the model frame holds the original variables
  fit2 <- fitter(local(Math ~ logBLL + Male + PIR), nhanes3lead)
  a2 <- adrf(fit2, treat = "logBLL", n = 3)
  expect_equal(get_est(a2), manual_adrf(fit2, nhanes3lead, "logBLL", get_values(a2)),
               tolerance = 1e-8)
})

test_that("weighted unconditional variance includes the unit weights", {
  set.seed(1)
  w <- rexp(nrow(nhanes3lead))^2
  fit <- lm(Math ~ logBLL + PIR + Male, data = nhanes3lead)

  N <- nrow(nhanes3lead)
  IFb <- sandwich::estfun(fit) %*% sandwich::bread(fit)
  tt <- delete.response(terms(fit))

  # Influence function of the weighted mean of predictions within the units
  # flagged by `b`: the coefficient term plus b_i * w_i * (m_i - mu) * N / sum(b * w)
  manual_vcov <- function(vals, b) {
    IF <- vapply(vals, function(v) {
      d <- transform(nhanes3lead, logBLL = v)
      X <- model.matrix(tt, data = d)
      m <- drop(X %*% coef(fit))
      mu <- sum(b * w * m) / sum(b * w)
      g <- colSums(b * w * X) / sum(b * w)
      drop(IFb %*% g) + b * w * (m - mu) * N / sum(b * w)
    }, numeric(N))

    crossprod(IF) / N^2
  }

  a <- adrf(fit, treat = "logBLL", n = 4, wts = w, range = c(0, 2))
  expect_equal(get_vcov(a), manual_vcov(get_values(a), rep(1, N)), tolerance = 1e-8)

  a_by <- adrf(fit, treat = "logBLL", n = 4, wts = w, range = c(0, 2), by = ~Male)
  expect_equal(get_vcov(a_by)[5:8, 5:8],
               manual_vcov(get_values(a), nhanes3lead$Male == 1),
               tolerance = 1e-8)
})

test_that("wts, subset, and by can refer to variables outside the model", {
  set.seed(1)
  d <- transform(nhanes3lead, w = rexp(nrow(nhanes3lead)))
  d$PIR[1:20] <- NA
  fit <- lm(Math ~ logBLL + PIR, data = d)
  used <- d[!is.na(d$PIR), ]

  a_w <- adrf(fit, treat = "logBLL", n = 3, wts = "w", range = c(0, 2))
  expect_equal(get_est(a_w), manual_adrf(fit, used, "logBLL", get_values(a_w), w = used$w),
               tolerance = 1e-8)

  a_s <- adrf(fit, treat = "logBLL", n = 3, subset = Race == "White", range = c(0, 2))
  expect_equal(get_est(a_s),
               manual_adrf(fit, used[used$Race == "White", ], "logBLL", get_values(a_s)),
               tolerance = 1e-8)

  a_b <- adrf(fit, treat = "logBLL", n = 3, by = ~Race, range = c(0, 2))
  expect_identical(levels(attr(a_b, ".by_grid")$Race), levels(d$Race))
  k <- which(attr(a_b, ".by_grid")$Race == "Black")
  expect_equal(get_est(a_b)[(k - 1) * 3 + 1:3],
               manual_adrf(fit, used[used$Race == "Black", ], "logBLL", get_values(a_b)),
               tolerance = 1e-8)

  a_c <- adrf(fit, treat = "logBLL", n = 3, cluster = ~Race)
  expect_true(attr(attr(a_c, ".vcov_type"), "clustered"))
})

test_that("unrecoverable variables named in subset or by give informative errors", {
  fitter <- function(formula, dat) lm(formula, data = dat)
  fit <- fitter(local(Math ~ logBLL + PIR), nhanes3lead)

  expect_error(adrf(fit, treat = "logBLL", subset = Race == "White"), "supply it to")
  expect_error(adrf(fit, treat = "logBLL", by = ~Race), "supply it to")

  # Supplying the data makes them available
  a <- adrf(fit, treat = "logBLL", n = 3, by = ~Race, data = nhanes3lead)
  expect_identical(nrow(attr(a, ".by_grid")), 4L)
})

test_that("models with offsets use the general path and give correct estimates", {
  set.seed(1)
  d <- transform(nhanes3lead, t = runif(nrow(nhanes3lead), 1, 3))
  d$cnt <- rpois(nrow(d), 2 * d$t)

  fits <- list(
    formula = glm(cnt ~ logBLL + PIR + offset(log(t)), data = d, family = poisson),
    argument = glm(cnt ~ logBLL + PIR, data = d, family = poisson, offset = log(t)),
    lm = lm(Math ~ logBLL + PIR + offset(PIR), data = d)
  )

  for (nm in names(fits)) {
    a <- adrf(fits[[nm]], treat = "logBLL", n = 3, range = c(0, 2))
    ref <- vapply(get_values(a), function(v) {
      mean(predict(fits[[nm]], newdata = transform(d, logBLL = v), type = "response"))
    }, numeric(1L))

    expect_equal(get_est(a), ref, tolerance = 1e-8, label = nm)
    expect_true(all(is.finite(get_se(a))), label = nm)
  }
})

test_that("models are routed to the fast path only when it reproduces their predictions", {
  route <- function(fit, type = "response", treat = "logBLL") {
    md <- adrftools:::process_model_data(fit)
    grid <- md[rep(1:20, 2), , drop = FALSE]
    grid[[treat]] <- rep(c(0, 1), each = 20)
    adrftools:::.get_model_type(fit, type = type, data_grid = grid, model_data = md)
  }

  d <- transform(nhanes3lead, t = runif(nrow(nhanes3lead), 1, 3))

  expect_identical(route(lm(Math ~ logBLL * Male + PIR, data = d)), "lm")
  expect_identical(route(glm(I(Block >= 12) ~ logBLL + PIR, data = d, family = binomial)), "glm")
  expect_identical(route(glm(I(Block >= 12) ~ logBLL + PIR, data = d, family = binomial),
                         type = "link"), "lm")
  # Offsets are added to the linear predictor, so these use the fast path
  expect_identical(route(glm(Block ~ logBLL + PIR + offset(log(t)), data = d,
                             family = poisson)), "glm")
  expect_identical(adrftools:::.get_model_type(lm(Math ~ logBLL, data = d), is_bayes = TRUE,
                                               data_grid = d), "default")

  skip_if_not_installed("betareg")
  expect_identical(route(betareg::betareg(I((Math + .5) / 31) ~ logBLL + PIR, data = d)),
                   "default")
})

test_that("fixest models are supported", {
  skip_if_not_installed("fixest")

  fit <- fixest::feols(Math ~ logBLL + PIR | Race, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 3)

  expect_equal(get_est(a), manual_adrf(fit, nhanes3lead, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))

  a_c <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional")
  expect_true(all(is.finite(get_se(a_c))))
})

test_that("a same-named dataset in the formula's environment is not used in place of the model's data", {
  # The formula is created here, next to an object named `d` that is not the
  # dataset the model is fit to
  f <- Math ~ poly(logBLL, 2) + PIR
  d <- nhanes3lead[1:200, ]

  fit_in <- function(dat) {
    d <- dat
    lm(f, data = d)
  }

  d_fit <- nhanes3lead[201:500, ]
  fit <- fit_in(d_fit)

  a <- adrf(fit, treat = "logBLL", n = 4, data = d_fit)
  expect_equal(get_est(a), manual_adrf(fit, d_fit, "logBLL", get_values(a)),
               tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))

  skip_if_not_installed("fwb")

  a_b <- adrf(fit, treat = "logBLL", n = 4, data = d_fit, vcov = "fwb",
              fwb.args = list(R = 20, verbose = FALSE))
  expect_equal(get_est(a_b), get_est(a), tolerance = 1e-8)
})

test_that("units with zero model weights do not change the unconditional variance", {
  d <- nhanes3lead
  set.seed(3)
  d$w <- runif(nrow(d))
  d$w[1:500] <- 0
  d2 <- d[d$w > 0, ]

  for (fitter in list(function(dd) lm(Math ~ logBLL * Age, data = dd, weights = w),
                      function(dd) suppressWarnings(glm(I(Math > 8) ~ logBLL * Age, data = dd,
                                                        weights = w, family = quasibinomial)))) {
    # glm() warns that zero-weight units are not used for the dispersion
    a_keep <- suppressWarnings(adrf(fitter(d), treat = "logBLL", n = 3, wts = "w", data = d))
    a_drop <- suppressWarnings(adrf(fitter(d2), treat = "logBLL", n = 3, wts = "w", data = d2))
    expect_equal(get_vcov(a_keep), get_vcov(a_drop), tolerance = 1e-8)
  }

  # The default target weights count every unit, including those with zero
  # model weights
  a <- adrf(lm(Math ~ logBLL * Age, data = d, weights = w), treat = "logBLL", n = 3)
  fit <- lm(Math ~ logBLL * Age, data = d, weights = w)
  expect_equal(get_est(a), manual_adrf(fit, d, "logBLL", get_values(a)), tolerance = 1e-8)
})

test_that("wts is validated", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  w <- rep(1, nrow(nhanes3lead))

  expect_error(adrf(fit, treat = "logBLL", wts = replace(w, 5, -2)), "negative")
  expect_error(adrf(fit, treat = "logBLL", wts = replace(w, 3, NA)), "missing or non-finite")
  expect_error(adrf(fit, treat = "logBLL", wts = 0 * w), "positive")
  expect_error(adrf(fit, treat = "logBLL", by = ~Male,
                    wts = ifelse(nhanes3lead$Male == 1, 0, 1)),
               "positive sum in every subgroup")
})

test_that("conditional variance with cluster uses sandwich::vcovCL() with processed clusters", {
  fit <- lm(Math ~ logBLL + Male + PIR, data = nhanes3lead)
  G <- function(a) cbind(1, get_values(a), mean(nhanes3lead$Male), mean(nhanes3lead$PIR))

  # A formula with a transformation was previously ignored
  a <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional",
            cluster = ~interaction(Race, Male))
  cl <- interaction(nhanes3lead$Race, nhanes3lead$Male)
  expect_equal(get_vcov(a), unname(G(a) %*% sandwich::vcovCL(fit, cluster = cl) %*% t(G(a))),
               tolerance = 1e-8)

  # A data frame previously gave an error
  a_df <- adrf(fit, treat = "logBLL", n = 3, vcov = "conditional", cluster = nhanes3lead["Race"])
  expect_equal(get_vcov(a_df),
               unname(G(a_df) %*% sandwich::vcovCL(fit, cluster = ~Race) %*% t(G(a_df))),
               tolerance = 1e-8)
})

test_that("clusters without units are dropped and single clusters give an error", {
  d <- nhanes3lead
  d$PIR[d$Race == levels(factor(d$Race))[1]] <- NA
  fit <- lm(Math ~ poly(logBLL, 2) + PIR, data = d)

  a <- adrf(fit, treat = "logBLL", n = 5, cluster = ~Race, data = d)
  expect_true(all(is.finite(get_se(a))))

  d$one <- 1
  expect_error(adrf(fit, treat = "logBLL", n = 5, cluster = ~one, data = d),
               "at least 2 clusters")
})

test_that("multiway cluster-robust variances are positive semidefinite", {
  set.seed(1)
  d <- nhanes3lead
  d$k <- sample(1:4, nrow(d), TRUE)
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = d)

  a <- adrf(fit, treat = "logBLL", n = 11, cluster = ~Race + k, data = d)
  ev <- eigen(get_vcov(a), symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-12 * max(ev))
})

test_that("character subsets refer to the original row names when rows were dropped", {
  d <- nhanes3lead
  d$PIR[1:30] <- NA
  fit <- lm(Math ~ logBLL + PIR, data = d)

  a <- adrf(fit, treat = "logBLL", n = 3, subset = as.character(31:40))
  expect_equal(get_est(a), manual_adrf(fit, d[31:40, ], "logBLL", get_values(a)),
               tolerance = 1e-8)

  expect_error(adrf(fit, treat = "logBLL", subset = as.character(1:5)), "not row names")
})

test_that("a supplied dataset is restricted to the rows used in the fit", {
  d <- nhanes3lead
  d$PIR[1:30] <- NA

  # The model's own data cannot be recovered
  fit <- lm(Math ~ poly(logBLL, 2) + PIR, data = d)
  fit$call$data <- quote(no_such_data)

  a <- adrf(fit, treat = "logBLL", n = 3, data = d)
  expect_equal(get_est(a), manual_adrf(fit, d[-(1:30), ], "logBLL", get_values(a)),
               tolerance = 1e-8)
})

test_that("models with offsets use the fast path and give correct estimates", {
  set.seed(4)
  d <- nhanes3lead
  d$cnt <- rpois(nrow(d), 3)
  d$ex <- runif(nrow(d), 1, 2)

  fits <- list(glm(cnt ~ poly(logBLL, 2) + Male + PIR, data = d, family = poisson,
                   offset = log(ex)),
               glm(cnt ~ poly(logBLL, 2) + Male + PIR + offset(log(ex)), data = d,
                   family = poisson))

  for (fit in fits) {
    expect_identical(adrftools:::.get_model_type(fit, data_grid = d[1:20, ], model_data = d),
                     "glm")

    a <- adrf(fit, treat = "logBLL", n = 4, vcov = "conditional")
    vals <- get_values(a)
    expect_equal(get_est(a), manual_adrf(fit, d, "logBLL", vals), tolerance = 1e-8)

    # Delta method with a numerical gradient of the average predictions, which
    # include the offset through predict()
    avg_pred <- function(b) {
      fit_b <- fit
      fit_b$coefficients <- b
      manual_adrf(fit_b, d, "logBLL", vals)
    }

    b <- coef(fit)
    G <- vapply(seq_along(b), function(j) {
      h <- 1e-6
      b1 <- b0 <- b
      b1[j] <- b1[j] + h
      b0[j] <- b0[j] - h
      (avg_pred(b1) - avg_pred(b0)) / (2 * h)
    }, numeric(length(vals)))

    expect_equal(get_vcov(a), unname(G %*% vcov(fit) %*% t(G)), tolerance = 1e-5)
  }
})

test_that("aliased coefficients are dropped for every variance", {
  d <- nhanes3lead
  d$PIR2 <- 2 * d$PIR

  fit_a <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR + PIR2,
               data = d, family = binomial)
  fit_b <- glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR,
               data = d, family = binomial)
  expect_true(anyNA(coef(fit_a)))

  for (v in c("conditional", "unconditional", "HC0")) {
    a_a <- suppressWarnings(adrf(fit_a, treat = "logBLL", n = 5, vcov = v))
    a_b <- adrf(fit_b, treat = "logBLL", n = 5, vcov = v)

    expect_equal(get_est(a_a), get_est(a_b), tolerance = 1e-8, label = v)
    expect_equal(get_vcov(a_a), get_vcov(a_b), tolerance = 1e-6, label = v)
  }

  # The same on the marginaleffects path
  local_mocked_bindings(.get_model_type = function(...) "default")

  a_a <- suppressWarnings(adrf(fit_a, treat = "logBLL", n = 5, vcov = "conditional"))
  a_b <- adrf(fit_b, treat = "logBLL", n = 5, vcov = "conditional")

  expect_equal(get_est(a_a), get_est(a_b), tolerance = 1e-8)
  expect_equal(get_vcov(a_a), get_vcov(a_b), tolerance = 1e-5)
})

test_that("models with a matrix response can be used", {
  d <- transform(nhanes3lead, s = Block, f = 20 - Block)
  fit <- glm(cbind(s, f) ~ poly(logBLL, 2) + Male + PIR, data = d, family = binomial)

  a <- adrf(fit, treat = "logBLL", n = 5, vcov = "conditional")
  vals <- get_values(a)

  expect_equal(get_est(a), manual_adrf(fit, d, "logBLL", vals), tolerance = 1e-8)

  mfx <- marginaleffects::avg_predictions(fit, variables = list(logBLL = vals),
                                          by = "logBLL")
  expect_equal(get_se(a), mfx$std.error, tolerance = 1e-4)

  expect_no_error(adrf(fit, treat = "logBLL", n = 5, data = d))
})

test_that("missing values in by, an empty subset, and a zero-width range give informative errors", {
  d <- nhanes3lead
  d$M2 <- d$Male
  d$M2[1:5] <- NA
  fit <- lm(Math ~ poly(logBLL, 2) + Male + PIR, data = d)

  expect_error(adrf(fit, treat = "logBLL", n = 5, by = ~M2), "missing values")
  expect_error(adrf(fit, treat = "logBLL", n = 5, subset = Male > 5), "at least one unit")
  expect_error(adrf(fit, treat = "logBLL", n = 5, subset = integer(0)), "at least one unit")
  expect_error(adrf(fit, treat = "logBLL", n = 5, range = c(1, 1)), "must differ")

  # Missing values outside the subset are fine
  expect_no_error(adrf(fit, treat = "logBLL", n = 5, by = ~M2, subset = !is.na(M2)))
})

test_that("models fit without data can be used", {
  logBLL <- nhanes3lead$logBLL
  Math <- nhanes3lead$Math
  Male2 <- nhanes3lead$Male
  Male2[1:5] <- NA
  other <- 1:3

  # glm() stores the environment of the formula as `data`, which here holds
  # objects of other lengths
  fit <- glm(Math ~ logBLL + Male2)
  expect_true(is.environment(fit$data))

  a <- expect_no_warning(adrf(fit, treat = "logBLL", n = 5, vcov = "conditional"))
  vals <- get_values(a)
  nd <- data.frame(logBLL, Male2)[!is.na(Male2), ]

  expect_equal(get_est(a), manual_adrf(fit, nd, "logBLL", vals), tolerance = 1e-8)
  expect_equal(get_vcov(a), get_vcov(adrf(glm(Math ~ logBLL + Male2, data = data.frame(Math, logBLL, Male2)),
                                          treat = "logBLL", n = 5, vcov = "conditional")),
               tolerance = 1e-10)

  # Inside a function
  f <- function() {
    x <- nhanes3lead$logBLL
    y <- nhanes3lead$Math
    z <- 1:2
    glm(y ~ x)
  }

  fit2 <- f()
  a2 <- adrf(fit2, treat = "x", n = 5)
  expect_equal(get_est(a2), manual_adrf(fit2, data.frame(x = nhanes3lead$logBLL), "x", get_values(a2)),
               tolerance = 1e-8)
})

test_that("treatments cannot share a name with the arguments or internals of the effect curve", {
  d <- transform(nhanes3lead, subset = logBLL, .est = logBLL)

  expect_error(adrf(lm(Math ~ subset + Male, data = d), treat = "subset", n = 5),
               "cannot be named")
  expect_error(adrf(lm(Math ~ .est + Male, data = d), treat = ".est", n = 5),
               "cannot be named")
})

test_that("character predictors with levels absent from subset use the fast path", {
  d <- transform(nhanes3lead, RaceC = as.character(Race))
  fit <- glm(I(Block >= 14) ~ poly(logBLL, 2) * RaceC + Male, data = d, family = binomial)

  s <- d$RaceC == "Black"
  md <- adrftools:::process_model_data(fit)
  grid <- md[s, , drop = FALSE][rep(1:20, 2), , drop = FALSE]
  expect_identical(adrftools:::.get_model_type(fit, data_grid = grid, model_data = md), "glm")

  a <- adrf(fit, treat = "logBLL", n = 5, subset = RaceC == "Black", vcov = "conditional")
  expect_equal(get_est(a), manual_adrf(fit, d[s, ], "logBLL", get_values(a)), tolerance = 1e-8)

  local_mocked_bindings(.get_model_type = function(...) "default")
  a_slow <- adrf(fit, treat = "logBLL", n = 5, subset = RaceC == "Black", vcov = "conditional")
  expect_equal(get_vcov(a), get_vcov(a_slow), tolerance = 1e-5)
})
