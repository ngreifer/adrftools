# Models of categorical outcomes give one prediction per category, and adrf()
# gives a curve for each category, stored as a column of `.by_grid`

d_cat <- transform(nhanes3lead,
                   M3 = cut(Math, 3, labels = c("lo", "mid", "hi")))

categories <- c("lo", "mid", "hi")

fit_categorical <- function(model, data = d_cat) {
  switch(model,
         polr = MASS::polr(ordered(M3) ~ logBLL + PIR + Male, data = data,
                           Hess = TRUE),
         multinom_weightit = WeightIt::multinom_weightit(M3 ~ logBLL + PIR + Male,
                                                         data = data),
         ordinal_weightit = WeightIt::ordinal_weightit(ordered(M3) ~ logBLL + PIR + Male,
                                                       data = data))
}

# Central-difference Jacobian of `f` at `b`
jacobian_fd <- function(f, b, eps = 1e-6) {
  vapply(seq_along(b), function(j) {
    h <- eps * max(abs(b[j]), 1)
    bp <- bm <- b
    bp[j] <- b[j] + h / 2
    bm[j] <- b[j] - h / 2
    (f(bp) - f(bm)) / h
  }, numeric(length(f(b))))
}

test_that("each category of the outcome gets an ADRF", {
  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "{model}",
    {
      skip_if_not_installed(pkg)

      fit <- fit_categorical(model)

      # No `type`: polr() models give probabilities, not "response"
      a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                                 vcov = stats::vcov(fit)))

      vals <- get_values(a)

      expect_identical(attr(a, ".outcome"), "outcome")
      expect_identical(as.character(attr(a, ".by_grid")[["outcome"]]), categories)

      ref <- marginaleffects::avg_predictions(fit, variables = list(logBLL = vals),
                                              by = c("group", "logBLL"),
                                              vcov = stats::vcov(fit))
      ref <- ref[order(match(ref$group, categories), ref$logBLL), ]

      expect_equal(get_est(a), ref$estimate, tolerance = 1e-8)
      expect_equal(get_se(a), ref$std.error, tolerance = 1e-5)

      # The probabilities sum to 1 at every treatment value
      expect_equal(rowSums(matrix(get_est(a), nrow = 3L)), rep(1, 3L),
                   tolerance = 1e-8)

      est <- a(logBLL = vals)
      expect_identical(as.character(est$outcome), rep(categories, each = 3L))
    },
    patrick::cases(
      list(model = "polr", pkg = "MASS"),
      list(model = "multinom_weightit", pkg = "WeightIt"),
      list(model = "ordinal_weightit", pkg = "WeightIt")
    )
  )
})

test_that("the unconditional variance of category ADRFs matches a hand computation", {
  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "{model}",
    {
      skip_if_not_installed(pkg)

      fit <- fit_categorical(model)

      a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                                 by = ~Male))
      vals <- get_values(a)

      # Estimates ordered by category, then subgroup, then treatment value
      grid <- expand.grid(v = vals, Male = 0:1, cat = categories,
                          stringsAsFactors = FALSE)

      preds_at <- function(b) {
        m <- marginaleffects::set_coef(fit, b)

        lapply(vals, function(v) {
          p <- predict(m, newdata = transform(d_cat, logBLL = v), type = pred_type)
          p[, categories]
        })
      }

      est_fun <- function(b) {
        P <- preds_at(b)

        vapply(seq_len(nrow(grid)), function(j) {
          mean(P[[match(grid$v[j], vals)]][d_cat$Male == grid$Male[j], grid$cat[j]])
        }, numeric(1L))
      }

      b0 <- marginaleffects::get_coef(fit)
      est <- est_fun(b0)

      expect_equal(get_est(a), unname(est), tolerance = 1e-8)

      G <- jacobian_fd(est_fun, b0)
      P0 <- preds_at(b0)
      N <- nrow(d_cat)

      # Matched to the coefficients by name: for polr() models, get_coef()
      # puts the thresholds first and estfun() puts them last
      infl <- suppressMessages(adrftools:::.coef_influence(fit))
      infl <- infl[, names(b0), drop = FALSE]

      S <- infl %*% t(G)

      for (j in seq_len(nrow(grid))) {
        g <- d_cat$Male == grid$Male[j]
        S[g, j] <- S[g, j] +
          (P0[[match(grid$v[j], vals)]][g, grid$cat[j]] - est[j]) * N / sum(g)
      }

      expect_equal(get_vcov(a), unname(crossprod(S) / N^2), tolerance = 1e-5)
    },
    patrick::cases(
      list(model = "polr", pkg = "MASS", pred_type = "probs"),
      list(model = "multinom_weightit", pkg = "WeightIt", pred_type = "response")
    )
  )
})

test_that("`outcome` selects categories of the outcome", {
  skip_if_not_installed("WeightIt")

  fit <- fit_categorical("ordinal_weightit")

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2)))
  vals <- get_values(a)

  # In adrf(), only the selected categories are computed, in the order of the
  # categories
  a_sel <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                                 outcome = c("hi", "lo")))

  expect_identical(as.character(attr(a_sel, ".by_grid")[["outcome"]]), c("lo", "hi"))
  expect_equal(get_est(a_sel), get_est(a)[c(1:3, 7:9)], tolerance = 1e-10)
  expect_equal(get_vcov(a_sel), get_vcov(a)[c(1:3, 7:9), c(1:3, 7:9)],
               tolerance = 1e-10)

  # In the effect curve and the functions that take it
  est_hi <- a(logBLL = vals, outcome = "hi")
  expect_identical(as.character(est_hi$outcome), rep("hi", 3L))
  expect_equal(est_hi$estimate, get_est(a)[7:9], tolerance = 1e-10)

  s <- summary(a)
  s_mid <- summary(a, outcome = "mid")
  expect_identical(nrow(s_mid), 1L)
  expect_equal(s_mid$p.value, s$p.value[2L], tolerance = 1e-10)

  p <- curve_projection(a, "linear")
  expect_identical(nrow(summary(p, outcome = "lo")), 2L)

  expect_s3_class(plot(a, outcome = c("lo", "mid")), "ggplot")

  r <- reference_curve(a, reference = 0)
  expect_equal(r(logBLL = 1, outcome = "hi")$estimate,
               get_est(a)[8L] - get_est(a)[7L],
               tolerance = 1e-10)

  # Errors
  expect_error(suppressMessages(adrf(fit, treat = "logBLL", n = 3, outcome = "top")),
               "categories of the outcome")
  expect_error(a(logBLL = 1, outcome = "top"), "categories of the outcome")
  expect_error(a(logBLL = 1, subset = Male == 1), "unused argument")

  fit_lm <- lm(Math ~ logBLL + PIR, data = nhanes3lead)
  expect_error(adrf(fit_lm, treat = "logBLL", n = 3, outcome = "lo"),
               "separate prediction for each category")
  expect_error(summary(adrf(fit_lm, treat = "logBLL", n = 3), outcome = "lo"),
               "separate prediction for each category")
})

test_that("categories of the outcome do not clash with a `by` variable named `outcome`", {
  skip_if_not_installed("MASS")

  d <- transform(d_cat, outcome = Male)
  fit <- fit_categorical("polr", data = d)

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                             by = ~outcome))

  expect_identical(attr(a, ".outcome"), ".outcome")
  expect_named(attr(a, ".by_grid"), c(".outcome", "outcome"))

  est <- a(logBLL = 1, outcome = "hi", subset = outcome == 1)

  expect_identical(nrow(est), 1L)
  expect_identical(as.character(est$.outcome), "hi")
  expect_identical(est$outcome, 1L)

  expect_output(print(a), "by: outcome")
  expect_output(print(a), "outcome categories: lo, mid, hi")
})

test_that("curve contrasts compare subgroups within each category", {
  skip_if_not_installed("MASS")

  fit <- fit_categorical("polr")

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                             by = ~Male))
  cc <- curve_contrast(a)

  # Estimates ordered by category, then subgroup (Male), then treatment value
  est <- matrix(get_est(a), nrow = 3L)
  expect_equal(get_est(cc), as.vector(est[, c(2, 4, 6)] - est[, c(1, 3, 5)]),
               tolerance = 1e-10)

  expect_identical(attr(cc, ".contrast"), rep("[Male = 1] - [Male = 0]", 3L))

  est_mid <- cc(logBLL = get_values(a), outcome = "mid")
  expect_identical(as.character(est_mid$outcome), rep("mid", 3L))
  expect_equal(est_mid$estimate, get_est(cc)[4:6], tolerance = 1e-10)

  expect_identical(nrow(summary(cc, outcome = "lo")), 1L)
  expect_s3_class(plot(cc), "ggplot")

  # Without `by`, the categories are not contrasted with each other
  a1 <- suppressMessages(adrf(fit, treat = "logBLL", n = 3))
  expect_error(curve_contrast(a1), "must have been called with `by`")
})

test_that("plots with subgroups give each category of the outcome its own panel", {
  skip_if_not_installed("MASS")

  fit <- fit_categorical("polr")

  a1 <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2)))
  a2 <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                              by = ~Male))

  # Without `by`, colors distinguish the categories in one panel
  p1 <- plot(a1)
  expect_s3_class(p1$facet, "FacetNull")
  expect_identical(levels(p1$data$.by_id), categories)

  # With `by`, colors distinguish the subgroups and panels the categories
  p2 <- plot(a2)
  expect_s3_class(p2$facet, "FacetWrap")
  expect_identical(levels(p2$data$.by_id), c("Male = 0", "Male = 1"))
  expect_identical(levels(p2$data$.facet), categories)

  # A subgroup keeps its color when the other is left out
  p3 <- plot(a2, outcome = "lo", subset = Male == 1)
  expect_identical(unique(as.character(p3$data$.facet)), "lo")
  expect_identical(unique(ggplot2::layer_data(p3)$colour),
                   unique(ggplot2::layer_data(p2)$colour[ggplot2::layer_data(p2)$group == 2L]))

  p4 <- plot(curve_contrast(a2))
  expect_s3_class(p4$facet, "FacetWrap")
  expect_identical(levels(p4$data$contrast), "[Male = 1] - [Male = 0]")
})

test_that("`type` can request a single prediction per unit", {
  skip_if_not_installed("WeightIt")

  fit <- fit_categorical("ordinal_weightit")
  scores <- c(lo = 1, mid = 2, hi = 3)

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2)))
  a_mean <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                                  type = "mean", values = scores))

  expect_null(attr(a_mean, ".outcome"))
  expect_null(attr(a_mean, ".by_grid"))

  # The mean score is the sum of the scores times the category probabilities
  expect_equal(get_est(a_mean),
               drop(matrix(get_est(a), nrow = 3L) %*% scores),
               tolerance = 1e-8)
})

test_that("shift curves and the ADE are computed for each category", {
  skip_if_not_installed("WeightIt")

  fit <- fit_categorical("multinom_weightit")

  sc <- suppressMessages(shift_curve(fit, treat = "logBLL", n = 5))
  ad <- ade(sc)

  expect_identical(as.character(ad$outcome), categories)

  # The derivatives of the probabilities sum to 0
  expect_equal(sum(ad$estimate), 0, tolerance = 1e-6)

  expect_equal(ade(sc, outcome = "mid")$estimate, ad$estimate[2L], tolerance = 1e-10)
})

test_that("bootstrapped category ADRFs are ordered like the estimates", {
  skip_if_not_installed("fwb")
  skip_if_not_installed("WeightIt")

  fit <- fit_categorical("multinom_weightit")

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                             vcov = "none"))

  set.seed(123)
  a_b <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2),
                               vcov = "fwb", fwb.args = list(R = 20, verbose = FALSE)))

  expect_equal(get_est(a_b), get_est(a), tolerance = 1e-8)
  expect_identical(dim(attr(a_b, ".boot")$t), c(20L, 9L))
  expect_equal(rowSums(matrix(attr(a_b, ".boot")$t[1L, ], nrow = 3L)), rep(1, 3L),
               tolerance = 1e-8)
})

test_that("category ADRFs from multiply imputed data are pooled", {
  skip_if_not_installed("mice")
  skip_if_not_installed("MASS")

  d <- d_cat[c("M3", "logBLL", "PIR", "Male")]
  set.seed(11)
  d$PIR[sample(nrow(d), 300L)] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 2, printFlag = FALSE, seed = 12))

  fits <- with(imp, MASS::polr(ordered(M3) ~ logBLL + PIR + Male, Hess = TRUE))

  a <- suppressMessages(suppressWarnings(adrf(fits, treat = "logBLL", n = 3,
                                              range = c(0, 2))))

  per_imp <- lapply(seq_len(imp$m), function(i) {
    di <- mice::complete(imp, i)
    fit_i <- MASS::polr(ordered(M3) ~ logBLL + PIR + Male, data = di, Hess = TRUE)
    suppressMessages(adrf(fit_i, treat = "logBLL", n = 3, range = c(0, 2), data = di))
  })

  expect_identical(as.character(attr(a, ".by_grid")[["outcome"]]), categories)
  expect_equal(get_est(a), colMeans(do.call("rbind", lapply(per_imp, get_est))),
               tolerance = 1e-8)
})

test_that("brms categorical models give an ADRF for each category", {
  skip_if_not_installed("brms")
  skip_if_not(identical(Sys.getenv("ADRFTOOLS_TEST_BRMS"), "true"),
              "Set ADRFTOOLS_TEST_BRMS=true to compile a brms model")

  fit <- suppressMessages(suppressWarnings(
    brms::brm(M3 ~ logBLL + Male + PIR, data = d_cat, family = brms::categorical(),
              chains = 1, iter = 500, refresh = 0, seed = 1, silent = 2)
  ))

  a <- adrf(fit, treat = "logBLL", n = 3, range = c(0, 2))
  vals <- get_values(a)

  expect_identical(as.character(attr(a, ".by_grid")[["outcome"]]), categories)

  # posterior_epred() gives draws x units x categories
  draws <- lapply(categories, function(k) {
    vapply(vals, function(v) {
      ep <- brms::posterior_epred(fit, newdata = transform(d_cat, logBLL = v))
      rowMeans(ep[, , k])
    }, numeric(250L))
  }) |>
    do.call(what = "cbind")

  expect_equal(unname(attr(a, ".draws")), unname(draws), tolerance = 1e-8)
  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)
})

test_that("print() and as.data.frame() treat the categories of the outcome apart from subgroups", {
  skip_if_not_installed("MASS")

  fit <- fit_categorical("polr")

  a <- suppressMessages(adrf(fit, treat = "logBLL", n = 3, range = c(0, 2), by = ~Male))

  # The categories are not counted as subgroups
  expect_output(print(a), "by: Male \\(2 subgroups\\)")

  d <- as.data.frame(a, simultaneous = FALSE)
  expect_named(d, c("outcome", "Male", "logBLL", "estimate", "conf.low", "conf.high"))
  expect_identical(nrow(d), 3L * 2L * 3L)
  expect_identical(as.character(unique(d$outcome)), categories)
  expect_equal(d$estimate, get_est(a), tolerance = 1e-12)

  d_hi <- as.data.frame(a, simultaneous = FALSE, outcome = "hi")
  expect_identical(as.character(unique(d_hi$outcome)), "hi")
  expect_equal(d_hi$estimate, get_est(a)[13:18], tolerance = 1e-12)
})
