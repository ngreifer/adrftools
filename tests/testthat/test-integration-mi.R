skip_if_not_installed("mice")

# Applies `fitter` and adrf() to each completed dataset and pools the results
# with Rubin's rules, as an independent check of adrf() on a mira object
pool_by_hand <- function(imp, fitter, ...) {
  per_imp <- lapply(seq_len(imp$m), function(i) {
    di <- mice::complete(imp, i)
    adrf(fitter(di), data = di, ...)
  })

  ests <- do.call("rbind", lapply(per_imp, get_est))
  out <- list(est = colMeans(ests))

  if (is.null(get_vcov(per_imp[[1L]]))) {
    return(out)
  }

  W <- Reduce("+", lapply(per_imp, get_vcov)) / imp$m
  out$vcov <- W + (1 + 1 / imp$m) * cov(ests)

  out
}

make_imputations <- function(m = 3) {
  d <- nhanes3lead[, c("Math", "Block", "logBLL", "PIR", "Male", "Age")]
  set.seed(11)
  d$PIR[sample(nrow(d), 300)] <- NA

  suppressWarnings(mice::mice(d, m = m, printFlag = FALSE, seed = 12))
}

test_that("ADRFs from mira objects are pooled with Rubin's rules", {
  imp <- make_imputations()
  fits <- with(imp, lm(Math ~ poly(logBLL, 2) * Male + PIR))

  a <- adrf(fits, treat = "logBLL", n = 4, range = c(0, 2))

  expect_output(print(a), sprintf("sample size: %s in each imputed dataset",
                                  format(nrow(nhanes3lead), big.mark = ",")))
  expect_output(print(a), "imputed datasets: 3")

  # Pool single-dataset ADRFs by hand
  per_imp <- lapply(seq_len(imp$m), function(i) {
    di <- mice::complete(imp, i)
    fi <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = di)
    adrf(fi, treat = "logBLL", n = 4, range = c(0, 2))
  })

  ests <- do.call("rbind", lapply(per_imp, get_est))
  W <- Reduce("+", lapply(per_imp, get_vcov)) / imp$m
  B <- cov(ests)

  expect_equal(get_est(a), colMeans(ests), tolerance = 1e-8)
  expect_equal(get_vcov(a), W + (1 + 1 / imp$m) * B, tolerance = 1e-8)

  expect_identical(as.character(attr(a, ".vcov_type")), "unconditional")

  # Degrees of freedom for joint inference (Li, Raghunathan, and Rubin, 1991)
  # with the average relative increase in variance computed with a
  # pseudo-inverse of the rank-deficient within-imputation covariance
  m <- imp$m
  sv <- svd(W)
  keep <- sv$d > 1e-8 * sv$d[1L]
  k <- sum(keep)
  W_inv <- sv$v[, keep] %*% (t(sv$u[, keep]) / sv$d[keep])
  r <- (1 + 1 / m) * sum(diag(B %*% W_inv)) / k
  t_ <- k * (m - 1)
  df_li <- if (t_ > 4) 4 + (t_ - 4) * (1 + (1 - 2 / t_) / r)^2 else t_ * (1 + 1 / k) * (1 + 1 / r)^2 / 2

  expect_equal(adrftools:::.mi_df(W, B, m), df_li, tolerance = 1e-8)

  # The small-sample degrees of freedom (Reiter, 2007) are below both those of
  # the complete data and the large-sample degrees of freedom
  df_com <- min(vapply(fits$analyses, df.residual, numeric(1L)))
  expect_equal(attr(a, ".df"), adrftools:::.mi_df(W, B, m, df_com), tolerance = 1e-10)
  expect_lt(attr(a, ".df"), df_com)
  expect_lt(attr(a, ".df"), df_li)
})

test_that("mira ADRFs from glm fits use the large-sample degrees of freedom", {
  imp <- make_imputations(m = 5)
  fits <- with(imp, glm(I(Block >= 14) ~ poly(logBLL, 2) + Male + PIR, family = binomial))

  a <- adrf(fits, treat = "logBLL", n = 11)

  # The complete-data degrees of freedom are infinite, but the pooled ones are
  # finite because of the between-imputation variance
  expect_true(is.finite(attr(a, ".df")))
  expect_gt(attr(a, ".df"), 4)

  s <- summary(a(logBLL = 1))
  expect_equal(attr(attr(s, "crit"), "df"), attr(a, ".df"))

  # Without between-imputation variance, the degrees of freedom are those of
  # the complete data
  V <- diag(2)
  expect_identical(adrftools:::.mi_df(V, 0 * V, 5, df_com = Inf), Inf)
  expect_identical(adrftools:::.mi_df(V, 0 * V, 5, df_com = 100), 100)

  # With few imputations and estimates, Reiter's formula does not apply, and
  # the result is still below the complete-data degrees of freedom
  df_small <- adrftools:::.mi_df(V, .05 * V, 2, df_com = 100)
  expect_true(df_small > 0 && df_small < 100)
})

test_that("mira ADRFs work for glm, by, wts, and downstream functions", {
  imp <- make_imputations()
  fits <- with(imp, glm(I(Block >= 14) ~ logBLL * Male + PIR, family = binomial))

  a <- adrf(fits, treat = "logBLL", n = 11, by = ~Male)
  expect_length(get_est(a), 22L)

  ests <- vapply(seq_len(imp$m), function(i) {
    di <- mice::complete(imp, i)
    fi <- glm(I(Block >= 14) ~ logBLL * Male + PIR, data = di, family = binomial)
    manual_adrf(fi, di[di$Male == 1, ], "logBLL", get_values(a))
  }, numeric(11L))
  expect_equal(get_est(a)[12:22], rowMeans(ests), tolerance = 1e-8)

  expect_s3_class(summary(a(logBLL = c(0, 1))), "summary.curve_est")
  expect_s3_class(plot(curve_contrast(a)), "ggplot")
  expect_s3_class(summary(amef(a), method = "satterthwaite"), "summary.effect_curve")

  # Conditional variance with the default model covariance
  a_c <- adrf(fits, treat = "logBLL", n = 5, vcov = "conditional")
  expect_identical(as.character(attr(a_c, ".vcov_type")), "conditional")
  expect_true(all(is.finite(get_se(a_c))))

  # Numeric weights of length n or n * m
  set.seed(1)
  w <- rexp(nrow(imp$data))
  a_w <- adrf(fits, treat = "logBLL", n = 3, wts = w, range = c(0, 2))
  a_wm <- adrf(fits, treat = "logBLL", n = 3, wts = rep(w, imp$m), range = c(0, 2))
  expect_equal(get_est(a_w), get_est(a_wm), tolerance = 1e-10)
  expect_error(adrf(fits, treat = "logBLL", wts = w[-1]), "length")
})

test_that("mira ADRFs reject bootstrapping and mismatched families", {
  imp <- make_imputations()
  fits <- with(imp, lm(Math ~ logBLL + PIR))

  expect_error(adrf(fits, treat = "logBLL", vcov = "boot"), "multiply imputed")
  expect_error(adrf(fits, treat = "logBLL", vcov = "fwb"), "multiply imputed")
  expect_error(adrf(fits, treat = "logBLL", vcov = list(1, 2)), "length equal")

  fits$analyses[[2]] <- glm(I(Block >= 14) ~ logBLL + PIR, data = mice::complete(imp, 2),
                            family = binomial)
  expect_error(adrf(fits, treat = "logBLL"), "same family")
})

test_that("mimira objects from MatchThem are supported", {
  skip_if_not_installed("MatchThem")
  skip_if_not_installed("WeightIt")

  imp <- make_imputations(m = 2)
  wimp <- suppressMessages(
    MatchThem::weightthem(logBLL ~ PIR + Male + Age, datasets = imp,
                          approach = "within", method = "glm")
  )
  fits <- with(wimp, WeightIt::glm_weightit(Math ~ logBLL * Male + PIR))

  a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2))

  ests <- vapply(seq_len(2), function(i) {
    di <- mice::complete(imp, i)
    manual_adrf(fits$analyses[[i]], di, "logBLL", get_values(a))
  }, numeric(3L))

  expect_equal(get_est(a), rowMeans(ests), tolerance = 1e-8)
  expect_true(all(is.finite(get_se(a))))
})

test_that("mira ADRFs use the supplied vcov and cluster", {
  imp <- make_imputations(m = 2)
  fits <- with(imp, lm(Math ~ logBLL + PIR + Male))

  a_c <- adrf(fits, treat = "logBLL", n = 3, vcov = "conditional", range = c(0, 2))
  a_hc3 <- adrf(fits, treat = "logBLL", n = 3, vcov = "HC3", range = c(0, 2))
  a_cl <- adrf(fits, treat = "logBLL", n = 3, vcov = "conditional", cluster = ~Male,
               range = c(0, 2))

  # Pool single-dataset HC3 ADRFs by hand
  per_imp <- lapply(seq_len(2), function(i) {
    fi <- lm(Math ~ logBLL + PIR + Male, data = mice::complete(imp, i))
    adrf(fi, treat = "logBLL", n = 3, vcov = "HC3", range = c(0, 2))
  })
  ests <- do.call("rbind", lapply(per_imp, get_est))
  V <- Reduce("+", lapply(per_imp, get_vcov)) / 2 + 1.5 * cov(ests)

  expect_equal(get_vcov(a_hc3), V, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(get_vcov(a_c), get_vcov(a_hc3))))
  expect_false(isTRUE(all.equal(get_vcov(a_c), get_vcov(a_cl))))
  expect_true(attr(attr(a_cl, ".vcov_type"), "clustered"))
})

test_that("mira ADRFs can use variables outside the model in by, wts, and cluster", {
  imp <- make_imputations(m = 2)
  fits <- with(imp, lm(Math ~ logBLL + PIR))

  a <- adrf(fits, treat = "logBLL", n = 3, by = ~Male, range = c(0, 2))

  ests <- vapply(seq_len(2), function(i) {
    di <- mice::complete(imp, i)
    fi <- lm(Math ~ logBLL + PIR, data = di)
    manual_adrf(fi, di[di$Male == 1, ], "logBLL", get_values(a))
  }, numeric(3L))
  expect_equal(get_est(a)[4:6], rowMeans(ests), tolerance = 1e-8)

  a_c <- adrf(fits, treat = "logBLL", n = 3, cluster = ~Male, range = c(0, 2))
  expect_true(attr(attr(a_c, ".vcov_type"), "clustered"))
})

test_that("mira ADRFs on the link scale match pooling by hand", {
  imp <- make_imputations()
  f <- I(Block >= 12) ~ logBLL * Male + PIR
  fits <- with(imp, glm(I(Block >= 12) ~ logBLL * Male + PIR, family = binomial))

  a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), type = "link")
  hand <- pool_by_hand(imp, function(d) glm(f, data = d, family = binomial),
                       treat = "logBLL", n = 3, range = c(0, 2), type = "link")

  expect_equal(get_est(a), hand$est, tolerance = 1e-8)
  expect_equal(get_vcov(a), hand$vcov, tolerance = 1e-6)
  expect_true(all(get_est(a) < 0))
})

test_that("mira ADRFs with subset and by match pooling by hand", {
  imp <- make_imputations()
  f <- Math ~ poly(logBLL, 2) * Male + PIR
  fits <- with(imp, lm(Math ~ poly(logBLL, 2) * Male + PIR))

  a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), subset = Age > 8)
  hand <- pool_by_hand(imp, function(d) lm(f, data = d),
                       treat = "logBLL", n = 3, range = c(0, 2), subset = Age > 8)
  expect_equal(get_est(a), hand$est, tolerance = 1e-8)
  expect_equal(get_vcov(a), hand$vcov, tolerance = 1e-8)

  a_by <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), subset = Age > 8,
               by = ~Male)
  hand_by <- pool_by_hand(imp, function(d) lm(f, data = d),
                          treat = "logBLL", n = 3, range = c(0, 2), subset = Age > 8,
                          by = ~Male)
  expect_equal(get_est(a_by), hand_by$est, tolerance = 1e-8)
  expect_equal(get_vcov(a_by), hand_by$vcov, tolerance = 1e-8)
})

test_that("mira ADRFs with vcov = 'none' pool the estimates", {
  imp <- make_imputations()
  fits <- with(imp, lm(Math ~ logBLL + PIR))

  a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), vcov = "none")
  hand <- pool_by_hand(imp, function(d) lm(Math ~ logBLL + PIR, data = d),
                       treat = "logBLL", n = 3, range = c(0, 2), vcov = "none")

  expect_equal(get_est(a), hand$est, tolerance = 1e-8)
  expect_null(get_vcov(a))
  expect_identical(as.character(attr(a, ".vcov_type")), "none")
})

test_that("mira ADRFs from models without a fast path match pooling by hand", {
  skip_if_not_installed("betareg")

  d <- nhanes3lead[, c("Math", "logBLL", "PIR", "Male")]
  d$y01 <- (d$Math + .5) / 31
  set.seed(11)
  d$PIR[sample(nrow(d), 300)] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 2, printFlag = FALSE, seed = 12))

  fits <- with(imp, betareg::betareg(y01 ~ logBLL + PIR + Male))

  for (vc in c("unconditional", "conditional")) {
    a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), vcov = vc)
    hand <- pool_by_hand(imp, function(dd) betareg::betareg(y01 ~ logBLL + PIR + Male, data = dd),
                         treat = "logBLL", n = 3, range = c(0, 2), vcov = vc)

    expect_equal(get_est(a), hand$est, tolerance = 1e-8, label = vc)
    expect_equal(get_vcov(a), hand$vcov, tolerance = 1e-5, label = vc)
  }
})

test_that("mira ADRFs from Bayesian models stack the posterior draws", {
  skip_on_cran()
  skip_if_not_installed("rstanarm")

  d <- nhanes3lead[1:500, c("Math", "logBLL", "PIR", "Male")]
  set.seed(1)
  d$PIR[sample(nrow(d), 50)] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 2, printFlag = FALSE, seed = 1))

  fitter <- function(dd) {
    suppressWarnings(rstanarm::stan_glm(Math ~ logBLL * Male + PIR, data = dd, chains = 1,
                                        iter = 300, refresh = 0, seed = 1))
  }
  fits <- list(analyses = lapply(seq_len(2), function(i) fitter(mice::complete(imp, i))))
  class(fits) <- c("mira", "matrix")

  a <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), by = ~Male)
  expect_identical(dim(attr(a, ".draws")), c(300L, 6L))

  per_imp <- lapply(fits$analyses, function(fi) {
    attr(adrf(fi, treat = "logBLL", n = 3, range = c(0, 2), by = ~Male), ".draws")
  })
  draws <- do.call("rbind", per_imp)

  expect_equal(unname(attr(a, ".draws")), unname(draws), tolerance = 1e-8)
  expect_equal(get_est(a), unname(colMeans(draws)), tolerance = 1e-8)

  a_none <- adrf(fits, treat = "logBLL", n = 3, range = c(0, 2), vcov = "none")
  expect_null(attr(a_none, ".draws"))
})
