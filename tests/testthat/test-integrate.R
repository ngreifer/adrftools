# Weighting the grid by the distribution of the treatment (`integrate`)

# Density of the treatment computed as adrf() should, from all units used and
# with the target weights, evaluated at `values`. adrf() interpolates from a
# 512-point grid, so its weights agree with these to about 1e-3.
manual_density <- function(A, w, values) {
  dens <- stats::density(A, weights = w / sum(w), bw = stats::bw.nrd0(A), n = 4096,
                         from = min(A) - 3 * stats::bw.nrd0(A),
                         to = max(A) + 3 * stats::bw.nrd0(A))
  stats::approx(dens$x, dens$y, xout = values)$y
}

manual_grid_w <- function(A, w, values) {
  wt <- adrftools:::get_trapezoidal_w(values) * manual_density(A, w, values)
  wt / sum(wt)
}

test_that("integrate = 'uniform' is the default and leaves results unchanged", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  expect_identical(summary(a, "flat", method = "satterthwaite"),
                   summary(a, "flat", method = "satterthwaite", integrate = "uniform"))
  expect_identical(coef(curve_projection(a, "linear")),
                   coef(curve_projection(a, "linear", integrate = "uniform")))
})

test_that("integrate = 'empirical' weights the grid by the density of all treatment values", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)
  v <- get_values(a)

  # The density uses every unit, not only those within the range of the grid
  expect_true(any(nhanes3lead$logBLL > max(v)))
  w <- adrftools:::get_grid_w(v, "empirical", attr(a, ".treat_density"))[[1L]]
  expect_equal(w, manual_grid_w(nhanes3lead$logBLL, rep(1, nrow(nhanes3lead)), v),
               tolerance = 5e-3)

  # The flat test is the Satterthwaite approximation to the weighted quadratic
  # form, computed by hand
  est <- get_est(a)
  V <- get_vcov(a)
  C <- diag(length(v)) - matrix(w, length(v), length(v), byrow = TRUE)
  stat <- sum(w * drop(C %*% est)^2)
  ev <- eigen(sqrt(w) * t(sqrt(w) * (C %*% V %*% t(C))), symmetric = TRUE,
              only.values = TRUE)$values
  df <- df.residual(fit)
  p_hand <- pf(stat / sum(ev), sum(ev)^2 / sum(ev^2), df, lower.tail = FALSE)

  s <- summary(a, "flat", method = "satterthwaite", integrate = "empirical")
  expect_equal(s$p.value, p_hand, tolerance = 1e-6)
  expect_identical(attr(s, "integrate"), "empirical")
  expect_output(print(s), "weighted by the distribution of logBLL")

  # The empirical and uniform tests differ (compared with a test whose
  # p-values are not near the lower bound of 1e-12)
  expect_false(isTRUE(all.equal(
    summary(a, "linear", method = "satterthwaite", integrate = "empirical")$p.value,
    summary(a, "linear", method = "satterthwaite")$p.value)))
})

test_that("target weights and subgroups define the density", {
  set.seed(5)
  wts <- rexp(nrow(nhanes3lead))
  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 11, by = ~Male, wts = wts)
  v <- get_values(a)
  td <- attr(a, ".treat_density")

  for (k in 1:2) {
    i <- nhanes3lead$Male == (k - 1)
    w_k <- adrftools:::get_grid_w(v, "empirical", td, k)[[1L]]
    expect_equal(w_k, manual_grid_w(nhanes3lead$logBLL[i], wts[i], v), tolerance = 5e-3)
  }

  # A contrast uses the pooled density of its two subgroups, each with its own
  # bandwidth, weighted by their total target weights
  cc <- curve_contrast(a)
  w_c <- adrftools:::get_grid_w(v, "empirical", attr(cc, ".treat_density"), 1L)[[1L]]
  f_pool <- Reduce("+", lapply(0:1, function(k) {
    i <- nhanes3lead$Male == k
    sum(wts[i]) * manual_density(nhanes3lead$logBLL[i], wts[i], v)
  })) / sum(wts)
  w_pool <- adrftools:::get_trapezoidal_w(v) * f_pool
  expect_equal(w_c, w_pool / sum(w_pool), tolerance = 5e-3)

  # Each subgroup is tested with its own weights
  s <- summary(a, "linear", method = "satterthwaite", integrate = "empirical")
  s1 <- summary(a, "linear", method = "satterthwaite", integrate = "empirical",
                subset = Male == 1)
  expect_equal(s1$p.value, s$p.value[2L], tolerance = 1e-10)

  expect_no_error(summary(cc, method = "satterthwaite", integrate = "empirical"))
})

test_that("an empirically weighted projection is weighted least squares on the grid", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)
  v <- get_values(a)

  p <- curve_projection(a, "linear", integrate = "empirical")
  w <- adrftools:::get_grid_w(v, "empirical", attr(a, ".treat_density"))[[1L]]

  expect_equal(unname(coef(p)), unname(coef(lm(get_est(a) ~ v, weights = w))),
               tolerance = 1e-8)
  expect_false(isTRUE(all.equal(coef(p), coef(curve_projection(a, "linear")))))
  expect_output(print(p), "weighted by the distribution of logBLL")

  # anova() requires the same weighting
  expect_error(anova(curve_projection(a, "quadratic"), p), "integrate")
  expect_no_error(anova(curve_projection(a, "quadratic", integrate = "empirical"), p))
})

test_that("modified curves keep the density, and objects without it give an error", {
  fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  # A reference value between grid points adds a grid point, whose density is
  # interpolated
  r <- reference_curve(amef(a), reference = 1.234)
  expect_identical(attr(r, ".treat_density"), attr(a, ".treat_density"))
  expect_true(is.finite(summary(r, method = "satterthwaite", integrate = "empirical")$p.value))

  a_old <- a
  attr(a_old, ".treat_density") <- NULL
  expect_error(summary(a_old, "flat", integrate = "empirical"), "rerun")
  expect_error(curve_projection(a_old, "linear", integrate = "empirical"), "rerun")
  expect_error(summary(a, "flat", integrate = "density"), "integrate")
})

test_that("multiply imputed data pool the densities across imputations", {
  skip_if_not_installed("mice")

  d <- nhanes3lead[, c("Math", "logBLL", "PIR", "Male")]
  set.seed(11)
  d$logBLL[sample(nrow(d), 200)] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 3, printFlag = FALSE, seed = 12))
  fits <- with(imp, lm(Math ~ poly(logBLL, 2) * Male + PIR))

  a <- adrf(fits, treat = "logBLL", n = 11)
  v <- get_values(a)

  A <- unlist(lapply(1:3, function(i) mice::complete(imp, i)$logBLL))
  w <- adrftools:::get_grid_w(v, "empirical", attr(a, ".treat_density"))[[1L]]

  # The bandwidth is that of the stacked treatment values
  expect_equal(w, manual_grid_w(A, rep(1, length(A)), v), tolerance = 5e-3)
  expect_true(is.finite(summary(a, "flat", method = "satterthwaite",
                                integrate = "empirical")$p.value))
})
