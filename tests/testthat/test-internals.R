test_that("get_kernel_w() returns normalized kernel weights", {
  v <- seq(0, 1, length.out = 11)
  x <- c(.05, .52, .93)

  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "kernel = {kernel}",
    {
      w <- adrftools:::get_kernel_w(x, v, kernel = kernel)

      expect_identical(dim(w), c(3L, 11L))
      expect_equal(rowSums(w), rep(1, 3), tolerance = 1e-10)
      expect_true(all(w >= 0))
    },
    kernel = c("gaussian", "epanechnikov", "triangular")
  )
})

test_that("get_kernel_w() handles grid points and extrapolation", {
  v <- seq(0, 1, length.out = 11)
  x <- c(.05, .52, .93)

  # Exact grid points get indicator weights
  expect_equal(adrftools:::get_kernel_w(v[c(2, 5)], v),
               rbind(replace(numeric(11), 2, 1), replace(numeric(11), 5, 1)))
  expect_equal(adrftools:::get_kernel_w(v, v), diag(11))

  # Points far outside the grid get zero weights under compact kernels
  w_far <- adrftools:::get_kernel_w(5, v, kernel = "epanechnikov")
  expect_true(all(w_far == 0))

  expect_error(adrftools:::get_kernel_w(x, v, kernel = "box"))
})

test_that("process_bw() handles grids, strings, functions, and numbers", {
  v <- seq(0, 1, length.out = 11)

  expect_equal(adrftools:::process_bw(NULL, v), .05)
  expect_equal(adrftools:::process_bw(NULL, v, constant = 1), .1)

  set.seed(1)
  v_irr <- sort(runif(20))
  expect_equal(adrftools:::process_bw(NULL, v_irr), adrftools:::bw.nrd(v_irr))

  # A mostly regular grid uses the modal spacing
  v_mostly <- c(seq(0, 1, by = .1), 1.05)
  expect_equal(adrftools:::process_bw(NULL, v_mostly), .05)

  expect_equal(adrftools:::process_bw("nrd0", v_irr), bw.nrd0(v_irr))
  expect_equal(adrftools:::process_bw(bw.SJ, v_irr), bw.SJ(v_irr))
  expect_equal(adrftools:::process_bw(.3, v), .3)

  expect_error(adrftools:::process_bw(-1, v))
  expect_error(adrftools:::process_bw(TRUE, v), "string, function, or number")
})

test_that("get_locpoly_w() reproduces polynomials up to degree p", {
  v <- seq(0, 1, length.out = 21)
  x <- c(.03, .51, .98)

  W <- adrftools:::get_locpoly_w(x, v)
  expect_equal(drop(W %*% v^3), x^3, tolerance = 1e-8)

  W0 <- adrftools:::get_locpoly_w(x, v, p = 0)
  expect_equal(W0, adrftools:::get_kernel_w(x, v))

  W1 <- adrftools:::get_locpoly_w(x, v, p = 1, kernel = "triangular")
  expect_equal(drop(W1 %*% (2 * v)), 2 * x, tolerance = 1e-8)

  expect_error(adrftools:::get_locpoly_w(5, v, kernel = "epanechnikov"),
               "extrapolation")
})

test_that(".gradient() computes numerical derivatives", {
  f <- function(b) c(sum(b^2), prod(b))
  b <- c(1, 2, 3)
  J <- rbind(2 * b, c(6, 3, 2))

  expect_equal(unname(adrftools:::.gradient(f, b)), J, tolerance = 1e-5)
  expect_equal(unname(adrftools:::.gradient(f, b, .direction = "left")), J, tolerance = 1e-5)
  expect_equal(unname(adrftools:::.gradient(f, b, .direction = "right")), J, tolerance = 1e-5)
  expect_equal(unname(adrftools:::.gradient(f, b, .method = "richardson")), J, tolerance = 1e-6)
  expect_equal(unname(adrftools:::.gradient(f, b, .method = "richardson", .direction = "left")),
               J, tolerance = 1e-5)
  expect_equal(unname(adrftools:::.gradient(f, b, .method = "richardson", .direction = "right")),
               J, tolerance = 1e-5)
})

test_that(".sim_p_value() agrees with exact p-values", {
  stat <- qchisq(.9, 3)

  set.seed(1)
  p <- adrftools:::.sim_p_value(stat, c(1, 1, 1), n = 1e5)
  expect_lt(abs(p - .1), 4 * sqrt(.1 * .9 / 1e5))

  # Splitting into chunks of unequal size gives the same kind of answer
  set.seed(1)
  p_chunk <- adrftools:::.sim_p_value(stat, c(1, 1, 1), n = 1e5 + 7, chunk = 3e4)
  expect_lt(abs(p_chunk - .1), 4 * sqrt(.1 * .9 / 1e5))

  # With finite df, Q / 3 has an F(3, 10) distribution
  set.seed(1)
  p_t <- adrftools:::.sim_p_value(stat, c(1, 1, 1), n = 1e5, df = 10)
  p_f <- pf(stat / 3, 3, 10, lower.tail = FALSE)
  expect_lt(abs(p_t - p_f), 4 * sqrt(p_f * (1 - p_f) / 1e5))

  # Unequal eigenvalues, against the saddlepoint approximation; zero
  # eigenvalues are dropped
  ev <- c(3, 1, .5, .1, 0)
  set.seed(1)
  p_w <- adrftools:::.sim_p_value(6, ev, n = 2e5)
  expect_equal(p_w, adrftools:::.quad_form_tail(6, ev, method = "saddlepoint"), tolerance = .03)

  expect_error(adrftools:::make_rmvt(c(0, 0), -diag(2)), "positive definite")
  expect_error(adrftools:::make_rmvt(c(0, 0), diag(3)), "Incompatible")
})

# test_that("simultaneous_conf_level() returns the pointwise level for a simultaneous band", {
#   V <- matrix(.5, 3, 3)
#   diag(V) <- 1
#
#   lvl <- adrftools:::simultaneous_conf_level(V, .95)
#   expect_gt(lvl, .95)
#   expect_lt(lvl, 1)
#
#   # Independent estimates follow the Sidak correction
#   set.seed(1)
#   expect_equal(adrftools:::simultaneous_conf_level(diag(3) * 4, .95), .95^(1 / 3),
#                tolerance = 1e-3)
#   expect_identical(adrftools:::simultaneous_conf_level(V, 0), 0)
#
#   expect_error(adrftools:::simultaneous_conf_level(matrix(1:4, 2), .95), "valid covariance")
# })

test_that(".to_psd() makes a matrix positive semidefinite", {
  m <- matrix(c(1, 1, 1, 1), 2)
  m_psd <- adrftools:::.to_psd(m)
  expect_true(all(eigen(m_psd, only.values = TRUE)$values > 0))

  pd <- diag(2)
  expect_identical(adrftools:::.to_psd(pd), pd)
})

test_that("pool_est() applies Rubin's rules", {
  est <- list(c(1, 2), c(3, 4))
  V <- list(diag(2), 3 * diag(2))

  p <- adrftools:::pool_est(est, V)
  expect_equal(p$est, c(2, 3))
  expect_equal(p$vcov, 2 * diag(2) + 1.5 * matrix(2, 2, 2))

  expect_null(adrftools:::pool_est(est)$vcov)
})

test_that("formatting utilities", {
  expect_identical(adrftools:::.format_p_value(c(.5, 1e-6, NA), digits = 3),
                   c("0.500", "< 0.001", NA))

  df <- data.frame(a = c(1.23456, -0.00001, NA), b = c("1", "2", "x"), c = c(Inf, 1, 2))
  out <- adrftools:::round_df_char(df, digits = 2, na_vals = ".")
  expect_identical(out$a, c("1.23", "-0.00", "."))
  expect_identical(out$c[1], "N/A")

  out_pad <- adrftools:::round_df_char(data.frame(a = c(1.5, 2.25)), digits = 2, pad = " ")
  expect_identical(nchar(out_pad$a), c(4L, 4L))

  expect_identical(adrftools:::round_df_char(data.frame(), 2), data.frame())
  expect_s3_class(adrftools:::round_df_char(matrix(1.234, 1, 1), 1), "data.frame")

  expect_identical(adrftools:::add_quotes("a"), '"a"')
  expect_identical(adrftools:::add_quotes("a", 1), "'a'")
  expect_identical(adrftools:::add_quotes("a", FALSE), "a")
  expect_identical(adrftools:::add_quotes("a", TRUE), '"a"')
  expect_identical(adrftools:::add_quotes("a", "`"), "`a`")
  expect_identical(adrftools:::add_quotes("a", "(["), "([a[(")
  expect_error(adrftools:::add_quotes("a", 3))

  expect_identical(adrftools:::center_just("ab", wrt = "abcdef"), "  ab")
  expect_type(adrftools:::center_just("ab"), "character")

  expect_identical(adrftools:::str2num(c("1", NA, "a")), c(1, NA, NA))
  expect_true(adrftools:::can_str2num(c("1", NA)))
  expect_false(adrftools:::can_str2num("a"))
  expect_true(adrftools:::can_str2num(1:2))

  expect_identical(adrftools:::get_by_grid_labels(data.frame(a = 1:2, b = c("x", "y"))),
                   c('a = 1, b = "x"', 'a = 2, b = "y"'))
})

test_that("miscellaneous helpers", {
  expect_identical(adrftools:::block_diag(matrix(1, 1, 2), 2),
                   rbind(c(1, 1, 0, 0), c(0, 0, 1, 1)))
  expect_identical(adrftools:::block_diag(1:2), matrix(1:2))

  expect_identical(adrftools:::get_varnames(quote(log(x) + y$z * f(w, 2))),
                   c("x", "y$z", "w"))

  df <- adrftools:::make_df(c("a", "b"), 2, types = c("integer", NA))
  expect_identical(names(df), c("a", "b"))
  expect_type(df$a, "integer")
  expect_identical(dim(adrftools:::make_df(NULL, NULL)), c(0L, 0L))
  expect_identical(rownames(adrftools:::make_df(1, c("r1", "r2"))), c("r1", "r2"))
  expect_error(adrftools:::make_df(2, 1, types = c("numeric", "numeric", "numeric")))
  expect_error(adrftools:::make_df(1, 1, types = "date"))

  expect_identical(adrftools:::make_list(2), list(NULL, NULL))
  expect_identical(adrftools:::make_list(NULL), list())
  expect_identical(names(adrftools:::make_list(c("a", "b"))), c("a", "b"))
  expect_error(adrftools:::make_list(list(1)))

  expect_identical(adrftools:::rep_with(1, c(a = 1, b = 2)), c(a = 1, b = 1))

  expect_true(adrftools:::all_the_same(c(1, 1, NA)))
  expect_false(adrftools:::all_the_same(c(1, NA), na.rm = FALSE))
  expect_true(adrftools:::all_the_same(c("a", "a")))
  expect_true(adrftools:::all_the_same(1))

  expect_identical(adrftools:::squish(c(-1, .5, 2), 0, 1), c(0, .5, 1))

  expect_true(adrftools:::is_count(3))
  expect_false(adrftools:::is_count(3.5))

  expect_true(adrftools:::any_apply(list(1, "a"), is.character))
  expect_false(adrftools:::any_apply(1:3, is.character))
  expect_true(adrftools:::all_apply(list(1, 2), is.numeric))
  expect_false(adrftools:::all_apply(factor("a"), is.numeric))

  f <- function(...) adrftools:::...get("b", 0)
  expect_identical(f(a = 1, b = 2), 2)
  expect_identical(f(a = 1), 0)

  g <- function(...) adrftools:::...mget(c("a", "c"))
  expect_identical(g(a = 1, b = 2), list(a = 1))
  expect_identical(g(b = 2), list())

  expect_s3_class(adrftools:::.it("a"), "cli_ansi_string")
})

test_that("process_df() returns Inf for chi-squared statistics and unsupported models", {
  fit <- glm(I(Block >= 14) ~ logBLL, data = nhanes3lead, family = binomial)
  expect_identical(adrftools:::process_df(fit), Inf)

  fit_lm <- lm(Math ~ logBLL, data = nhanes3lead)
  expect_equal(adrftools:::process_df(fit_lm), df.residual(fit_lm))

  expect_identical(adrftools:::process_df(structure(list(), class = "not_a_model")), Inf)
})

test_that("process_family() reduces families to their link functions", {
  fit <- glm(I(Block >= 14) ~ logBLL, data = nhanes3lead, family = binomial)

  f <- adrftools:::process_family(fit)
  expect_named(f, c("link", "linkfun", "linkinv", "mu.eta"), ignore.order = TRUE)

  expect_identical(adrftools:::process_family(fit, type = "link")$link, "identity")
  expect_identical(adrftools:::process_family(structure(list(), class = "not_a_model"))$link,
                   "identity")
})

test_that("process_ci.type() validates the interval type", {
  expect_null(adrftools:::process_ci.type("perc", "none"))
  expect_identical(adrftools:::process_ci.type("perc", "unconditional"), "wald")
  expect_identical(adrftools:::process_ci.type("perc", "posterior"), "perc")
  expect_error(adrftools:::process_ci.type("bc", "posterior"))
  expect_identical(adrftools:::process_ci.type("bc", "bootstrap", FALSE), "bc")
})

test_that("simultaneous posterior p-values invert the simultaneous intervals", {
  set.seed(1)
  S <- matrix(.7, 5, 5)
  diag(S) <- 1
  dr <- mvtnorm::rmvnorm(4000, mean = c(1.6, 1.8, 2, 2.1, 2.6) / 1.2, sigma = S)

  pp <- adrftools:::posterior_p_value(dr, parm = 1:5, null = 0)
  ps <- adrftools:::posterior_p_value(dr, parm = 1:5, null = 0, simultaneous = TRUE)

  expect_true(all(ps >= pp))

  # The smallest alpha at which the simultaneous band excludes 0, by search
  alphas <- seq(.0005, .5, by = .0005)
  lo <- vapply(alphas, function(a) {
    adrftools:::posterior_ci(dr, parm = 1:5, level = 1 - a, simultaneous = TRUE)[, 1]
  }, numeric(5L))
  inverted <- vapply(1:5, function(j) min(alphas[lo[j, ] > 0]), numeric(1L))

  expect_equal(ps, inverted, tolerance = .005)
})

test_that("posterior p-values invert the intervals and skip estimates with no variance", {
  # The pointwise p-value is the smallest alpha at which the type 7 percentile
  # interval excludes the null
  dr <- matrix(as.numeric(1:100), ncol = 1L)
  p <- adrftools:::posterior_p_value(dr, null = 3.2)
  ci_in <- adrftools:::posterior_ci(dr, level = 1 - p * 1.0001)
  ci_out <- adrftools:::posterior_ci(dr, level = 1 - p * .9999)
  expect_true(ci_out[1L] <= 3.2 && ci_in[1L] > 3.2)

  # A column with no variance gets NA and does not affect the others
  set.seed(1)
  D <- cbind(rnorm(200), rnorm(200, 1), 0)
  pp <- adrftools:::posterior_p_value(D, null = 0)
  ps <- adrftools:::posterior_p_value(D, null = 0, simultaneous = TRUE)

  expect_true(is.na(pp[3L]) && is.na(ps[3L]))
  expect_equal(ps[1:2], adrftools:::posterior_p_value(D[, 1:2], null = 0, simultaneous = TRUE))
  expect_true(all(pp[1:2] <= 1 & ps[1:2] >= pp[1:2]))
})

test_that("simultaneous posterior intervals follow Montiel Olea and Plagborg-Moller's Algorithm 2", {
  set.seed(1)
  S <- matrix(.95, 11, 11)
  diag(S) <- 1

  covers <- function(x, ci) {
    rowSums(x >= rep(ci[, 1], each = nrow(x)) & x <= rep(ci[, 2], each = nrow(x))) == ncol(x)
  }

  # The largest tail probability whose equal-tailed intervals contain at least
  # 90% of the draws, found by bisection as in the paper
  bisect <- function(dr, level) {
    covg <- function(z) {
      q <- apply(dr, 2L, quantile, probs = c(z, 1 - z), type = 7)
      mean(covers(dr, t(q)))
    }

    lo <- 0
    hi <- (1 - level) / 2
    for (i in 1:60) {
      mid <- (lo + hi) / 2
      if (covg(mid) >= level) lo <- mid else hi <- mid
    }

    t(apply(dr, 2L, quantile, probs = c(lo, 1 - lo), type = 7))
  }

  for (R in c(100, 4000)) {
    dr <- mvtnorm::rmvnorm(R, sigma = S)
    ci <- adrftools:::posterior_ci(dr, level = .9, simultaneous = TRUE)

    expect_identical(dim(ci), c(11L, 2L))
    expect_equal(unname(ci), unname(bisect(dr, .9)), tolerance = 1e-8)
    expect_gte(mean(covers(dr, ci)), .9)
  }

  ci_p <- adrftools:::posterior_ci(dr, level = .9)
  expect_true(all(ci[, 2] - ci[, 1] > ci_p[, 2] - ci_p[, 1]))
})

test_that("process_family() rebuilds link functions from the link name", {
  fam <- structure(list(family = "custom", link = "logit"), class = "family")
  local_mocked_bindings(get_family = function(...) fam, .package = "insight")

  f <- adrftools:::process_family(structure(list(), class = "some_model"))
  expect_equal(f$linkinv(0), .5)
  expect_equal(f$linkfun(.5), 0)
  expect_equal(f$mu.eta(0), .25)

  fam$link <- "not_a_link"
  f2 <- adrftools:::process_family(structure(list(), class = "some_model"))
  expect_identical(f2$link, "identity")
})

test_that("simultaneous Wald critical values and p-values are accurate and consistent", {
  # With independent estimates, the exact values follow from the Sidak formula
  k <- 5
  crit_exact <- qnorm((1 + .95^(1 / k)) / 2)

  set.seed(1)
  expect_equal(adrftools:::.mvt_crit(.95, diag(k)), crit_exact, tolerance = 1e-3,
               ignore_attr = TRUE)

  z <- c(-2.5, 1, 2.5, 3, 6)
  p_exact <- 1 - (1 - 2 * pnorm(-abs(z)))^k

  set.seed(1)
  p <- adrftools:::.mvt_p_values(z, diag(k))
  expect_equal(p, p_exact, tolerance = 1e-4)

  # Equal |z| give equal p-values, which decrease in |z|
  expect_identical(p[1L], p[3L])
  expect_true(all(diff(p[order(abs(z))]) <= 0))

  # Far from the null, the Bonferroni bound is used, which is within 1e-4
  expect_equal(p[5L], k * 2 * pnorm(-6))

  # One estimate gives the pointwise p-value
  expect_equal(adrftools:::.mvt_p_values(2, matrix(1), df = 10), 2 * pt(-2, 10))

  # The p-value at the critical value is 1 - conf_level
  R <- matrix(.6, 4, 4)
  diag(R) <- 1
  set.seed(2)
  cr <- adrftools:::.mvt_crit(.9, R, df = 30)
  set.seed(3)
  expect_equal(adrftools:::.mvt_p_values(cr, R, df = 30), .1, tolerance = 2e-3)
})

test_that("summary() intervals and plot() bands use the same simultaneous critical value", {
  fit <- lm(Math ~ poly(logBLL, 3) + Male + PIR, data = nhanes3lead)
  a <- adrf(fit, treat = "logBLL", n = 21)

  set.seed(10)
  s <- summary(a(logBLL = get_values(a)))
  set.seed(10)
  rb <- ribbon_data(plot(a))

  expect_equal(rb$ymin, s$conf.low, tolerance = 1e-12)
  expect_equal(rb$ymax, s$conf.high, tolerance = 1e-12)
})
