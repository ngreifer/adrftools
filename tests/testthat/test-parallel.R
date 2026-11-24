skip_if_not_installed("future")
skip_if_not_installed("future.apply")

# These tests use `multisession` workers, which load the installed build of the
# package rather than the one loaded with `load_all()`, so install the working
# copy first (e.g., to a library put first on `.libPaths()`).

# Sets a plan with two workers, returning whether that worked here (a single
# core or a sandbox makes it fail)
start_two_workers <- function() {
  tryCatch({
    future::plan(future::multisession, workers = 2L)
    isTRUE(future::nbrOfWorkers() >= 2L)
  }, error = function(e) FALSE)
}

# The estimates of an effect curve and whatever its inference rests on, with
# the next draws from the random number stream after it was computed
run_curve <- function(make) {
  a <- make()

  list(est = get_est(a),
       vcov = get_vcov(a),
       boot = unname(attr(a, ".boot")$t),
       draws = unname(attr(a, ".draws")),
       after = stats::runif(3L))
}

test_that("effect curves do not depend on the plan", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  patrick::with_parameters_test_that(
    "{case}",
    {
      skip_if_not_installed(pkg)

      old <- future::plan(future::sequential)
      on.exit(future::plan(old), add = TRUE)

      set.seed(123)
      one <- run_curve(make)

      skip_if_not(start_two_workers(), "no second worker available")

      set.seed(123)
      many <- run_curve(make)

      expect_equal(many[names(many) != "after"], one[names(one) != "after"])

      # The code after the call must not depend on the plan either
      expect_identical(many[["after"]], one[["after"]])
    },
    .cases = patrick::cases(
      list(case = "the columns of the numerical gradient", pkg = "betareg",
           make = function() {
             d <- transform(nhanes3lead, y = (Block + .5) / (max(Block) + 1))
             fit <- betareg::betareg(y ~ splines::ns(logBLL, 3) * Male + Age,
                                     data = d)
             adrf(fit, treat = "logBLL", n = 11)
           }),
      list(case = "the gradient for each category of the outcome", pkg = "MASS",
           make = function() {
             d <- transform(nhanes3lead, y = cut(Block, 3))
             fit <- MASS::polr(y ~ splines::ns(logBLL, 3) + Male + Age,
                               data = d, Hess = TRUE)
             suppressMessages(adrf(fit, treat = "logBLL", n = 11))
           }),
      list(case = "the fractional weighted bootstrap", pkg = "fwb",
           make = function() {
             fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
             adrf(fit, treat = "logBLL", n = 5, vcov = "fwb",
                  fwb.args = list(R = 20, verbose = FALSE))
           }),
      list(case = "the traditional bootstrap", pkg = "fwb",
           make = function() {
             fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
             adrf(fit, treat = "logBLL", n = 5, vcov = "boot",
                  fwb.args = list(R = 20, verbose = FALSE))
           }),
      list(case = "the posterior draws at each treatment value", pkg = "bartisan",
           make = function() {
             fit <- bartisan::bartisan(
               Math ~ logBLL + Male + PIR + Age, data = nhanes3lead[1:600, ],
               family = gaussian(),
               control = bartisan::bartisan_control(num_burn = 20,
                                                    num_draws = 50,
                                                    num_trees = 10))
             adrf(fit, treat = "logBLL", n = 5)
           }),
      list(case = "a dbarts fit, whose trees must be stored to be sent", pkg = "dbarts",
           make = function() {
             fit <- dbarts::bart2(Math ~ logBLL + Male + PIR + Age,
                                  data = nhanes3lead[1:600, ], keepTrees = TRUE,
                                  verbose = FALSE, n.chains = 1, n.samples = 50,
                                  n.burn = 20)
             adrf(fit, treat = "logBLL", n = 5)
           }),
      list(case = "a stan4bart fit, which is kept in this session", pkg = "stan4bart",
           make = function() {
             fit <- suppressWarnings(stan4bart::stan4bart(
               Math ~ bart(logBLL + Age) + Male, data = nhanes3lead[1:300, ],
               verbose = -1, chains = 1, iter = 40, warmup = 20,
               bart_args = list(keepTrees = TRUE)))
             adrf(fit, treat = "logBLL", n = 5)
           }),
      list(case = "the imputed datasets", pkg = "mice",
           make = function() {
             d <- nhanes3lead[, c("Math", "logBLL", "PIR", "Male")]
             d$PIR[seq(1, nrow(d), by = 10)] <- NA
             imp <- mice::mice(d, m = 3, printFlag = FALSE, seed = 12)
             fits <- with(imp, lm(Math ~ poly(logBLL, 2) * Male + PIR))
             adrf(fits, treat = "logBLL", n = 5)
           })
    )
  )
})

test_that("effect curves do not depend on the plan under another RNG kind", {
  skip_on_cran()
  skip_if_not_installed("betareg")
  skip_if_not_installed("fwb")

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  old_kind <- RNGkind()
  on.exit(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]), add = TRUE)

  make <- function() {
    d <- transform(nhanes3lead, y = (Block + .5) / (max(Block) + 1))
    fit <- betareg::betareg(y ~ splines::ns(logBLL, 3) + Male, data = d)
    a <- adrf(fit, treat = "logBLL", n = 5)

    fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)
    b <- adrf(fit, treat = "logBLL", n = 5, vcov = "boot",
              fwb.args = list(R = 20, verbose = FALSE))

    list(get_vcov(a), attr(b, ".boot")$t, stats::runif(3L))
  }

  set.seed(123, kind = "Knuth-TAOCP-2002")
  one <- make()

  skip_if_not(start_two_workers(), "no second worker available")

  set.seed(123, kind = "Knuth-TAOCP-2002")
  many <- make()

  expect_equal(many[1:2], one[1:2])
  expect_identical(many[[3L]], one[[3L]])
  expect_identical(RNGkind()[1L], "Knuth-TAOCP-2002")
})

test_that("two calls in a row draw different bootstrap replicates", {
  skip_on_cran()
  skip_if_not_installed("fwb")

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  skip_if_not(start_two_workers(), "no second worker available")

  fit <- lm(Math ~ poly(logBLL, 2) * Male + PIR, data = nhanes3lead)

  set.seed(123)
  b1 <- adrf(fit, treat = "logBLL", n = 5, vcov = "boot",
             fwb.args = list(R = 20, verbose = FALSE))
  b2 <- adrf(fit, treat = "logBLL", n = 5, vcov = "boot",
             fwb.args = list(R = 20, verbose = FALSE))

  expect_false(identical(attr(b1, ".boot")$t, attr(b2, ".boot")$t))
})

test_that("units given streams draw the same numbers with or without a plan", {
  skip_on_cran()

  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)

  run <- function() {
    set.seed(123)
    out <- par_lapply(1:5, function(i) stats::runif(2L),
                      seeds = parallel_streams(5L))
    list(out = out, after = stats::runif(3L))
  }

  one <- run()

  skip_if_not(start_two_workers(), "no second worker available")

  many <- run()

  expect_identical(many, one)

  # Kept in this session, the units use the same streams
  set.seed(123)
  local <- par_lapply(1:5, function(i) stats::runif(2L),
                      seeds = parallel_streams(5L), local = TRUE)

  expect_identical(local, one[["out"]])
})

test_that("a refit finds the global objects its call names", {
  e <- new.env(parent = globalenv())
  assign("local_data", data.frame(), envir = e)

  g <- .call_globals(quote(median(local_data, na.rm = T)), e)

  # `median()` is in an attached package, which a new R session may not have
  # attached; `local_data` travels with `e`, and `T` is in base
  expect_named(g, "median")
  expect_identical(g[["median"]], stats::median)
})
