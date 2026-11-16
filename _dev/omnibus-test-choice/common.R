# Shared code for the comparison of omnibus test statistics (run-oracle.R,
# run-tests.R, analyze.R).
#
# The aim is a test to recommend for each hypothesis, so only the statistic
# varies: every test in a replicate uses the same covariance matrix of the
# ADRF estimates. The primary covariance is an oracle, the covariance of the
# estimates across independent replicates of the same scenario (run-oracle.R);
# the package's default HC2 unconditional covariance is recorded as a
# secondary, practical lens.
#
# In every scenario the conditional dose-response function differs from the
# ADRF: the heterogeneity term
#   (2 X1 - 1) {(A^2 - 1) + .3 A^3} + .5 X2 A
# averages to 0 at every value of A (X1 is independent of A with
# P(X1 = 1) = .5, and E[X2] = 0), so each unit's curve is cubic while the ADRF
# has the shape added by `shape` (flat, linear, quadratic, or cubic). The
# outcome model, Y ~ poly(A, 3) * (X1 + X2), is correctly specified.

suppressMessages(devtools::load_all(Sys.getenv("ADRF_ROOT", "."), quiet = TRUE))

n_grid <- 21L
sigma <- 3

fmla <- Y ~ poly(A, 3) * (X1 + X2)

# Treatment distributions: uniform and independent of the covariates (as in
# the motivating example), or skewed (4 Beta(2, 4) - 2) and confounded by X2,
# where weighting by the empirical treatment distribution matters
gen_A <- function(n, dist, X2) {
  switch(dist,
         uniform = stats::runif(n, -2, 2),
         skewed = 4 * stats::qbeta(stats::pnorm(.6 * X2 + .8 * stats::rnorm(n)), 2, 4) - 2)
}

q_A <- function(p, dist) {
  switch(dist,
         uniform = stats::qunif(p, -2, 2),
         skewed = 4 * stats::qbeta(p, 2, 4) - 2)
}

# The grid is fixed for each distribution so the oracle covariance is well
# defined: (-1.9, 1.9) for the uniform treatment, as in the motivating
# example, and the population .025 and .975 quantiles for the skewed one
ranges <- list(uniform = c(-1.9, 1.9),
               skewed = q_A(c(.025, .975), "skewed"))

grids <- lapply(ranges, function(r) seq(r[1L], r[2L], length.out = n_grid))

# Shapes of the ADRF and their sizes at n = 1000 (scaled by sqrt(1000 / n)
# so power is comparable across sample sizes), by treatment distribution;
# calibrated in a pilot (80 replicates per value) for roughly 50% power of the
# uniform-weighted flat test with the HC2 covariance. The skewed treatment's
# grid is narrower, so its effects are larger.
shape_fun <- list(
  flat = function(a) 0 * a,
  linear = function(a) a,
  quadratic = function(a) a^2,
  cubic = function(a) a^3
)

effect_sizes <- list(
  uniform = c(flat = 0, linear = .21, quadratic = .21, cubic = .095),
  skewed = c(flat = 0, linear = .45, quadratic = .38, cubic = .19)
)

gen_data <- function(n, dist, shape) {
  X1 <- stats::rbinom(n, 1, .5)
  X2 <- stats::rnorm(n)
  A <- gen_A(n, dist, X2)

  het <- (2 * X1 - 1) * ((A^2 - 1) + .3 * A^3) + .5 * X2 * A
  g <- effect_sizes[[dist]][[shape]] * sqrt(1000 / n) * shape_fun[[shape]](A)

  Y <- g + het + X1 + X2 + sigma * stats::rnorm(n)

  data.frame(Y, A, X1, X2)
}

# Share of the units within the range of the grid that fall in each grid
# point's bin (bins split at the midpoints between grid points), the
# empirical treatment distribution on the grid
empirical_w <- function(A, grid) {
  br <- c(-Inf, (grid[-1L] + grid[-length(grid)]) / 2, Inf)
  inside <- A >= grid[1L] & A <= grid[length(grid)]
  counts <- tabulate(findInterval(A[inside], br), nbins = length(grid))
  counts / sum(counts)
}

one_fit <- function(n, dist, shape, vcov = "unconditional") {
  d <- gen_data(n, dist, shape)
  fit <- lm(fmla, data = d)

  a <- adrf(fit, treat = "A", n = n_grid, range = ranges[[dist]], data = d, vcov = vcov)

  list(est = unname(attr(a, ".est")),
       V = if (vcov != "none") unname(attr(a, ".vcov")),
       wF = empirical_w(d$A, grids[[dist]]),
       a = a)
}

# ---- Linear maps and quadratic forms on the grid --------------------------

# Matrix of an AMEF evaluated at the grid points as a linear map of the ADRF
# estimates at the grid points; it depends only on the grid
amef_map <- function(grid) {
  k <- length(grid)
  d0 <- data.frame(A = seq(grid[1L], grid[k], length.out = 50), Y = stats::rnorm(50))
  a0 <- adrf(lm(Y ~ A, d0), treat = "A", n = k, range = range(grid), vcov = "none")
  stopifnot(isTRUE(all.equal(unname(attr(a0, ".values")), grid)))

  vapply(seq_len(k), function(j) {
    e <- numeric(k)
    e[j] <- 1
    attr(a0, ".est") <- e
    unname(attr(amef(a0), ".est"))
  }, numeric(k))
}

# Deviations from the w-weighted least squares fit of a polynomial of the
# given degree (degree 0: deviations from the weighted mean)
resid_maker <- function(grid, w, degree) {
  X <- outer(grid, 0:degree, `^`)
  diag(length(grid)) - X %*% solve(crossprod(X, X * w), t(X * w))
}

# Westling-type primitive of the centered ADRF under the empirical treatment
# distribution, Omega(a_k) = sum_{j <= k} wF_j (theta_j - gamma), with gamma the
# wF-weighted mean of theta
primitive_map <- function(wF) {
  k <- length(wF)
  lower.tri(diag(k), diag = TRUE) %*% (diag(wF) %*% (diag(k) - matrix(wF, k, k, byrow = TRUE)))
}

# P-value of the quadratic form est' Q est with est ~ N(., V): the
# eigenvalues of B' V B, with Q = B B', weight independent chi-squared(1)
# variables
qf_p <- function(est, V, Q) {
  stat <- drop(crossprod(est, Q %*% est))
  eq <- eigen((Q + t(Q)) / 2, symmetric = TRUE)
  keep <- eq$values > 1e-12 * max(eq$values)
  B <- eq$vectors[, keep, drop = FALSE] %*% diag(sqrt(eq$values[keep]), sum(keep))
  ev <- eigen(crossprod(B, V %*% B), symmetric = TRUE, only.values = TRUE)$values
  adrftools:::.omnibus_p_value(stat, ev, scale = max(ev), method = "imhof")
}

# Wald test of L est = 0 with the generalized inverse of L V L' on its
# leading `rank` eigenvalues
wald_p <- function(est, V, L, rank) {
  dv <- drop(L %*% est)
  e <- eigen(L %*% V %*% t(L), symmetric = TRUE)
  U <- e$vectors[, seq_len(rank), drop = FALSE]
  W <- sum(drop(crossprod(U, dv))^2 / e$values[seq_len(rank)])
  stats::pchisq(W, rank, lower.tail = FALSE)
}

# Maximum test of a vector est ~ N(., V): the share of `nsim` draws whose
# largest absolute value is at least the observed one
max_p <- function(x, V, nsim = 5000) {
  e <- eigen(V, symmetric = TRUE)
  S <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(x))
  Z <- matrix(stats::rnorm(nsim * length(x)), nrow = length(x))
  (1 + sum(apply(abs(S %*% Z), 2L, max) >= max(abs(x)))) / (nsim + 1)
}

# Reference values: the population 0, .25, .5, .75, and 1 quantiles of the
# treatment within the range of the grid, at the nearest grid point
ref_idx <- function(dist) {
  r <- ranges[[dist]]
  pr <- c(0, .25, .5, .75, 1)
  p_A <- function(x) {
    if (dist == "uniform") stats::punif(x, -2, 2) else stats::pbeta((x + 2) / 4, 2, 4)
  }
  lo <- p_A(r[1L])
  hi <- p_A(r[2L])
  vals <- q_A(lo + pr * (hi - lo), dist)
  setNames(vapply(vals, function(v) which.min(abs(grids[[dist]] - v)), integer(1L)),
           c("ref_min", "ref_q25", "ref_med", "ref_q75", "ref_max"))
}

# All tests for one replicate with covariance V. Uniform weighting uses the
# trapezoidal weights; empirical weighting uses the replicate's empirical
# treatment distribution on the grid. Under a cubic outcome model the
# deviations of the ADRF from flat span 3 dimensions and from linear 2.
all_tests <- function(est, V, wF, dist, D) {
  grid <- grids[[dist]]
  k <- length(grid)
  w_u <- adrftools:::get_trapezoidal_w(grid)
  W <- list(u = w_u, e = wF)
  refs <- ref_idx(dist)

  p <- c()

  for (wn in names(W)) {
    w <- W[[wn]]
    C0 <- resid_maker(grid, w, 0)
    C1 <- resid_maker(grid, w, 1)

    p[paste0("adrf_flat_", wn)] <- qf_p(est, V, t(C0) %*% (w * C0))
    p[paste0("adrf_linear_", wn)] <- qf_p(est, V, t(C1) %*% (w * C1))

    for (rn in names(refs)) {
      R <- diag(k)
      R[, refs[[rn]]] <- R[, refs[[rn]]] - 1
      p[paste0(rn, "_", wn)] <- qf_p(est, V, t(R) %*% (w * R))
    }

    p[paste0("amef_zero_", wn)] <- qf_p(est, V, t(D) %*% (w * D))
    p[paste0("amef_flat_", wn)] <- qf_p(est, V, t(C0 %*% D) %*% (w * (C0 %*% D)))
  }

  L <- primitive_map(wF)
  p["prim_l2"] <- qf_p(est, V, t(L) %*% (wF * L))
  p["prim_max"] <- max_p(drop(L %*% est), L %*% V %*% t(L))

  p["wald_flat"] <- wald_p(est, V, resid_maker(grid, w_u, 0), rank = 3L)
  p["wald_linear"] <- wald_p(est, V, resid_maker(grid, w_u, 1), rank = 2L)

  p
}

# Scenario cells
cells <- expand.grid(dist = c("uniform", "skewed"),
                     shape = c("flat", "linear", "quadratic", "cubic"),
                     n = c(300, 1000), stringsAsFactors = FALSE)
cell_id <- function(dist, shape, n) paste(dist, shape, n, sep = "_")

# Runs a study function over units in parallel with the live-progress
# helpers; each worker keeps the study code in an environment stored in an
# option, because future clears the worker's global environment between
# chunks
# (the progress helpers must be sourced at top level by the calling script)
run_units <- function(units, fun_name, title, out_base, workers) {
  pr <- prog_init(total = length(units), title = title, unit = "replicate",
                  kind = "simulation", workers = workers)
  on.exit(prog_end(pr, "failed"), add = TRUE)

  root <- normalizePath(".")

  worker <- function(u) {
    env <- getOption("adrf_test_choice_env")

    if (is.null(env)) {
      setwd(root)
      Sys.setenv(ADRF_ROOT = root)
      env <- new.env(parent = globalenv())
      sys.source(file.path(root, "_dev/omnibus-test-choice/common.R"), envir = env)
      options(adrf_test_choice_env = env)
    }

    get(fun_name, envir = env)(u)
  }

  res <- prog_parallel(pr, units, worker, workers = workers, schedule = "dynamic",
                       chunk_size = 25L,
                       label = function(u) sprintf("%s, %s, n = %d", u$dist, u$shape, u$n),
                       checkpoint = out_base)

  saveRDS(res, paste0(out_base, ".rds"))

  on.exit()
  prog_end(pr, "done")

  invisible(res)
}

# Stage 1: the ADRF estimates only, for the oracle covariance of each cell
oracle_unit <- function(u) {
  f <- one_fit(u$n, u$dist, u$shape, vcov = "none")
  list(cell = cell_id(u$dist, u$shape, u$n), est = f$est)
}

# Stage 2: every test with the oracle covariance and with the HC2 covariance
test_unit <- function(u) {
  if (is.null(.amef_maps)) {
    .amef_maps <<- lapply(grids, amef_map)
    f_or <- if (nzchar(Sys.getenv("TEST_CHOICE_SMOKE"))) "oracle-cov-smoke.rds" else "oracle-cov.rds"
    .oracle <<- readRDS(file.path(Sys.getenv("ADRF_ROOT", "."), "_dev/omnibus-test-choice", f_or))
  }

  id <- cell_id(u$dist, u$shape, u$n)
  f <- one_fit(u$n, u$dist, u$shape)
  D <- .amef_maps[[u$dist]]

  list(cell = id, dist = u$dist, shape = u$shape, n = u$n,
       version = as.character(getNamespaceVersion("adrftools")),
       p_oracle = all_tests(f$est, .oracle[[id]], f$wF, u$dist, D),
       p_hc2 = all_tests(f$est, f$V, f$wF, u$dist, D))
}

.amef_maps <- NULL
.oracle <- NULL
