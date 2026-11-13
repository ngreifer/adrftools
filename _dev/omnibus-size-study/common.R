# Shared code for the omnibus test size study (audit issue M19)
#
# The data-generating process is the "flat" scenario of
# _dev/omnibus-pvalue-study/calibration.R: the ADRF of A on Y is exactly flat,
# so the omnibus test of flatness has a true null. The outcome model includes
# the X1:X2 term of the data-generating process and has 13 coefficients. (The
# calibration study's model omitted it; because A is correlated with X1 and X2,
# the omitted term made the pseudo-true ADRF of that model depart from flat by
# up to .039, so its null was false and the rejection rate grew with n. Results
# from that model are kept in misspecified/.)

suppressMessages(devtools::load_all(Sys.getenv("ADRF_ROOT", "."), quiet = TRUE))

gen_data <- function(n) {
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, .5)
  A <- .6 * X1 + .4 * X2 + rnorm(n)
  Y <- 1 + .8 * X1 - .5 * X2 + .3 * X1 * X2 + rnorm(n)
  data.frame(Y, A, X1, X2)
}

# Population .025 and .975 quantiles of A (a mixture of two normals with
# means 0 and .4 and SD sqrt(1.36)), so the grid is the same in every sample
# and the covariance of the estimates across replicates is well defined
pA <- function(q) .5 * pnorm(q, 0, sqrt(1.36)) + .5 * pnorm(q, .4, sqrt(1.36))
a_range <- c(uniroot(function(q) pA(q) - .025, c(-5, 5))$root,
             uniroot(function(q) pA(q) - .975, c(-5, 5))$root)

n_grid <- 21L

fmla <- Y ~ poly(A, 3) * (X1 + X2) + X1:X2

# Matrix of the flat test: Q = est' M est, the trapezoid-weighted sum of
# squared deviations from the weighted mean
grid <- seq(a_range[1], a_range[2], length.out = n_grid)
w_trap <- adrftools:::get_trapezoidal_w(grid)
C_flat <- diag(n_grid) - matrix(w_trap, n_grid, n_grid, byrow = TRUE)
M_flat <- crossprod(C_flat * sqrt(w_trap))

# Influence values of the ADRF estimates for an lm fit, one row per unit,
# scaled so that crossprod(S) / N^2 is the unconditional covariance. `adj` is
# applied to the rows of the coefficient part: 1 (HC0), 1 / sqrt(1 - h)
# (HC2), or 1 / (1 - h) (HC3). `cov_term` adds the term for sampling the
# covariates. `Xv` is the list of counterfactual model matrices from
# grid_mm(), computed once per replicate.
grid_mm <- function(fit, d) {
  tt <- delete.response(terms(fit))
  lapply(grid, function(v) {
    d$A <- v
    model.matrix(tt, data = d)
  })
}

uncond_if <- function(fit, d, adj = 1, cov_term = TRUE, Xv = grid_mm(fit, d)) {
  ef <- sandwich::estfun(fit) * adj
  infl <- ef %*% t(sandwich::bread(fit))

  P <- vapply(Xv, function(X) drop(X %*% coef(fit)), numeric(nrow(d)))
  G <- t(vapply(Xv, colMeans, numeric(ncol(Xv[[1L]]))))

  S <- infl %*% t(G)

  if (cov_term) {
    S <- S + sweep(P, 2L, colMeans(P))
  }

  list(S = S, est = colMeans(P), G = G)
}

# Satterthwaite degrees of freedom for the scale of the test statistic,
# tr(M V) with V = crossprod(S) / N^2, from the influence values alone: the
# variance of tr(M V) = sum(a) / N^2 is estimated by sum(a^2) / N^4. Because
# the a are heavy tailed, this overstated the degrees of freedom in the
# diagnostic study (by a factor of about 2.7 at n = 250 and 1.3 at n = 8000).
if_df <- function(S, M = M_flat) {
  a <- rowSums((S %*% M) * S)
  2 * sum(a)^2 / sum(a^2)
}

# P-value of the flat test for given estimates and covariance, using the
# package's own code path
flat_p <- function(a, est, V, df = Inf, method = "imhof") {
  attr(a, ".est") <- est
  attr(a, ".vcov") <- V
  summary(a, hypothesis = "flat", method = method, df = df)$p.value
}

one_rep <- function(n) {
  d <- gen_data(n)
  fit <- lm(fmla, data = d)
  a <- adrf(fit, treat = "A", n = n_grid, range = a_range, data = d)

  N <- n
  p <- length(coef(fit))
  h <- hatvalues(fit)

  u0 <- uncond_if(fit, d)
  u2 <- uncond_if(fit, d, adj = 1 / sqrt(1 - h))
  u3 <- uncond_if(fit, d, adj = 1 / (1 - h))

  V <- list(
    u0 = crossprod(u0$S) / N^2,
    u1 = crossprod(uncond_if(fit, d, adj = sqrt(N / (N - p)))$S) / N^2,
    u2 = crossprod(u2$S) / N^2,
    u3 = crossprod(u3$S) / N^2,
    c0 = u0$G %*% sandwich::vcovHC(fit, type = "HC0") %*% t(u0$G),
    c3 = u0$G %*% sandwich::vcovHC(fit, type = "HC3") %*% t(u0$G)
  )

  # The package's unconditional covariance should equal the hand-built one
  stopifnot(isTRUE(all.equal(unname(attr(a, ".vcov")), unname(V$u0), tolerance = 1e-8)))

  est <- attr(a, ".est")

  dfs <- c(u0 = if_df(u0$S), u2 = if_df(u2$S), u3 = if_df(u3$S))

  pv <- c(
    vapply(V, function(v) flat_p(a, est, v), numeric(1L)),
    u0_df = flat_p(a, est, V$u0, df = dfs[["u0"]]),
    u2_df = flat_p(a, est, V$u2, df = dfs[["u2"]]),
    u3_df = flat_p(a, est, V$u3, df = dfs[["u3"]])
  )

  list(n = n, version = as.character(getNamespaceVersion("adrftools")),
       est = est, p = pv,
       tr = vapply(V, function(v) sum(M_flat * v), numeric(1L)),
       diag_u0 = diag(V$u0), diag_u3 = diag(V$u3),
       df = dfs, max_h = max(h))
}

# Flat test with arbitrary weights `w` over the grid, written out so the
# weighting can be changed; with trapezoidal weights it equals flat_p()
flat_test <- function(est, V, w, df = Inf) {
  k <- length(est)
  C <- diag(k) - matrix(w, k, k, byrow = TRUE)
  stat <- sum(w * drop(C %*% est)^2)
  ev <- eigen(sqrt(w) * t(sqrt(w) * (C %*% V %*% t(C))), symmetric = TRUE,
              only.values = TRUE)$values
  scale <- max(eigen(sqrt(w) * t(sqrt(w) * V), symmetric = TRUE, only.values = TRUE)$values)
  adrftools:::.omnibus_p_value(stat, ev, scale, method = "imhof", df = df)
}

flat_M <- function(w) {
  k <- length(w)
  C <- diag(k) - matrix(w, k, k, byrow = TRUE)
  crossprod(C * sqrt(w))
}

# Candidate fixes: the bootstrap covariance, the bootstrap distribution of the
# statistic, weighting the grid by the estimated density of the treatment, the
# influence-function degrees of freedom, and the grouped jackknife. The covariance estimates are
# returned so other corrections can be evaluated afterward.
one_rep_methods <- function(n, R = 499) {
  d <- gen_data(n)
  fit <- lm(fmla, data = d)
  N <- n
  h <- hatvalues(fit)

  Xv <- grid_mm(fit, d)
  u0 <- uncond_if(fit, d, Xv = Xv)
  u3 <- uncond_if(fit, d, adj = 1 / (1 - h), Xv = Xv)
  V0 <- crossprod(u0$S) / N^2
  V3 <- crossprod(u3$S) / N^2
  est <- u0$est

  ab <- adrf(fit, treat = "A", n = n_grid, range = a_range, data = d,
             vcov = "fwb", fwb.args = list(R = R, verbose = FALSE))
  stopifnot(isTRUE(all.equal(unname(attr(ab, ".est")), unname(est), tolerance = 1e-8)))

  boot_t <- attr(ab, ".boot")[["t"]]
  Vb <- cov(boot_t)

  # Density weights: trapezoidal weights times a kernel density estimate of
  # the treatment at the grid points
  f_hat <- stats::approx(stats::density(d$A), xout = grid)$y
  w_dens <- w_trap * f_hat / sum(w_trap * f_hat)
  M_dens <- flat_M(w_dens)

  # Bootstrap distribution of the centered statistic
  Q <- drop(est %*% M_flat %*% est)
  dev <- sweep(boot_t, 2L, est)
  Q_star <- rowSums((dev %*% M_flat) * dev)

  Q_d <- drop(est %*% M_dens %*% est)
  Q_star_d <- rowSums((dev %*% M_dens) * dev)

  df0 <- if_df(u0$S)
  df0_d <- if_df(u0$S, M_dens)

  Vg50 <- gj_vcov(fit, d, K = 50, Xv = Xv)
  Vg20 <- gj_vcov(fit, d, K = 20, Xv = Xv)

  pv <- c(
    u0 = flat_test(est, V0, w_trap),
    u3 = flat_test(est, V3, w_trap),
    u0_df = flat_test(est, V0, w_trap, df = df0),
    boot_cov = flat_test(est, Vb, w_trap),
    boot_dist = (1 + sum(Q_star >= Q)) / (R + 1),
    u0_dens = flat_test(est, V0, w_dens),
    u3_dens = flat_test(est, V3, w_dens),
    u0_dens_df = flat_test(est, V0, w_dens, df = df0_d),
    boot_dist_dens = (1 + sum(Q_star_d >= Q_d)) / (R + 1),
    gj50 = flat_test(est, Vg50, w_trap),
    gj50_df = flat_test(est, Vg50, w_trap, df = 49),
    gj20_df = flat_test(est, Vg20, w_trap, df = 19)
  )

  list(n = n, version = as.character(getNamespaceVersion("adrftools")),
       est = est, p = pv, V0 = V0, V3 = V3, Vb = Vb, Vg50 = Vg50, w_dens = w_dens,
       df = c(u0 = df0, u0_dens = df0_d))
}

# Grouped (delete-a-group) jackknife covariance of the ADRF estimates: units
# are split at random into K groups, and the model is refit and the ADRF
# recomputed without each group, which deletes the group from both the fit and
# the average over the covariates. Only refitting is needed, so this applies to
# any model. It has K - 1 degrees of freedom.
gj_vcov <- function(fit, d, K = 50, Xv = grid_mm(fit, d)) {
  N <- nrow(d)
  g <- sample(rep_len(seq_len(K), N))

  # The model matrix with a fixed basis spans the same space on any subset of
  # units, so refitting on its rows is the same as refitting the formula
  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))

  th <- t(vapply(seq_len(K), function(k) {
    keep <- g != k
    b <- stats::lm.fit(X[keep, , drop = FALSE], y[keep])$coefficients
    vapply(Xv, function(Xa) mean(Xa[keep, , drop = FALSE] %*% b), numeric(1L))
  }, numeric(length(grid))))

  (K - 1) / K * crossprod(sweep(th, 2L, colMeans(th)))
}
