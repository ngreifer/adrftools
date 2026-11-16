# Shared code for the comparison of the omnibus flat test with Wald tests
# (run-wald.R, analyze-wald.R). Builds on common.R.
#
# For the linear outcome model of common.R, the ADRF is a cubic polynomial in
# the treatment, theta(a) = G(a) beta, so its deviations from its weighted mean
# on the grid, C theta, lie in a three-dimensional space. Flatness of the ADRF
# is then the rank-3 linear hypothesis C G beta = 0, which a Wald test can test
# directly. The omnibus test measures the same deviations with a fixed metric
# (the trapezoid-weighted L2 norm over the treatment range) rather than the
# inverse of their covariance.

source(file.path(Sys.getenv("ADRF_ROOT", "."), "_dev/omnibus-size-study/common.R"), local = TRUE)

# Effects of the treatment added to the outcome for the power comparison; the
# ADRF is then delta * f(a) plus a constant
effects <- list(
  none = function(a) 0 * a,
  linear = function(a) a,
  quadratic = function(a) a^2,
  sigmoid = function(a) stats::plogis(4 * (a - 1)),
  cubic = function(a) a^3
)

# Wald test of C theta = 0 using the generalized inverse of C V C' on its
# leading eigenvalues; `rank` is fixed at 3 (the dimension of the deviations of
# a cubic from its mean) unless given. Returns the statistic, the rank, and the
# chi-squared and F p-values.
wald_flat <- function(est, V, df = Inf, rank = 3L) {
  d <- drop(C_flat %*% est)
  Vd <- C_flat %*% V %*% t(C_flat)
  e <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE)
  U <- e$vectors[, seq_len(rank), drop = FALSE]
  W <- sum(drop(crossprod(U, d))^2 / e$values[seq_len(rank)])

  c(stat = W,
    p = if (is.finite(df)) stats::pf(W / rank, rank, df, lower.tail = FALSE)
        else stats::pchisq(W, rank, lower.tail = FALSE),
    cond = e$values[1L] / e$values[rank])
}

# Wald test that all coefficients involving the treatment are 0, which is the
# stronger null of no effect for any covariate pattern
wald_coef <- function(b, V, idx, df = Inf) {
  bb <- b[idx]
  W <- drop(crossprod(bb, solve(V[idx, idx], bb)))
  q <- length(idx)

  if (is.finite(df)) stats::pf(W / q, q, df, lower.tail = FALSE)
  else stats::pchisq(W, q, lower.tail = FALSE)
}

one_rep_wald <- function(n, effect = "none", delta = 0) {
  d <- gen_data(n)
  d$Y <- d$Y + delta * effects[[effect]](d$A)

  fit <- lm(fmla, data = d)
  N <- n
  p <- length(coef(fit))
  df_res <- N - p
  h <- hatvalues(fit)

  Xv <- grid_mm(fit, d)
  u0 <- uncond_if(fit, d, Xv = Xv)
  u2 <- uncond_if(fit, d, adj = 1 / sqrt(1 - h), Xv = Xv)
  G <- u0$G
  est <- u0$est
  b <- coef(fit)

  Vb <- list(
    cl = stats::vcov(fit),
    hc0 = sandwich::vcovHC(fit, type = "HC0"),
    hc2 = sandwich::vcovHC(fit, type = "HC2")
  )

  V <- list(
    u0 = crossprod(u0$S) / N^2,
    u2 = crossprod(u2$S) / N^2,
    c_cl = G %*% Vb$cl %*% t(G),
    c0 = G %*% Vb$hc0 %*% t(G),
    c2 = G %*% Vb$hc2 %*% t(G)
  )

  # Omnibus test; with the classical covariance, also with the residual df
  p_omni <- c(vapply(V, function(v) flat_test(est, v, w_trap), numeric(1L)),
              c_cl_F = flat_test(est, V$c_cl, w_trap, df = df_res))

  # Wald test of the same hypothesis
  wf <- lapply(V, wald_flat, est = est)
  p_wald <- c(vapply(wf, `[[`, numeric(1L), "p"),
              c_cl_F = wald_flat(est, V$c_cl, df = df_res)[["p"]])

  # Wald test that every coefficient involving the treatment is 0
  idx <- grep("poly(A", names(b), fixed = TRUE)
  p_coef <- c(cl_F = wald_coef(b, Vb$cl, idx, df = df_res),
              hc0 = wald_coef(b, Vb$hc0, idx),
              hc2 = wald_coef(b, Vb$hc2, idx))

  list(n = n, effect = effect, delta = delta,
       version = as.character(getNamespaceVersion("adrftools")),
       est = est,
       p_omni = p_omni, p_wald = p_wald, p_coef = p_coef,
       stat_omni = drop(est %*% M_flat %*% est),
       stat_wald = vapply(wf, `[[`, numeric(1L), "stat"),
       cond_wald = vapply(wf, `[[`, numeric(1L), "cond"),
       tr = vapply(V, function(v) sum(M_flat * v), numeric(1L)))
}
