# Shared code for the check of the HC2-adjusted unconditional variance with a
# logistic outcome model and with heteroskedastic errors (run-glm-het.R,
# analyze-glm-het.R). Builds on common.R, whose data-generating process and
# outcome model are used throughout; the ADRF is flat in every scenario.

source(file.path(Sys.getenv("ADRF_ROOT", "."), "_dev/omnibus-size-study/common.R"), local = TRUE)

# Error SDs for the linear-model scenarios. "increasing": larger errors at the
# extremes of the treatment, where the leverages are largest; "decreasing":
# smaller errors there; "covariate": heteroskedastic in X1, unrelated to A.
# Each has an average variance near 1, as in the homoskedastic study.
het_sd <- list(
  increasing = function(d) .43 + .52 * abs(d$A - .2),
  decreasing = function(d) 1.6 / (1 + abs(d$A - .2)),
  covariate = function(d) exp(.5 * d$X1 - .25)
)

# Linear model with heteroskedastic errors; the covariance estimators are
# built from the influence values as in the earlier studies (u2 equals what
# adrf() now computes)
one_rep_het <- function(n, pattern) {
  d <- gen_data(n)
  d$Y <- 1 + .8 * d$X1 - .5 * d$X2 + .3 * d$X1 * d$X2 + het_sd[[pattern]](d) * rnorm(n)

  fit <- lm(fmla, data = d)
  h <- hatvalues(fit)

  Xv <- grid_mm(fit, d)
  u0 <- uncond_if(fit, d, Xv = Xv)
  u2 <- uncond_if(fit, d, adj = 1 / sqrt(1 - h), Xv = Xv)
  u3 <- uncond_if(fit, d, adj = 1 / (1 - h), Xv = Xv)
  est <- u0$est

  V <- list(
    u0 = crossprod(u0$S) / n^2,
    u2 = crossprod(u2$S) / n^2,
    u3 = crossprod(u3$S) / n^2,
    c_cl = u0$G %*% stats::vcov(fit) %*% t(u0$G)
  )

  list(model = "lm", pattern = pattern, n = n,
       version = as.character(getNamespaceVersion("adrftools")),
       est = est,
       p = vapply(V, function(v) flat_test(est, v, w_trap), numeric(1L)),
       tr = vapply(V, function(v) sum(M_flat * v), numeric(1L)))
}

# Logistic outcome model with a flat ADRF (the outcome does not depend on A);
# the HC0 variance is obtained by switching off the package's HC2 adjustment.
# Tests are run with and without the default logit transformation.
one_rep_glm <- function(n) {
  d <- gen_data(n)
  d$Y <- stats::rbinom(n, 1, stats::plogis(-.5 + .8 * d$X1 - .5 * d$X2 + .3 * d$X1 * d$X2))

  fit <- suppressWarnings(glm(fmla, data = d, family = binomial))

  a2 <- adrf(fit, treat = "A", n = n_grid, range = a_range, data = d)
  a0 <- testthat::with_mocked_bindings(
    adrf(fit, treat = "A", n = n_grid, range = a_range, data = d),
    .hc2_adjust = function(infl, model) infl,
    .package = "adrftools"
  )

  p_of <- function(a, tr) summary(a, hypothesis = "flat", method = "imhof", transform = tr)$p.value

  list(model = "glm", pattern = "binomial", n = n,
       version = as.character(getNamespaceVersion("adrftools")),
       est = attr(a2, ".est"),
       p = c(u0 = p_of(a0, FALSE), u2 = p_of(a2, FALSE),
             u0_logit = p_of(a0, TRUE), u2_logit = p_of(a2, TRUE)),
       tr = c(u0 = sum(M_flat * attr(a0, ".vcov")), u2 = sum(M_flat * attr(a2, ".vcov"))),
       max_h = max(hatvalues(fit)),
       converged = fit$converged)
}
