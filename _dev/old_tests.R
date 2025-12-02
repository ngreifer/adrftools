#Old functions for hypothesis tests
.narrowest_crit_containing_flat <- function(est, vcov, df = Inf, nsim = 1e5) {

  se <- sqrt(diag(vcov))

  zeros <- se < 1e-10

  v <- ss(vcov, !zeros, !zeros) |>
    cov2cor() |>
    .to_psd()

  obj <- function(z) {
    abs(min(est + z * se) - max(est - z * se))
  }

  t_opt <- optimize(obj,
                    lower = 0,
                    upper = 10,
                    tol = 1e-10)$minimum

  1 - mvtnorm::pmvt(lower = alloc(-abs(t_opt), ncol(v)),
                    upper = alloc(abs(t_opt), ncol(v)),
                    df = df, corr = v,
                    keepAttr = FALSE,
                    abseps = 1e-5)
}

.narrowest_crit_containing_zero <- function(est, vcov, df = Inf, nsim = 1e5) {

  se <- sqrt(diag(vcov))

  zeros <- se < 1e-10

  v <- ss(vcov, !zeros, !zeros) |>
    cov2cor() |>
    .to_psd()

  obj_u <- function(z) {
    abs(min(est + z * se))
  }

  t_opt_u <- optimize(obj_u,
                      lower = 0,
                      upper = 10,
                      tol = 1e-10)$minimum

  #Find z such that top of lower CI touches 0
  obj_l <- function(z) {
    abs(max(est - z * se))
  }

  t_opt_l <- optimize(obj_l,
                      lower = 0,
                      upper = 10,
                      tol = 1e-10)$minimum

  t_opt <- max(t_opt_u, t_opt_l)

  1 - mvtnorm::pmvt(lower = alloc(-abs(t_opt), ncol(v)),
                    upper = alloc(abs(t_opt), ncol(v)),
                    df = df, corr = v,
                    keepAttr = FALSE,
                    abseps = 1e-5)
}

.narrowest_crit_containing_flat_boot <- function(boot, parm, ci.type = "perc") {
  #Find CI level such that bottom of upper CI touches top of lower CI
  obj <- function(l) {
    ci <- confint(boot, parm = parm, level = l, ci.type = ci.type,
                  simultaneous = TRUE)
    abs(min(ci[, 2L]) - max(ci[, 1L]))
  }

  ci_opt <- optimize(obj,
                     lower = 1 / nrow(boot[["t"]]),
                     upper = 1 - 1 / nrow(boot[["t"]]),
                     tol = 1e-10)$minimum

  1 - ci_opt
}

.narrowest_crit_containing_zero_boot <- function(boot, parm, hypothesis = 0, ci.type = "perc") {

  #Find CI level such that bottom of upper CI touches hypothesis
  obj_u <- function(l) {
    ci <- confint(boot, parm = parm, level = l, ci.type = ci.type,
                  simultaneous = TRUE)
    abs(min(ci[, 2L]) - hypothesis)
  }

  l_opt_u <- optimize(obj_u,
                      lower = 1 / nrow(boot[["t"]]),
                      upper = 1 - 1 / nrow(boot[["t"]]),
                      tol = 1e-10)$minimum

  #Find CI level such that top of lower CI touches hypothesis
  obj_l <- function(l) {
    ci <- confint(boot, parm = parm, level = l, ci.type = ci.type,
                  simultaneous = TRUE)
    abs(hypothesis - max(ci[, 1L]))
  }

  l_opt_l <- optimize(obj_l,
                      lower = 1 / nrow(boot[["t"]]),
                      upper = 1 - 1 / nrow(boot[["t"]]),
                      tol = 1e-10)$minimum

  l_opt <- max(l_opt_u, l_opt_l)

  1 - l_opt
}

.max_abs_t_dist <- function(v, n = 1e5, max_size = 1e7, df = Inf, abs = TRUE) {
  zeros <- diag(v) < 1e-10
  v <- ss(v, !zeros, !zeros) |>
    cov2cor()

  m <- nrow(v)

  mu <- alloc(0.0, m)

  max_n <- max_size %/% m

  rmvt <- make_rmvt(mu, v, df)

  if (isTRUE(abs)) {
    .abs <- base::abs
  }
  else {
    .abs <- force
  }

  if (n <= max_n) {
    return(fmax(.abs(rmvt(n))))
  }

  iters <- n %/% max_n
  r <- n %% max_n

  z <- unlist(lapply(seq_len(iters), function(i) {
    fmax(.abs(rmvt(max_n)))
  }))

  if (r > 0) {
    z <- c(z, fmax(.abs(rmvt(r))))
  }

  z
}
