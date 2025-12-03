#Posterior inference
posterior_ci <- function(draws, parm, level = .95, simultaneous = FALSE) {
  if (!missing(parm)) {
    draws <- ss(draws, j = parm)
  }

  alpha <- 1 - level

  if (simultaneous && length(parm) > 1L) {
    alpha <- .pointwise_p_to_simul_p(alpha, draws)
  }

  probs <- c(alpha / 2, 1 - alpha / 2)

  dapply(draws, fquantile, probs = probs,
         return = "matrix", drop = FALSE) |>
    t()
}

posterior_p_value <- function(draws, parm, null = 0, simultaneous = FALSE) {
  if (!missing(parm)) {
    draws <- ss(draws, j = parm)
  }

  p_val <- 2 * pmin(colMeans(draws <= null),
                    colMeans(draws >= null))

  if (simultaneous && length(parm) > 1L) {
    p_val[] <- vapply(p_val, .pointwise_p_to_simul_p, numeric(1L),
                      draws = draws)
  }

  p_val
}

.pointwise_p_to_simul_p <- function(p, draws) {
  level <- 1 - p

  k <- ncol(draws)

  fun <- function(q) {
    interval <- dapply(draws, fquantile, probs = c(q, 1 - q),
                       return = "matrix", drop = FALSE)

    all_est_above <- rowSums(TRA(draws, interval[1L, ]) >= 0) == k
    all_est_below <- rowSums(TRA(draws, interval[2L, ]) <= 0) == k

    coverage <- mean(all_est_above & all_est_below)

    coverage - level
  }

  a <- uniroot(fun, interval = c(p / (2.5 * k), max(p / 2, level / 2)),
               tol = 1e-10)$root

  2 * a
}
