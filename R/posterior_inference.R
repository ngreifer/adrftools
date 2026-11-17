#Posterior inference
posterior_ci <- function(draws, parm, level = .95, simultaneous = FALSE) {
  if (!missing(parm)) {
    draws <- ss(draws, j = parm)
  }

  alpha <- 1 - level

  if (simultaneous && fncol(draws) > 1L) {
    alpha <- 2 * .sup_t_zeta(draws, level)
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

  p_val <- rep.int(NA_real_, fncol(draws))

  # Estimates with no variance (e.g., a reference curve at its reference
  # point) have no p-value, as with Wald inference
  ok <- !.zero_var(draws)

  if (!any(ok)) {
    return(p_val)
  }

  # The largest tail probability at which each pointwise equal-tailed interval
  # contains the null; the pointwise p-value is the smallest alpha at which the
  # interval excludes it, which is twice that
  zeta0 <- vapply(which(ok), function(j) {
    .null_zeta(draws[, j], null)
  }, numeric(1L))

  if (simultaneous && sum(ok) > 1L) {
    # The simultaneous p-value is one minus the largest level at which the
    # simultaneous band excludes the null. The band at level L uses the
    # largest tail probability at which a proportion L of the draws are
    # inside, so it excludes a null that leaves its column's interval above
    # tail probability zeta0 when more than a proportion L of the draws have
    # a larger largest tail probability (see .sup_t_zeta())
    zeta <- .draw_zeta(ss(draws, j = ok))

    p_val[ok] <- vapply(zeta0, function(z0) mean(zeta <= z0), numeric(1L))
  }
  else {
    p_val[ok] <- pmin(1, 2 * pmax(0, zeta0))
  }

  p_val
}

# Columns of draws with no variance, judged relative to the scale of the draws
.zero_var <- function(draws) {
  ranges <- unname(fmax(draws, na.rm = FALSE) - fmin(draws, na.rm = FALSE))

  ranges <= sqrt(.Machine$double.eps) * max(abs(draws), .Machine$double.xmin)
}

# Bayesian and bootstrap sup-t bands of Montiel Olea and Plagborg-Moller
# (2019, Algorithm 2): the band is the product of the equal-tailed percentile
# intervals [Q_j(zeta), Q_j(1 - zeta)] for each column j, with the tail
# probability zeta the largest for which the band contains a proportion
# `level` of the draws. Rather than searching for zeta by bisection, it is
# found exactly. With the type 7 quantiles used by fquantile(), a draw of rank
# r in a column of R draws is inside that column's interval exactly when
# zeta <= min(r - 1, R - r) / (R - 1); a draw is inside the band when this
# holds in every column, so zeta is an order statistic of each draw's smallest
# such bound.
.sup_t_zeta <- function(draws, level) {
  zeta <- .draw_zeta(draws)

  m <- ceiling(level * length(zeta) - sqrt(.Machine$double.eps))

  sort(zeta, decreasing = TRUE)[min(max(m, 1L), length(zeta))]
}

# For each draw, the largest tail probability at which the draw is inside the
# equal-tailed percentile intervals of every column
.draw_zeta <- function(draws) {
  R <- fnrow(draws)

  lo <- apply(draws, 2L, rank, ties.method = "min")
  hi <- apply(draws, 2L, rank, ties.method = "max")

  # The minimum of each row, as column minima of the transpose
  unname(fmin(t(pmin(hi - 1L, R - lo)))) / (R - 1)
}

# The largest tail probability at which `null` is inside the equal-tailed
# percentile interval of `x`, inverting the type 7 quantile function; negative
# when `null` is outside the range of `x`
.null_zeta <- function(x, null) {
  R <- length(x)
  x <- sort(x)

  if (null < x[1L] || null > x[R]) {
    return(-1)
  }

  # Positions h, from 1 to R, at which the quantile function crosses `null`
  h_lo <- stats::approx(x, seq_len(R), xout = null, ties = max)$y
  h_hi <- stats::approx(x, seq_len(R), xout = null, ties = min)$y

  min(h_lo - 1, R - h_hi) / (R - 1)
}
