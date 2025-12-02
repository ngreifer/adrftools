#Compute gradient numerically using centered difference
.gradient <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", .method = "fd", ...) {
  .method <- match_arg(.method, c("fd", "richardson"))

  if (.method == "fd") {
    .gradientFD(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction, ...)
  }
  else if (.method == "richardson") {
    .gradientRich(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction, ...)
  }
}

#Finite difference gradient
.gradientFD <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", ...) {

  .direction <- match_arg(.direction, c("center", "left", "right"))

  if (is_null(.parm)) {
    .parm <- seq_along(.x)
  }

  .x0 <- .x

  .eps <- squish(abs(.x) * .eps, lo = .eps, hi = Inf)

  if (.direction != "center") {
    .f0 <- .f(.x0, ...)
  }

  for (jj in seq_along(.parm)) {
    j <- .parm[jj]

    if (.direction == "center") {
      .x[j] <- .x0[j] + .eps[j] / 2

      f_new_r <- .f(.x, ...)
    }
    else if (.direction == "left") {
      f_new_r <- .f0
    }
    else if (.direction == "right") {
      .x[j] <- .x0[j] + .eps[j]

      f_new_r <- .f(.x, ...)
    }

    if (j == 1L) {
      jacob <- matrix(0, nrow = length(f_new_r), ncol = length(.parm),
                      dimnames = list(names(f_new_r), names(.x)[.parm]))
    }

    if (.direction == "center") {
      .x[j] <- .x0[j] - .eps[j] / 2

      f_new_l <- .f(.x, ...)
    }
    else if (.direction == "left") {
      .x[j] <- .x0[j] - .eps[j]

      f_new_l <- .f(.x, ...)
    }
    else if (.direction == "right") {
      f_new_l <- .f0
    }

    jacob[, jj] <- (f_new_r - f_new_l) / .eps[j]

    .x[j] <- .x0[j]
  }

  jacob
}

#Using Richardson extrapolation
.gradientRich <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", ...) {

  .direction <- match_arg(.direction, c("center", "left", "right"))

  if (is_null(.parm)) {
    .parm <- seq_along(.x)
  }

  if (.direction != "center") {
    .f0 <- .f(.x, ...)
  }

  n <- length(.x)

  d <- 1e-4
  r <- 4
  v <- 2
  a <- NULL
  h <- abs(d * .x) + .eps * (abs(.x) < 1e-5)

  for (k in seq_len(r)) {
    eps_i <- alloc(0.0, length(.parm))

    for (ii in seq_along(.parm)) {
      i <- .parm[ii]

      eps_i[i] <- h[i]

      a_k_ii <- switch(.direction,
                       "center" = (.f(.x + eps_i, ...) - .f(.x - eps_i, ...)) / (2 * h[i]),
                       "right" = (.f(.x + 2 * eps_i, ...) - .f0) / (2 * h[i]),
                       "left" = (.f0 - .f(.x - 2 * eps_i, ...)) / (2 * h[i]))

      if (is_null(a)) {
        a <- array(NA_real_, dim = c(length(a_k_ii), r, n))
      }

      a[, k, ii] <- a_k_ii

      eps_i[i] <- 0
    }

    h <- h / v
  }

  for (m in seq_len(r - 1L)) {
    a <- (a[, 1L + seq_len(r - m), , drop = FALSE] * (4^m) - a[, seq_len(r - m), , drop = FALSE]) / (4^m - 1)
  }

  array(a, dim = dim(a)[c(1L, 3L)],
        dimnames = list(names(a_k_ii), names(.x)[.parm]))
}
