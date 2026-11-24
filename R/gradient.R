#Compute gradient numerically using centered difference
.gradient <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", .method = "fd",
                      .packages = NULL, .chunk_size = 1L, ...) {
  .method <- arg::match_arg(.method, c("fd", "richardson"))

  if (.method == "fd") {
    .gradientFD(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction,
                .packages = .packages, .chunk_size = .chunk_size, ...)
  }
  else if (.method == "richardson") {
    .gradientRich(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction, ...)
  }
}

#Finite difference gradient. Each column is computed from its own evaluations
#of `.f`, so the columns run on the workers of a future plan when there is one
#(see par_lapply()). `.f` should be a function defined at the top level of the
#package, with what it needs passed through `...`, so that a worker is sent
#those objects and not the frame `.f` was written in; `.packages` and
#`.chunk_size` are passed to par_lapply().
.gradientFD <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center",
                        .packages = NULL, .chunk_size = 1L, ...) {

  .direction <- arg::match_arg(.direction, c("center", "left", "right"))

  if (is_null(.parm)) {
    .parm <- seq_along(.x)
  }

  .eps <- squish(abs(.x) * .eps, lo = .eps, hi = Inf)

  .f0 <- {
    if (.direction != "center") .f(.x, ...)
  }

  cols <- par_lapply(.parm, .fd_column, .f = .f, .x0 = .x, .eps = .eps,
                     .direction = .direction, .f0 = .f0, ...,
                     seeds = parallel_streams(length(.parm)),
                     chunk_size = .chunk_size, packages = .packages)

  matrix(unlist(cols, use.names = FALSE),
         nrow = length(cols[[1L]]), ncol = length(.parm),
         dimnames = list(names(cols[[1L]]), names(.x)[.parm]))
}

#One column of the finite difference gradient, for the `j`th element of `.x0`
.fd_column <- function(j, .f, .x0, .eps, .direction, .f0, ...) {
  .x_r <- .x_l <- .x0

  switch(.direction,
         center = {
           .x_r[j] <- .x0[j] + .eps[j] / 2
           .x_l[j] <- .x0[j] - .eps[j] / 2
           (.f(.x_r, ...) - .f(.x_l, ...)) / .eps[j]
         },
         left = {
           .x_l[j] <- .x0[j] - .eps[j]
           (.f0 - .f(.x_l, ...)) / .eps[j]
         },
         right = {
           .x_r[j] <- .x0[j] + .eps[j]
           (.f(.x_r, ...) - .f0) / .eps[j]
         })
}

#Using Richardson extrapolation
.gradientRich <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", ...) {

  .direction <- arg::match_arg(.direction, c("center", "left", "right"))

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
