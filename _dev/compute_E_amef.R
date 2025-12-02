.compute_E_amef <- function(treat_var, values, est, vcov, S, by_id, wts, vcov_type) {
  n <- length(values)
  ni <- length(treat_var)

  b_split <- gsplit(seq_len(ni), by_id)

  kw_mat <- matrix(0, nrow = ni, ncol = length(est))
  ik_mat <- matrix(0, nrow = length(b_split), ncol = length(est))

  for (k in seq_along(b_split)) {
    .est_k <- seq_len(n) + (k - 1L) * n
    b <- b_split[[k]]

    kw_k <- get_kernel_w(treat_var[b], values)

    kw_mat[b, .est_k] <- kw_k
    ik_mat[k, .est_k] <- fsum(kw_k) / sum(kw_k)
  }

  E_est <- drop(tcrossprod(est, ik_mat))

  if (vcov_type == "unconditional") {
    S_Eest <- tcrossprod(S, ik_mat)

    for (k in seq_along(b_split)) {
      .est_k <- seq_len(n) + (k - 1L) * n

      S_Eest[, k] <- S_Eest[, k] + (kw_mat[, .est_k, drop = FALSE] %*% est[.est_k] - E_est[k])
    }

    V_est <- crossprod(S_Eest / ni)
  }
  else {
    V_est <- quad_mult(ik_mat, vcov)
  }

  list(E_est = E_est,
       V_est = V_est)
}
