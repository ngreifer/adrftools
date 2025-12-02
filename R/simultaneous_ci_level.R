simultaneous_conf_level <- function(vcov, conf_level = .95, df = Inf,
                                    maxiter = 1e5, ptol = 1e-4, ...) {
  conf_level <- process_conf_level(conf_level)

  if (conf_level == 0) {
    return(0)
  }

  chk::chk_not_missing(vcov, "`vcov`")

  if (!is.matrix(vcov) || !is.numeric(vcov) || !all(is.finite(vcov)) ||
      !all(check_if_zero(as.vector(vcov - t(vcov)))) ||
      any(diag(vcov) < 0)) {
    .err("`vcov` must be a valid covariance matrix")
  }

  chk::chk_number(df)
  chk::chk_gt(df, 0)

  if (!all(check_if_zero(diag(vcov) - 1)) ||
      any(vcov < -1) || any(vcov > 1)) {
    zeros <- check_if_zero(diag(vcov))

    vcov <- vcov |>
      ss(!zeros, !zeros) |>
      cov2cor()
  }

  withCallingHandlers({
    t_crit <- mvtnorm::qmvt(conf_level,
                            tail = "both.tails",
                            df = df,
                            corr = vcov,
                            maxiter = maxiter,
                            ptol = ptol,
                            ...)
  },
  error = function(e) {
    .err(sprintf("There was an error computing simultaneous confidence intervals:\n%s",
                 conditionMessage(e)),
         tidy = FALSE)
  })

  1 - 2 * pt(-t_crit$quantile, df)
}
