#' @exportS3Method stats::coef adrf_est
coef.adrf_est <- function(x, ...) {
  x[["est"]]
}

#' @exportS3Method stats::vcov adrf_est
vcov.adrf_est <- function(x, ...) {
  x[["vcov"]]
}

#' @exportS3Method stats::coef amef_est
coef.amef_est <- coef.adrf_est

#' @exportS3Method stats::vcov amef_est
vcov.amef_est <- vcov.adrf_est
