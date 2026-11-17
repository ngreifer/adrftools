data("nhanes3lead", package = "adrftools")

# Average counterfactual predictions computed directly with predict(), used as
# an independent check on adrf() point estimates
manual_adrf <- function(fit, data, treat, values, w = NULL, type = "response") {
  if (is.null(w)) {
    w <- rep(1, nrow(data))
  }

  vapply(values, function(v) {
    data[[treat]] <- v
    stats::weighted.mean(stats::predict(fit, newdata = data, type = type), w)
  }, numeric(1L))
}

# Unconditional (Hansen et al.) variance computed by hand for a linear model:
# influence function of the coefficients plus the sampling term for the
# covariate distribution. With `type = "HC2"` (what adrf() uses), each unit's
# influence on the coefficients is divided by sqrt(1 - h); with "HC1", it is
# multiplied by sqrt(n / (n - p)).
manual_uncond_vcov_lm <- function(fit, data, treat, values, type = "HC2") {
  N <- nrow(data)
  IFb <- sandwich::estfun(fit) %*% sandwich::bread(fit)

  if (type == "HC2") {
    IFb <- IFb / sqrt(1 - stats::hatvalues(fit))
  }
  else if (type == "HC1") {
    IFb <- IFb * sqrt(N / (N - length(stats::coef(fit))))
  }

  tt <- stats::delete.response(stats::terms(fit))

  pieces <- lapply(values, function(v) {
    data[[treat]] <- v
    X <- stats::model.matrix(tt, data = data)
    p <- drop(X %*% stats::coef(fit))
    list(g = colMeans(X), p = p)
  })

  G <- do.call("rbind", lapply(pieces, `[[`, "g"))
  P <- do.call("cbind", lapply(pieces, `[[`, "p"))

  IF <- tcrossprod(IFb, G) + sweep(P, 2L, colMeans(P))

  crossprod(IF) / N^2
}

# The same for a shift curve: each unit's treatment is shifted by each value of
# `shifts` (an additive shift), so the sampling term involves the treatment as
# well as the covariates
manual_uncond_shift_vcov_lm <- function(fit, data, treat, shifts) {
  N <- nrow(data)
  IFb <- sandwich::estfun(fit) %*% sandwich::bread(fit) /
    sqrt(1 - stats::hatvalues(fit))

  tt <- stats::delete.response(stats::terms(fit))

  pieces <- lapply(shifts, function(v) {
    nd <- data
    nd[[treat]] <- data[[treat]] + v
    X <- stats::model.matrix(tt, stats::model.frame(tt, nd))
    list(g = colMeans(X), p = drop(X %*% stats::coef(fit)))
  })

  G <- do.call("rbind", lapply(pieces, `[[`, "g"))
  P <- do.call("cbind", lapply(pieces, `[[`, "p"))

  IF <- tcrossprod(IFb, G) + sweep(P, 2L, colMeans(P))

  crossprod(IF) / N^2
}

get_est <- function(x) unname(attr(x, ".est", exact = TRUE))
get_vcov <- function(x) unname(attr(x, ".vcov", exact = TRUE))
get_values <- function(x) attr(x, ".values", exact = TRUE)
get_se <- function(x) sqrt(diag(get_vcov(x)))

# Layers of a plot() result: the curve, and the confidence band
plot_data <- function(p) ggplot2::layer_data(p, 1L)
ribbon_data <- function(p) {
  geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  ggplot2::layer_data(p, which(geoms == "GeomRibbon"))
}
