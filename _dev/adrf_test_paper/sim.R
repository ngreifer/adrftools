# Simulation studies

## Sim 1: Comparison of curve_test(Q), curve_test(S), and Wald in setting
## of parametric ADRF
##
## - N: 250, 500, 1000
## - DGP: flat, linear, cubic; confounding present
## - Estimation: spline lm with covariates
## - tests: all tx terms 0 with Wald, flat ADRF (all points equal to midpoint), AMEF = 0 (all slopes 0)

#n: sample size; p: polynomial of
make_data <- function(n = 250, p = 0, k = 10, r = .2, A_r2 = .3, Y_r2 = .3) {
  mu <- rep(0, k)
  Sigma <- r + (1 - r) * diag(k)
  X <- MASS::mvrnorm(n, mu, Sigma)

  #Coefs give appropriate R2 while ensuring conditional var of A is 1
  A_coef_X <- rep(1, k) * sqrt(A_r2 / ((1 - A_r2) * (r * k^2 + (1 - r)*k)))
  A_lp <- drop(X %*% A_coef_X)

  A <- rnorm(n, A_lp, 1)

  if (p == 0) {
    Y_lp_A <- rep(0, n)
  }
  else if (p == 1) {
    Y_lp_A <- 2 * A
  }
  else if (p == 2) {
    Y_lp_A <- - (A + 2.5) * (A - 3) / 10
  }
  else if (p == 3) {
    Y_lp_A <- (A + 3.5) * (A  + .5) * (A - 3) / 10
  }
  else {
    stop("p must be between 0 and 3")
  }

  Y_coef_X <- rep(1, k) * sqrt(Y_r2 / ((1 - Y_r2) * (r * k^2 + (1 - r)*k)))

  Y <- rnorm(n, 4 + X %*% Y_coef_X + Y_lp_A, exp(X %*% Y_coef_X/20))

  cbind(A = A, Y = Y,
        as.data.frame(X))
}

fit_curves <- function(d, np = 41) {
  # f <- sprintf("Y ~ ns(A, df = 4) + (%s)",
  #              paste(c(1, names(d)[-(1:2)]), collapse = " + "))
  f <- sprintf("Y ~ poly(A, 3) * (%s)",
               paste(c(1, names(d)[-(1:2)]), collapse = " + "))
  fit <- lm(as.formula(f), data = d)
  # V <- sandwich::vcovHC(fit, type = "HC0")
  V <- sandwich::sandwich(fit)

  # a <- adrf(fit, treat = "A", vcov = V,
  #           values = seq(-3, 3, length.out = np))
  #
  # m <- amef(fit, treat = "A", vcov = V,
  #           values = seq(-3, 3, length.out = np))

  au <- adrf(fit, treat = "A", vcov = "unconditional", n = np)
  ac <- adrf(fit, treat = "A", vcov = V, n = np)

  mu <- amef(fit, treat = "A", vcov = "unconditional", n = np)
  mc <- amef(fit, treat = "A", vcov = V, n = np)

  list(adrf_u = au,
       adrf_c = ac,
       amef_u = mu,
       amef_c = mc,
       fit = fit,
       data = d,
       V = V)
}

.wald_test <- function(fit, parm, vcov, test = "chisq") {

  m <- 0 * vcov
  m[cbind(parm, parm)] <- 1
  m <- m[rowSums(m != 0) > 0,,drop = FALSE]

  b <- matrix(coef(fit), nrow = 1)

  stat <- (b %*% t(m)) %*% solve(m %*% vcov %*% t(m), m %*% t(b))

  if (test == "chisq")
    pchisq(drop(stat), length(parm), lower.tail = FALSE)
  else
    pf(drop(stat) / length(parm), length(parm), df.residual(fit), lower.tail = FALSE)
}

flat_test <- function(fits, method = "adrf_u", approx = "sim", ...) {
  switch(method,
         "adrf_u" = curve_test(fits$adrf_u, "flat", approx, ...)$p,
         "adrf_c" = curve_test(fits$adrf_c, "flat", approx, ...)$p,
         "amef_u" = curve_test(fits$amef_u, 0, approx, ...)$p,
         "amef_c" = curve_test(fits$amef_c, 0, approx, ...)$p,
         "wald_chi" = .wald_test(fits$fit, 1 + 1:3, vcov = fits$V, test = "chisq"),
         "wald_f" = .wald_test(fits$fit, 1 + 1:3, vcov = fits$V, test = "f"),
         stop(sprintf("%s is not an acceptable method", method)))
}

do_sim <- function(i, grid) {
  dat <- make_data(n = grid$n[i],
                   p = grid$p[i],
                   k = grid$k[i],
                   r = grid$r[i],
                   A_r2 = grid$A_r2[i],
                   Y_r2 = grid$Y_r2[i])

  fits <- fit_curves(dat, np = 51)

  data.frame(i = i
             , adrf_c_pat = flat_test(fits, "adrf_c", "pat")
             , adrf_u_pat_ = flat_test(fits, "adrf_u", "pat")
             # , adrf_c_pt = flat_test(fits, "adrf_c", "pat", df = Inf)
             , adrf_c_imh = flat_test(fits, "adrf_c", "imhof")
             , adrf_u_imh = flat_test(fits, "adrf_u", "imhof")
             , amef_c_pat = flat_test(fits, "amef_c", "pat")
             , amef_u_pat = flat_test(fits, "amef_u", "pat")
             # , amef_c_pat = flat_test(fits, "amef_c", "pat", df = Inf)
             , amef_c_imh = flat_test(fits, "amef_c", "imhof")
             , amef_u_imh = flat_test(fits, "amef_u", "imhof")
             # , wald_chi = flat_test(fits, "wald_chi")
             # , wald_f = flat_test(fits, "wald_f"))
  )
}

grid <- expand.grid(i = 1:2000,
                    n = c(200, 1000),
                    p = c(0, 2),
                    k = 3,
                    r = .2,
                    A_r2 = .3,
                    Y_r2 = .3)

library(splines); library(dplyr)
library(parallel)
cl <- makeCluster(3, type = "FORK")
clusterSetRNGStream(cl, 9876)
# set.seed(123, kind = "L'Ecuyer-CMRG")
res <- do.call("rbind", pbapply::pblapply(seq_len(nrow(grid)), do_sim, grid = grid, cl = cl)); stopCluster(cl)

res_l <- bind_cols(grid, res[-1]) |> tidyr::pivot_longer(cols = names(res)[-1], names_to = "test", values_to = "pval")

bind_cols(grid, res[-1]) |> group_by(n, p) |> summarize(across(names(res)[-1], ~mean(. < .1)))

n <- length(unique(grid$i))
binom_ribbon <- data.frame(x = seq(0, 1, length = 1e3)) |>
  mutate(p = (x * n + 1.96^2/2) / (n + 1.96^2),
         se = sqrt(p * (1 - p) / (n + 1.96^2)))

ggplot(res_l |> filter(p == 0)) + stat_ecdf(aes(x = pval, color = test), geom = "step") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  facet_grid(cols = vars(n)) +
  geom_ribbon(data = binom_ribbon, aes(x = p, ymax = p + 1.96 * se, ymin = p - 1.96 * se), alpha = .2) +
  coord_cartesian(xlim=c(0, .15), ylim = c(0, .2))

## Sim 2: Comparison of curve_test(Q) and curve_test(S) in setting
## of general ADRF
##
## - N: 250, 500, 1000
## - DGP: Cubic CDRF but flat marginal; confounding present
## - Estimation: spline lm with covariates + interactions
## - tests: flat ADRF (all points equal to midpoint), AMEF = 0 (all slopes 0)
## - Notes: need to think about variance; unconditional is more correct
