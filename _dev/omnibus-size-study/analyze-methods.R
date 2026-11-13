# Summaries of run-methods.R: rejection rates of the candidate fixes, plus
# tests evaluated afterward from the stored covariance estimates:
#   oracle, oracle_dens: the covariance of the estimates across replicates
#   u0_unbiased: HC0 rescaled so that its average equals the oracle's along the
#     statistic (removes the average bias, leaving only its variability)
#   u3_df: HC3 with the influence-function degrees of freedom of HC0
#   avg_u0_u3: the average of HC0 and HC3, which equals HC2 to first order in
#     the leverages (1 / (1 - h) is about the average of 1 and 1 / (1 - h)^2)
#   avg_u0_gj: the average of HC0 and the grouped jackknife (K = 50), which
#     needs only refits of the model, as an approximation to HC2 for any model
#
# Run from the package root: Rscript _dev/omnibus-size-study/analyze-methods.R

source("_dev/omnibus-size-study/common.R")

res <- readRDS("_dev/omnibus-size-study/methods.rds")
res <- Filter(Negate(is.null), res)
stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

n_rep <- vapply(res, `[[`, numeric(1L), "n")
P <- do.call("rbind", lapply(res, `[[`, "p"))
EST <- do.call("rbind", lapply(res, `[[`, "est"))
sizes <- sort(unique(n_rep))

extra <- matrix(NA_real_, length(res), 6L,
                dimnames = list(NULL, c("oracle", "oracle_dens", "u0_unbiased", "u3_df",
                                        "avg_u0_u3", "avg_u0_gj")))

for (n in sizes) {
  i <- which(n_rep == n)
  V_mc <- cov(EST[i, ])

  c_scale <- sum(M_flat * V_mc) / mean(vapply(res[i], function(r) sum(M_flat * r$V0), numeric(1L)))

  extra[i, ] <- t(vapply(res[i], function(r) {
    c(oracle = flat_test(r$est, V_mc, w_trap),
      oracle_dens = flat_test(r$est, V_mc, r$w_dens),
      u0_unbiased = flat_test(r$est, c_scale * r$V0, w_trap),
      u3_df = flat_test(r$est, r$V3, w_trap, df = r$df[["u0"]]),
      avg_u0_u3 = flat_test(r$est, (r$V0 + r$V3) / 2, w_trap),
      avg_u0_gj = flat_test(r$est, (r$V0 + r$Vg50) / 2, w_trap))
  }, numeric(6L)))
}

P <- cbind(P, extra)

rate <- function(alpha) {
  do.call("rbind", lapply(sizes, function(n) {
    i <- n_rep == n
    data.frame(n = n, reps = sum(i), t(round(colMeans(P[i, ] < alpha), 3)))
  }))
}

cat("Rejection rate at the .05 level (Monte Carlo SE about",
    round(sqrt(.05 * .95 / sum(n_rep == sizes[1])), 4), ")\n")
print(rate(.05), row.names = FALSE)

cat("\nRejection rate at the .01 level\n")
print(rate(.01), row.names = FALSE)

saveRDS(list(p05 = rate(.05), p01 = rate(.01), P = P, n = n_rep),
        "_dev/omnibus-size-study/methods-summary.rds")
