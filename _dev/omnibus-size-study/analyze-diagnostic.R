# Summaries of run-diagnostic.R: rejection rates by sample size and
# covariance estimator, an oracle test using the covariance of the estimates
# across replicates, and the bias and variability of the estimated variance of
# the test statistic.
#
# Run from the package root: Rscript _dev/omnibus-size-study/analyze-diagnostic.R

source("_dev/omnibus-size-study/common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

res <- readRDS("_dev/omnibus-size-study/diagnostic.rds")
res <- Filter(Negate(is.null), res)

stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

n_rep <- vapply(res, `[[`, numeric(1L), "n")
P <- do.call("rbind", lapply(res, `[[`, "p"))
EST <- do.call("rbind", lapply(res, `[[`, "est"))
TR <- do.call("rbind", lapply(res, `[[`, "tr"))
DF <- do.call("rbind", lapply(res, `[[`, "df"))
D0 <- do.call("rbind", lapply(res, `[[`, "diag_u0"))
D3 <- do.call("rbind", lapply(res, `[[`, "diag_u3"))

sizes <- sort(unique(n_rep))

# Oracle: the flat test using the covariance of the estimates across
# replicates at each sample size in place of an estimate
set.seed(1)
oracle <- numeric(length(res))
for (n in sizes) {
  i <- which(n_rep == n)
  V_mc <- cov(EST[i, ])
  d <- gen_data(n)
  a <- adrf(lm(fmla, data = d), treat = "A", n = n_grid, range = a_range, data = d)
  oracle[i] <- vapply(i, function(j) flat_p(a, EST[j, ], V_mc), numeric(1L))
}
P <- cbind(P, oracle = oracle)

size_tab <- do.call("rbind", lapply(sizes, function(n) {
  i <- n_rep == n
  data.frame(n = n, reps = sum(i), t(round(colMeans(P[i, ] < .05), 3)))
}))

cat("Rejection rate at the .05 level (Monte Carlo SE about",
    round(sqrt(.05 * .95 / min(size_tab$reps)), 4), ")\n")
print(size_tab, row.names = FALSE)

# Bias of the estimated variance of the statistic, E[tr(M V_hat)] / tr(M V_mc),
# and its variability, expressed as Satterthwaite df 2 / CV^2, against the
# average influence-function df
bias_tab <- do.call("rbind", lapply(sizes, function(n) {
  i <- n_rep == n
  tr_mc <- sum(M_flat * cov(EST[i, ]))
  cv2 <- apply(TR[i, ], 2L, function(x) var(x) / mean(x)^2)
  data.frame(n = n,
             t(round(colMeans(TR[i, ]) / tr_mc, 3)),
             df_true_u0 = round(2 / cv2[["u0"]], 1),
             df_true_u3 = round(2 / cv2[["u3"]], 1),
             df_if_u0 = round(mean(DF[i, "u0"]), 1),
             df_if_u3 = round(mean(DF[i, "u3"]), 1))
}))

cat("\nE[tr(M V_hat)] / tr(M V_mc) by estimator, and degrees of freedom of the variance estimate\n")
print(bias_tab, row.names = FALSE)

# Bias of the pointwise variances along the grid
pt_tab <- do.call("rbind", lapply(sizes, function(n) {
  i <- n_rep == n
  v_mc <- apply(EST[i, ], 2L, var)
  data.frame(n = n, point = c("first", "middle", "last"),
             u0 = round((colMeans(D0[i, ]) / v_mc)[c(1, 11, 21)], 3),
             u3 = round((colMeans(D3[i, ]) / v_mc)[c(1, 11, 21)], 3))
}))

cat("\nE[V_hat] / Var(est) at the first, middle, and last grid points\n")
print(pt_tab, row.names = FALSE)

saveRDS(list(size = size_tab, bias = bias_tab, pointwise = pt_tab),
        "_dev/omnibus-size-study/diagnostic-summary.rds")
