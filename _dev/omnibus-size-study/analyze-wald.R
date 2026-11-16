# Summaries of run-wald.R: size of the omnibus and Wald tests by covariance
# estimator, the bias of each covariance estimate for the omnibus statistic,
# power under four alternatives, and how often the two tests disagree.
#
# Run from the package root: Rscript _dev/omnibus-size-study/analyze-wald.R

source("_dev/omnibus-size-study/wald-common.R")

res <- readRDS("_dev/omnibus-size-study/wald.rds")
res <- Filter(Negate(is.null), res)
stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

eff <- vapply(res, `[[`, "", "effect")
nn <- vapply(res, `[[`, numeric(1L), "n")
null <- eff == "none"

stack_p <- function(field, i) do.call("rbind", lapply(res[i], `[[`, field))

# Oracle tests: the covariance of the estimates across the null replicates at
# each sample size
sizes <- sort(unique(nn[null]))
oracle <- matrix(NA_real_, length(res), 2L, dimnames = list(NULL, c("omni", "wald")))
V_or <- list()

for (n in sizes) {
  i <- which(null & nn == n)
  EST <- do.call("rbind", lapply(res[i], `[[`, "est"))
  V_or[[as.character(n)]] <- V <- cov(EST)

  oracle[i, ] <- t(vapply(res[i], function(r) {
    c(flat_test(r$est, V, w_trap), wald_flat(r$est, V)[["p"]])
  }, numeric(2L)))
}

rate <- function(P, i, alpha) round(colMeans(P[i, , drop = FALSE] < alpha), 3)

size_table <- function(alpha) {
  do.call("rbind", lapply(sizes, function(n) {
    i <- which(null & nn == n)
    po <- cbind(stack_p("p_omni", i), oracle = oracle[i, "omni"])
    pw <- cbind(stack_p("p_wald", i), oracle = oracle[i, "wald"])
    pc <- stack_p("p_coef", i)
    data.frame(n = n, reps = length(i),
               t(setNames(rate(po, seq_along(i), alpha), paste0("omni_", colnames(po)))),
               t(setNames(rate(pw, seq_along(i), alpha), paste0("wald_", colnames(pw)))),
               t(setNames(rate(pc, seq_along(i), alpha), paste0("coef_", colnames(pc)))),
               check.names = FALSE)
  }))
}

cat("Size at the .05 level (Monte Carlo SE about .005)\n")
s05 <- size_table(.05)
print(t(s05))

cat("\nSize at the .01 level (Monte Carlo SE about .0022)\n")
s01 <- size_table(.01)
print(t(s01))

# Average estimated variance of the omnibus statistic, tr(M V), relative to
# that of the oracle covariance
bias <- do.call("rbind", lapply(sizes, function(n) {
  i <- which(null & nn == n)
  tr <- do.call("rbind", lapply(res[i], `[[`, "tr"))
  round(colMeans(tr) / sum(M_flat * V_or[[as.character(n)]]), 3)
}))
rownames(bias) <- sizes
cat("\nAverage estimated tr(M V) relative to the truth\n")
print(bias)

# Average Wald statistic under the null; its expectation is 3 with the exact
# covariance (3 * df / (df - 2) for the F version)
wstat <- do.call("rbind", lapply(sizes, function(n) {
  i <- which(null & nn == n)
  round(colMeans(do.call("rbind", lapply(res[i], `[[`, "stat_wald"))), 3)
}))
rownames(wstat) <- sizes
cat("\nAverage Wald statistic under the null (3 if the covariance were exact)\n")
print(wstat)

# Power at n = 500
alts <- setdiff(unique(eff), "none")
power <- do.call("rbind", lapply(alts, function(e) {
  i <- which(eff == e)
  po <- stack_p("p_omni", i)
  pw <- stack_p("p_wald", i)
  pc <- stack_p("p_coef", i)
  data.frame(effect = e, delta = res[[i[1L]]]$delta, reps = length(i),
             omni_c_cl_F = mean(po[, "c_cl_F"] < .05), wald_c_cl_F = mean(pw[, "c_cl_F"] < .05),
             omni_u0 = mean(po[, "u0"] < .05), wald_u0 = mean(pw[, "u0"] < .05),
             coef_cl_F = mean(pc[, "cl_F"] < .05),
             both = mean(po[, "c_cl_F"] < .05 & pw[, "c_cl_F"] < .05),
             only_omni = mean(po[, "c_cl_F"] < .05 & pw[, "c_cl_F"] >= .05),
             only_wald = mean(po[, "c_cl_F"] >= .05 & pw[, "c_cl_F"] < .05))
}))
cat("\nPower at the .05 level, n = 500\n")
print(power, digits = 3)

# Agreement under the null, with the classical covariance
i <- which(null)
po <- stack_p("p_omni", i)[, "c_cl_F"]
pw <- stack_p("p_wald", i)[, "c_cl_F"]
agree <- c(cor_p = cor(po, pw), spearman_p = cor(po, pw, method = "spearman"),
           disagree_05 = mean((po < .05) != (pw < .05)),
           both_reject_given_one = sum(po < .05 & pw < .05) / sum(po < .05 | pw < .05))
cat("\nAgreement of the two tests under the null (classical covariance, all n)\n")
print(round(agree, 3))

# Spread of the eigenvalues of the covariance of the deviations from flat; the
# omnibus statistic equals the Wald statistic up to scale when they are equal
ev <- do.call("rbind", lapply(sizes, function(n) {
  V <- V_or[[as.character(n)]]
  e <- eigen(C_flat %*% V %*% t(C_flat), symmetric = TRUE, only.values = TRUE)$values[1:3]
  em <- eigen(sqrt(w_trap) * t(sqrt(w_trap) * (C_flat %*% V %*% t(C_flat))), symmetric = TRUE,
              only.values = TRUE)$values[1:3]
  c(cov = e / e[1L], weighted = em / em[1L])
}))
rownames(ev) <- sizes
cat("\nRelative eigenvalues of the oracle covariance of the deviations (unweighted and trapezoid-weighted)\n")
print(round(ev, 3))

saveRDS(list(size05 = s05, size01 = s01, bias = bias, wald_stat = wstat, power = power,
             agree = agree, eigen = ev),
        "_dev/omnibus-size-study/wald-summary.rds")
