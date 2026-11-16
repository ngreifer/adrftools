# Summaries of run-glm-het.R: size of the omnibus flat test with HC0, HC2,
# HC3, the classical covariance, and an oracle (the covariance of the
# estimates across replicates), for linear models with heteroskedastic errors
# and for a logistic model, and the bias of each covariance estimate for the
# test statistic.
#
# Run from the package root: Rscript _dev/omnibus-size-study/analyze-glm-het.R

source("_dev/omnibus-size-study/glm-het-common.R")

res <- readRDS("_dev/omnibus-size-study/glm-het.rds")
res <- Filter(Negate(is.null), res)
stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

pat <- vapply(res, `[[`, "", "pattern")
nn <- vapply(res, `[[`, numeric(1L), "n")
cells <- unique(data.frame(pattern = pat, n = nn))

# One row per cell; the logistic and linear cells record different tests, so
# they are summarized separately
summ <- lapply(seq_len(nrow(cells)), function(k) {
  i <- which(pat == cells$pattern[k] & nn == cells$n[k])
  EST <- do.call("rbind", lapply(res[i], `[[`, "est"))
  V <- cov(EST)

  P <- do.call("rbind", lapply(res[i], `[[`, "p"))
  P <- cbind(P, oracle = vapply(res[i], function(r) flat_test(r$est, V, w_trap), numeric(1L)))

  TR <- do.call("rbind", lapply(res[i], `[[`, "tr"))
  bias <- colMeans(TR) / sum(M_flat * V)

  data.frame(pattern = cells$pattern[k], n = cells$n[k], reps = length(i),
             t(round(colMeans(P < .05), 3)),
             t(setNames(round(colMeans(P < .01), 3), paste0(colnames(P), "_01"))),
             t(setNames(round(bias, 3), paste0("tr_", names(bias)))),
             check.names = FALSE)
})

summ <- lapply(split(summ, cells$pattern), function(x) do.call("rbind", x))

for (pt in names(summ)) {
  cat("\n==", pt, "==\n")
  print(t(summ[[pt]][, -1]))
}

# For the logistic model, whether the logit transformation changed the flat
# test, and whether any fit failed to converge
g <- which(pat == "binomial")
if (length(g)) {
  P <- do.call("rbind", lapply(res[g], `[[`, "p"))
  cat("\nLogistic model: max |p(transformed) - p(untransformed)|:",
      signif(max(abs(P[, "u2_logit"] - P[, "u2"])), 3),
      "  fits not converged:", sum(!vapply(res[g], `[[`, logical(1L), "converged")),
      "  median max leverage:", round(stats::median(vapply(res[g], `[[`, numeric(1L), "max_h")), 3), "\n")
}

saveRDS(summ, "_dev/omnibus-size-study/glm-het-summary.rds")
