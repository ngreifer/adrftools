# Summaries of run-equiv.R: rejection rates of omnibus tests with equivalent
# null hypotheses, and how often tests of the same null disagree.
#
# Run from the package root: Rscript _dev/omnibus-size-study/analyze-equiv.R

res <- readRDS("_dev/omnibus-size-study/equiv.rds")
res <- Filter(Negate(is.null), res)
stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

eff <- vapply(res, `[[`, "", "effect")
P <- do.call("rbind", lapply(res, `[[`, "p"))

rates <- do.call("rbind", lapply(unique(eff), function(e) {
  round(colMeans(P[eff == e, , drop = FALSE] < .05), 3)
}))
rownames(rates) <- unique(eff)
cat("Rejection rate at the .05 level (n = 500, 1000 replicates per effect; SE about .007 at .05)\n")
print(t(rates))

# Tests of the same null: ADRF flat vs. reference = 0, flat reference, and AMEF
# = 0; ADRF linear vs. AMEF flat
pairs <- list(c("adrf_flat", "ref0_low"), c("adrf_flat", "ref0_mid"), c("adrf_flat", "ref0_high"),
              c("adrf_flat", "ref_flat_low"), c("adrf_flat", "amef_0"),
              c("adrf_linear", "amef_flat"))

agree <- do.call("rbind", lapply(pairs, function(pr) {
  data.frame(pair = paste(pr, collapse = " vs. "),
             max_abs_diff = signif(max(abs(P[, pr[1L]] - P[, pr[2L]])), 3),
             spearman = round(stats::cor(P[, pr[1L]], P[, pr[2L]], method = "spearman"), 3),
             disagree_05 = round(mean((P[, pr[1L]] < .05) != (P[, pr[2L]] < .05)), 3))
}))
cat("\nAgreement across all replicates\n")
print(agree, row.names = FALSE)

saveRDS(list(rates = rates, agree = agree), "_dev/omnibus-size-study/equiv-summary.rds")
