# Summaries of run-tests.R: for each test statistic, its size (in the cells
# where its null hypothesis is true) and power, with the oracle covariance
# (primary) and the HC2 covariance, and how often tests of the same null
# disagree.
#
# Run from the package root: Rscript _dev/omnibus-test-choice/analyze.R

res <- readRDS("_dev/omnibus-test-choice/tests.rds")
res <- Filter(Negate(is.null), res)
stopifnot(all(vapply(res, `[[`, "", "version") == "0.1.0.9001"))

info <- data.frame(dist = vapply(res, `[[`, "", "dist"),
                   shape = vapply(res, `[[`, "", "shape"),
                   n = vapply(res, `[[`, numeric(1L), "n"))

P <- list(oracle = do.call("rbind", lapply(res, `[[`, "p_oracle")),
          hc2 = do.call("rbind", lapply(res, `[[`, "p_hc2")))

flat_tests <- c("adrf_flat_u", "adrf_flat_e",
                paste0(c("ref_min", "ref_q25", "ref_med", "ref_q75", "ref_max"), "_u"),
                paste0(c("ref_min", "ref_q25", "ref_med", "ref_q75", "ref_max"), "_e"),
                "amef_zero_u", "amef_zero_e", "prim_l2", "prim_max", "wald_flat")
linear_tests <- c("adrf_linear_u", "adrf_linear_e", "amef_flat_u", "amef_flat_e", "wald_linear")

# Rejection rate of each test in each cell
rate_table <- function(Pm, tests, alpha = .05) {
  cells <- unique(info)
  out <- do.call("cbind", lapply(seq_len(nrow(cells)), function(k) {
    i <- info$dist == cells$dist[k] & info$shape == cells$shape[k] & info$n == cells$n[k]
    round(colMeans(Pm[i, tests, drop = FALSE] < alpha), 3)
  }))
  colnames(out) <- with(cells, paste(dist, shape, n, sep = "_"))
  out
}

tabs <- list()

for (cv in names(P)) {
  for (fam in c("flat", "linear")) {
    tests <- if (fam == "flat") flat_tests else linear_tests
    tb <- rate_table(P[[cv]], tests)
    tabs[[paste(cv, fam, sep = "_")]] <- tb

    for (dist in c("uniform", "skewed")) for (n in c(300, 1000)) {
      cols <- paste(dist, c("flat", "linear", "quadratic", "cubic"), n, sep = "_")
      cat(sprintf("\n== %s covariance, %s family, %s treatment, n = %d ==\n", cv, fam, dist, n))
      x <- tb[, cols]
      colnames(x) <- c("flat", "linear", "quadratic", "cubic")
      print(x)
    }
  }
}

# Disagreement with the uniform-weighted flat ADRF test (oracle covariance,
# all cells pooled)
dis <- vapply(setdiff(flat_tests, "adrf_flat_u"), function(t) {
  mean((P$oracle[, "adrf_flat_u"] < .05) != (P$oracle[, t] < .05))
}, numeric(1L))
cat("\nShare of replicates in which each flat-family test disagrees with adrf_flat_u at .05 (oracle):\n")
print(round(sort(dis), 3))

dis_lin <- vapply(setdiff(linear_tests, "adrf_linear_u"), function(t) {
  mean((P$oracle[, "adrf_linear_u"] < .05) != (P$oracle[, t] < .05))
}, numeric(1L))
cat("\nShare of replicates in which each linear-family test disagrees with adrf_linear_u at .05 (oracle):\n")
print(round(sort(dis_lin), 3))

saveRDS(list(tables = tabs, disagree_flat = dis, disagree_linear = dis_lin),
        "_dev/omnibus-test-choice/summary.rds")
