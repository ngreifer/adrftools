# Stage 1: ADRF estimates from 4000 replicates per scenario cell, whose
# covariance is the oracle used by every test in stage 2 (run-tests.R).
#
# Run from the package root: Rscript _dev/omnibus-test-choice/run-oracle.R [workers]
# Set TEST_CHOICE_SMOKE=1 to run a few units instead.

source("_dev/omnibus-test-choice/common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L
smoke <- nzchar(Sys.getenv("TEST_CHOICE_SMOKE"))

reps <- if (smoke) 30 else 4000

units <- unlist(lapply(seq_len(nrow(cells)), function(i) {
  rep(list(as.list(cells[i, ])), reps)
}), recursive = FALSE)

out_base <- if (smoke) "_dev/omnibus-test-choice/oracle-smoke" else "_dev/omnibus-test-choice/oracle"

set.seed(20261117)

res <- run_units(units, "oracle_unit", "Test choice: oracle covariances", out_base, workers)

# Oracle covariance for each cell
ids <- vapply(res, `[[`, "", "cell")
est <- do.call("rbind", lapply(res, `[[`, "est"))
oracle <- lapply(split(seq_along(ids), ids), function(i) stats::cov(est[i, , drop = FALSE]))

saveRDS(oracle, if (smoke) "_dev/omnibus-test-choice/oracle-cov-smoke.rds" else "_dev/omnibus-test-choice/oracle-cov.rds")
