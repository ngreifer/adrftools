# Stage 2: every test statistic, with the oracle covariance from stage 1
# (run-oracle.R, which must have finished) and with the package's HC2
# covariance, 1000 replicates per scenario cell.
#
# Run from the package root: Rscript _dev/omnibus-test-choice/run-tests.R [workers]
# Set TEST_CHOICE_SMOKE=1 to run a few units instead (with the smoke oracle).

source("_dev/omnibus-test-choice/common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L
smoke <- nzchar(Sys.getenv("TEST_CHOICE_SMOKE"))

stopifnot(file.exists(if (smoke) "_dev/omnibus-test-choice/oracle-cov-smoke.rds" else "_dev/omnibus-test-choice/oracle-cov.rds"))

reps <- if (smoke) 1 else 1000

units <- unlist(lapply(seq_len(nrow(cells)), function(i) {
  rep(list(as.list(cells[i, ])), reps)
}), recursive = FALSE)

out_base <- if (smoke) "_dev/omnibus-test-choice/tests-smoke" else "_dev/omnibus-test-choice/tests"

set.seed(20261118)

run_units(units, "test_unit", "Test choice: tests", out_base, workers)
