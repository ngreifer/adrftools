# Size of the omnibus flat test by sample size and covariance estimator
#
# Follow-up: 4000 more replicates at n = 2000, with a different seed, because
# the two runs at that size (2000 replicates in run-diagnostic.R and 1000 in
# run-methods.R) differed by about 3 Monte Carlo SEs.
#
# Each replicate draws data with a flat ADRF, fits the 12-coefficient outcome
# model, and records the flat-test p-value under several covariance estimators
# and degrees-of-freedom corrections (see one_rep() in common.R), along with
# the estimates, so an oracle test using the covariance of the estimates
# across replicates can be computed afterward (analyze-diagnostic.R).
#
# Run from the package root: Rscript _dev/omnibus-size-study/run-diagnostic-2000.R [workers]

source("_dev/omnibus-size-study/common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 2L

units <- rep(2000, 4000)

set.seed(20261113)

pr <- prog_init(total = length(units), title = "Omnibus test size at n = 2000 (follow-up)",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

# Each worker loads the development version of adrftools once (the installed
# version would otherwise be used for adrf() in the workers)
root <- normalizePath(".")

worker_rep <- function(n) {
  # future clears the worker's global environment between chunks, so the
  # study code is kept in an environment stored in an option
  env <- getOption("adrf_size_study_env")

  if (is.null(env)) {
    setwd(root)
    Sys.setenv(ADRF_ROOT = root)
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "_dev/omnibus-size-study/common.R"), envir = env)
    options(adrf_size_study_env = env)
  }

  get("one_rep", envir = env)(n)
}

res <- prog_parallel(pr, units, worker_rep, workers = workers, schedule = "dynamic",
                     chunk_size = 20L,
                     label = function(n) paste0("n = ", n),
                     checkpoint = "_dev/omnibus-size-study/diagnostic-2000")

saveRDS(res, "_dev/omnibus-size-study/diagnostic-2000.rds")

on.exit()
prog_end(pr, "done")
