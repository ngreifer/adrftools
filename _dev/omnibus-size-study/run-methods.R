# Candidate fixes for the size of the omnibus flat test
#
# Each replicate records the flat-test p-value under the candidate fixes in
# one_rep_methods() (common.R): the fractional weighted bootstrap covariance,
# the bootstrap distribution of the statistic, density weights over the grid,
# and influence-function degrees of freedom. The covariance estimates are kept
# so other corrections can be evaluated afterward (analyze-methods.R).
#
# Run from the package root: Rscript _dev/omnibus-size-study/run-methods.R [workers]

source("_dev/omnibus-size-study/common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 2L

# Fewer replicates at n = 2000, where each takes longest
sizes <- c(250, 500, 2000)
reps <- c(1500, 1500, 1000)

units <- rep(sizes, times = reps)

set.seed(20261112)

pr <- prog_init(total = length(units), title = "Omnibus test size: candidate fixes",
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

  get("one_rep_methods", envir = env)(n)
}

res <- prog_parallel(pr, units, worker_rep, workers = workers, schedule = "dynamic",
                     chunk_size = 10L,
                     label = function(n) paste0("n = ", n),
                     checkpoint = "_dev/omnibus-size-study/methods")

saveRDS(res, "_dev/omnibus-size-study/methods.rds")

on.exit()
prog_end(pr, "done")
