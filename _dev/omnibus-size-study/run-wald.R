# Omnibus flat test against Wald tests of the same and a stronger hypothesis
#
# Null replicates (2000 per sample size) give the size of each test under
# the current unconditional covariance, the HC2 version, and the classical
# covariance of the coefficients, which is exact for this data-generating
# process; power replicates (1000 per effect, n = 500) compare the two
# statistics with the classical covariance. See one_rep_wald() in
# wald-common.R.
#
# Run from the package root: Rscript _dev/omnibus-size-study/run-wald.R [workers]

source("_dev/omnibus-size-study/wald-common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L

null_units <- lapply(rep(c(250, 500, 1000, 2000), each = 2000),
                     function(n) list(n = n, effect = "none", delta = 0))

alt <- list(linear = .1, quadratic = .1, sigmoid = .3, cubic = .03)
alt_units <- unlist(lapply(names(alt), function(e) {
  rep(list(list(n = 500, effect = e, delta = alt[[e]])), 1000)
}), recursive = FALSE)

units <- c(null_units, alt_units)

# A few units of each kind, to check the script before the full run
smoke <- nzchar(Sys.getenv("WALD_SMOKE"))
if (smoke) {
  units <- units[c(1:3, 8001:8003)]
}

out_base <- if (smoke) "_dev/omnibus-size-study/wald-smoke" else "_dev/omnibus-size-study/wald"

set.seed(20261116)

pr <- prog_init(total = length(units), title = "Omnibus vs. Wald tests",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

root <- normalizePath(".")

worker_rep <- function(u) {
  # future clears the worker's global environment between chunks, so the
  # study code is kept in an environment stored in an option
  env <- getOption("adrf_wald_study_env")

  if (is.null(env)) {
    setwd(root)
    Sys.setenv(ADRF_ROOT = root)
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "_dev/omnibus-size-study/wald-common.R"), envir = env)
    options(adrf_wald_study_env = env)
  }

  get("one_rep_wald", envir = env)(u$n, u$effect, u$delta)
}

res <- prog_parallel(pr, units, worker_rep, workers = workers, schedule = "dynamic",
                     chunk_size = 50L,
                     label = function(u) sprintf("n = %d, %s", u$n, u$effect),
                     checkpoint = out_base)

saveRDS(res, paste0(out_base, ".rds"))

on.exit()
prog_end(pr, "done")
