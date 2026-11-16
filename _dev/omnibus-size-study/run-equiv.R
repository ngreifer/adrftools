# Size and power of omnibus tests whose null hypotheses are equivalent
# (1000 replicates per effect, n = 500; see equiv-common.R).
#
# Run from the package root: Rscript _dev/omnibus-size-study/run-equiv.R [workers]
# Set EQUIV_SMOKE=1 to run a few units instead.

source("_dev/omnibus-size-study/equiv-common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L

units <- unlist(lapply(names(equiv_effects), function(e) {
  rep(list(list(n = 500, effect = e)), 1000)
}), recursive = FALSE)

smoke <- nzchar(Sys.getenv("EQUIV_SMOKE"))
if (smoke) {
  units <- units[c(1, 1001, 2001)]
}

out_base <- if (smoke) "_dev/omnibus-size-study/equiv-smoke" else "_dev/omnibus-size-study/equiv"

set.seed(20261117)

pr <- prog_init(total = length(units), title = "Equivalent omnibus tests",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

root <- normalizePath(".")

worker_rep <- function(u) {
  # future clears the worker's global environment between chunks, so the
  # study code is kept in an environment stored in an option
  env <- getOption("adrf_equiv_env")

  if (is.null(env)) {
    setwd(root)
    Sys.setenv(ADRF_ROOT = root)
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "_dev/omnibus-size-study/equiv-common.R"), envir = env)
    options(adrf_equiv_env = env)
  }

  get("one_rep_equiv", envir = env)(u$n, u$effect)
}

res <- prog_parallel(pr, units, worker_rep, workers = workers, schedule = "dynamic",
                     chunk_size = 20L,
                     label = function(u) sprintf("%s, n = %d", u$effect, u$n),
                     checkpoint = out_base)

saveRDS(res, paste0(out_base, ".rds"))

on.exit()
prog_end(pr, "done")
