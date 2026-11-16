# Size of the omnibus flat test with the HC2-adjusted unconditional variance
# for a logistic outcome model and for linear models with heteroskedastic
# errors (2000 replicates per cell; see glm-het-common.R).
#
# Run from the package root: Rscript _dev/omnibus-size-study/run-glm-het.R [workers]
# Set GLM_HET_SMOKE=1 to run a few units of each kind instead.

source("_dev/omnibus-size-study/glm-het-common.R")
source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

args <- commandArgs(trailingOnly = TRUE)
workers <- if (length(args)) as.integer(args[1]) else 1L

reps <- 2000

het_units <- unlist(lapply(names(het_sd), function(pt) {
  lapply(rep(c(250, 500, 1000), each = reps), function(n) list(kind = "het", n = n, pattern = pt))
}), recursive = FALSE)

glm_units <- lapply(rep(c(500, 1000, 2000), each = reps),
                    function(n) list(kind = "glm", n = n, pattern = "binomial"))

units <- c(het_units, glm_units)

smoke <- nzchar(Sys.getenv("GLM_HET_SMOKE"))
if (smoke) {
  units <- units[c(1:2, length(het_units) + 1:2)]
}

out_base <- if (smoke) "_dev/omnibus-size-study/glm-het-smoke" else "_dev/omnibus-size-study/glm-het"

set.seed(20261116)

pr <- prog_init(total = length(units), title = "HC2 size: GLM and heteroskedasticity",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

root <- normalizePath(".")

worker_rep <- function(u) {
  # future clears the worker's global environment between chunks, so the
  # study code is kept in an environment stored in an option
  env <- getOption("adrf_glm_het_env")

  if (is.null(env)) {
    setwd(root)
    Sys.setenv(ADRF_ROOT = root)
    env <- new.env(parent = globalenv())
    sys.source(file.path(root, "_dev/omnibus-size-study/glm-het-common.R"), envir = env)
    options(adrf_glm_het_env = env)
  }

  if (u$kind == "glm") get("one_rep_glm", envir = env)(u$n)
  else get("one_rep_het", envir = env)(u$n, u$pattern)
}

res <- prog_parallel(pr, units, worker_rep, workers = workers, schedule = "dynamic",
                     chunk_size = 25L,
                     label = function(u) sprintf("%s, n = %d", u$pattern, u$n),
                     checkpoint = out_base)

saveRDS(res, paste0(out_base, ".rds"))

on.exit()
prog_end(pr, "done")
