# Calibration of the omnibus test in summary.effect_curve() under known truth
#
# Data are generated with a known ADRF so that each tested null hypothesis is
# true. Under a valid test, p-values are uniform: about 5% below .05 and 1%
# below .01. Every p-value is also checked for impossible values (outside
# [0, 1] or NaN). The flat-ADRF scenario is repeated with the outcome multiplied
# by 1e-6 to check that results do not depend on its scale.
#
# Run from the package root: Rscript _dev/omnibus-pvalue-study/calibration.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

progress_r <- path.expand("~/.claude/skills/live-progress/assets/progress.R")
use_progress <- file.exists(progress_r)
if (use_progress) {
  source(progress_r)
}

out_file <- "_dev/omnibus-pvalue-study/calibration.rds"

n_obs <- 500
reps <- 400
methods <- c("imhof", "davies", "saddlepoint", "satterthwaite")

gen_data <- function(truth, binary = FALSE) {
  X1 <- rnorm(n_obs)
  X2 <- rbinom(n_obs, 1, .5)
  A <- .6 * X1 + .4 * X2 + rnorm(n_obs)
  f <- switch(truth, flat = 0, linear = .5 * A)
  eta <- 1 + f + .8 * X1 - .5 * X2 + .3 * X1 * X2

  Y <- {
    if (binary) rbinom(n_obs, 1, plogis(eta - 1))
    else eta + rnorm(n_obs)
  }

  data.frame(Y, A, X1, X2)
}

scenarios <- list(
  flat = list(truth = "flat", test = function(a) list(a, "flat")),
  flat_scaled = list(truth = "flat", scale = 1e-6, test = function(a) list(a, "flat")),
  linear = list(truth = "linear", test = function(a) list(a, "linear")),
  binary_flat_transformed = list(truth = "flat", binary = TRUE,
                                 test = function(a) list(a, "flat")),
  binary_flat_untransformed = list(truth = "flat", binary = TRUE, transform = FALSE,
                                   test = function(a) list(a, "flat")),
  reference = list(truth = "flat",
                   test = function(a) list(reference_curve(a, reference = 0), 0)),
  amef = list(truth = "flat", test = function(a) list(amef(a), 0))
)

tasks <- expand.grid(rep = seq_len(reps), scenario = names(scenarios),
                     stringsAsFactors = FALSE)

done <- if (file.exists(out_file)) readRDS(out_file) else NULL
results <- if (is.null(done)) list() else split(done, seq_len(nrow(done)))

if (use_progress) {
  pr <- prog_init(total = nrow(tasks), title = "Omnibus test calibration",
                  unit = "replicate", kind = "simulation")
}

# The seed depends only on the replicate, so scenarios that share a truth use
# the same data (in particular, "flat" and "flat_scaled" differ only in the
# scale of the outcome, so their p-values should be identical)
set.seed(20261112)
seeds <- sample.int(1e8, reps)[tasks$rep]

for (i in seq_len(nrow(tasks))) {
  sc <- scenarios[[tasks$scenario[i]]]
  key <- paste(tasks$scenario[i], tasks$rep[i])

  if (!is.null(done) && key %in% paste(done$scenario, done$rep)) {
    if (use_progress) prog_tick(pr, label = key)
    next
  }

  set.seed(seeds[i])
  d <- gen_data(sc$truth, binary = isTRUE(sc$binary))

  if (!is.null(sc$scale)) {
    d$Y <- d$Y * sc$scale
  }

  fit <- {
    if (isTRUE(sc$binary)) glm(Y ~ poly(A, 3) * (X1 + X2), data = d, family = binomial)
    else lm(Y ~ poly(A, 3) * (X1 + X2), data = d)
  }

  a <- adrf(fit, treat = "A", n = 21, data = d)
  tst <- sc$test(a)

  p <- vapply(methods, function(mt) {
    tryCatch(summary(tst[[1]], hypothesis = tst[[2]], method = mt,
                     transform = if (is.null(sc$transform)) TRUE else sc$transform)$p.value,
             error = function(e) NA_real_)
  }, numeric(1L))

  results[[length(results) + 1L]] <- data.frame(scenario = tasks$scenario[i],
                                                rep = tasks$rep[i],
                                                method = methods, p = unname(p))

  saveRDS(do.call("rbind", results), out_file)

  if (use_progress) prog_tick(pr, label = key)
}

res <- do.call("rbind", results)
saveRDS(res, out_file)

res$impossible <- is.na(res$p) | res$p < 0 | res$p > 1

summ <- do.call("rbind", lapply(split(res, list(res$scenario, res$method), drop = TRUE), function(r) {
  data.frame(scenario = r$scenario[1], method = r$method[1], n = nrow(r),
             impossible = sum(r$impossible),
             reject_05 = mean(r$p < .05, na.rm = TRUE),
             reject_01 = mean(r$p < .01, na.rm = TRUE),
             ks_p = suppressWarnings(ks.test(r$p[!r$impossible], "punif")$p.value))
}))

cat("Monte Carlo SE of a .05 rejection rate with", reps, "replicates:",
    round(sqrt(.05 * .95 / reps), 4), "\n")
print(summ[order(summ$scenario, summ$method), ], row.names = FALSE, digits = 3)

if (use_progress) prog_end(pr, "done")
