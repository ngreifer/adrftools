# Accuracy of the p-value methods for the omnibus test in summary.effect_curve()
#
# The omnibus test statistic is a quadratic form Q = sum_k ev_k * chi^2_1 (or,
# with finite df, compared with stat * chi^2_df / df). This compares each
# method's upper-tail probability, as computed by adrftools:::.quad_form_tail(),
# with an exact reference, across eigenvalue patterns, p-values from .9 down to
# 1e-14, and df = Inf, 2500, and 30.
#
# Exact reference: Ruben's series in upper-tail form,
#   P(Q > q) = sum_k a_k P(chi^2_{n + 2k} > q / beta),
# where beta = min(ev) and a_k are the coefficients of
# prod_j (beta / ev_j)^(1/2) (1 - (1 - beta / ev_j) z)^(-1/2). Every term is
# positive, so there is no cancellation in the tails. With finite df, each term
# is an F tail probability. Equal eigenvalues reduce to exact chi^2 and F
# probabilities, which check the implementation.
#
# Run from the package root: Rscript _dev/omnibus-pvalue-study/accuracy.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

progress_r <- path.expand("~/.claude/skills/live-progress/assets/progress.R")
use_progress <- file.exists(progress_r)
if (use_progress) {
  source(progress_r)
}

ruben_upper <- function(q, ev, df = Inf, tol = 1e-18, max_terms = 2e4) {
  beta <- min(ev)
  gam <- 1 - beta / ev
  n <- length(ev)
  log_a0 <- sum(.5 * log(beta / ev))

  # Number of terms for the coefficient tail to fall below tol
  gmax <- max(gam)
  K <- if (gmax == 0) 0 else min(max_terms, ceiling(log(tol) / log(gmax)) + 50)

  a <- numeric(K + 1)
  a[1] <- exp(log_a0)

  if (K > 0) {
    g <- vapply(seq_len(K), function(r) .5 * sum(gam^r), numeric(1L))

    for (k in seq_len(K)) {
      a[k + 1] <- sum(g[seq_len(k)] * a[k:1]) / k
    }
  }

  dfs <- n + 2 * (0:K)

  tail <- {
    if (is.finite(df)) pf(q / (beta * dfs), dfs, df, lower.tail = FALSE)
    else pchisq(q / beta, dfs, lower.tail = FALSE)
  }

  sum(a * tail)
}

# Eigenvalue patterns. "adrf_*" are taken from actual omnibus tests on
# nhanes3lead, truncated to eigenvalues at least .01 times the largest so the
# exact series converges.
data("nhanes3lead")
fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
a <- adrf(fit, treat = "logBLL", n = 21)
vals <- attr(a, ".values")
w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
w <- w / sum(w)
V <- attr(a, ".vcov")
ev_flat <- eigen(sqrt(w) * t(sqrt(w) * ((diag(21) - matrix(w, 21, 21, byrow = TRUE)) %*%
                                         V %*% t(diag(21) - matrix(w, 21, 21, byrow = TRUE)))),
                 symmetric = TRUE, only.values = TRUE)$values
ev_flat <- ev_flat[ev_flat > .01 * max(ev_flat)]

patterns <- list(
  equal_1 = 1,
  equal_5 = rep(1, 5),
  equal_20 = rep(1, 20),
  equal_50 = rep(1, 50),
  two = c(1, .5),
  two_far = c(1, .02),
  dominant = c(1, rep(.02, 9)),
  geometric = .75^(0:15),
  geometric_fast = .5^(0:6),
  adrf_flat = ev_flat / max(ev_flat)
)

target_p <- c(.9, .5, .1, .01, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14)
dfs <- c(Inf, 2500, 30)

methods <- list(
  imhof = list(method = "imhof"),
  imhof_tight = list(method = "imhof", epsabs = 1e-15, epsrel = 1e-10, limit = 1e5),
  davies = list(method = "davies"),
  davies_tight = list(method = "davies", acc = 1e-14, lim = 1e6),
  liu = list(method = "liu"),
  satterthwaite = list(method = "satterthwaite"),
  saddlepoint = list(method = "saddlepoint")
)

grid <- expand.grid(pattern = names(patterns), df = dfs, target = target_p,
                    stringsAsFactors = FALSE)

if (use_progress) {
  pr <- prog_init(total = nrow(grid), title = "Omnibus p-value accuracy",
                  unit = "scenario", kind = "study")
}

rows <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  ev <- patterns[[grid$pattern[i]]]
  df <- grid$df[i]
  target <- grid$target[i]

  # The statistic whose exact upper-tail probability is the target
  stat <- uniroot(function(s) log(ruben_upper(s, ev, df)) - log(target),
                  lower = 1e-8, upper = 1e4 * sum(ev), tol = 1e-12)$root
  exact <- ruben_upper(stat, ev, df)

  if (length(unique(ev)) == 1L) {
    k <- length(ev)
    closed <- {
      if (is.finite(df)) pf(stat / (ev[1] * k), k, df, lower.tail = FALSE)
      else pchisq(stat / ev[1], k, lower.tail = FALSE)
    }
    stopifnot(isTRUE(all.equal(exact, closed, tolerance = 1e-10)))
  }

  rows[[i]] <- do.call("rbind", lapply(names(methods), function(nm) {
    args <- methods[[nm]]
    p <- tryCatch(suppressWarnings(do.call(adrftools:::.quad_form_tail,
                                           c(list(stat = stat, ev = ev, df = df), args))),
                  error = function(e) NA_real_)
    data.frame(pattern = grid$pattern[i], df = df, target = target, exact = exact,
               method = nm, p = p)
  }))

  if (use_progress) {
    prog_tick(pr, label = sprintf("%s, df = %s, p = %g", grid$pattern[i], df, target))
  }
}

res <- do.call("rbind", rows)

# Hybrid candidates: an integration method when its result is valid and at
# least `cut`, otherwise the saddlepoint approximation
add_hybrid <- function(res, base, cut) {
  b <- res[res$method == base, ]
  sp <- res[res$method == "saddlepoint", ]
  stopifnot(identical(b$pattern, sp$pattern), identical(b$target, sp$target),
            identical(b$df, sp$df))
  use_b <- !is.na(b$p) & b$p >= cut & b$p <= 1
  h <- b
  h$method <- sprintf("%s_saddle_%g", base, cut)
  h$p <- ifelse(use_b, b$p, sp$p)
  rbind(res, h)
}

for (cut in c(1e-2, 1e-3, 1e-4)) {
  res <- add_hybrid(res, "imhof", cut)
  res <- add_hybrid(res, "davies", cut)
}

res$impossible <- is.na(res$p) | res$p < 0 | res$p > 1
res$log10_ratio <- ifelse(res$impossible | res$p <= 0, NA, log10(res$p / res$exact))

out_dir <- "_dev/omnibus-pvalue-study"
saveRDS(res, file.path(out_dir, "accuracy.rds"))

# Summaries: share of impossible values, and the largest error in log10 units
# (1 = off by a factor of 10), overall and by the size of the exact p-value
res$p_range <- cut(res$exact, c(0, 1e-9, 1e-5, 1e-3, 1),
                   labels = c("< 1e-9", "1e-9 to 1e-5", "1e-5 to .001", "> .001"))

summ <- aggregate(cbind(impossible = res$impossible,
                        max_abs_log10 = abs(res$log10_ratio)) ~ method + p_range,
                  data = transform(res, log10_ratio = ifelse(is.na(log10_ratio), 99, log10_ratio)),
                  FUN = function(x) round(max(x), 3))
summ_imp <- aggregate(impossible ~ method, data = res, FUN = mean)

print(summ_imp)
print(summ[order(summ$method, summ$p_range), ], row.names = FALSE)

if (use_progress) {
  prog_end(pr, "done")
}
