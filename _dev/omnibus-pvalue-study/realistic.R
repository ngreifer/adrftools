# Accuracy of the omnibus p-value methods for eigenvalues from actual tests,
# without truncation, compared with a Monte Carlo reference (1e7 draws).
#
# Run from the package root: Rscript _dev/omnibus-pvalue-study/realistic.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

data("nhanes3lead")

# Eigenvalues of the weighted covariance of the tested vector, as used by
# summary.effect_curve()
test_ev <- function(a, hyp) {
  vals <- attr(a, ".values")
  n <- length(vals)
  w <- (c(0, diff(vals)) + c(diff(vals), 0)) / 2
  w <- w / sum(w)
  V <- attr(a, ".vcov")[seq_len(n), seq_len(n)]

  X <- switch(hyp,
              flat = matrix(1, n, 1),
              linear = cbind(1, vals),
              zero = NULL)

  M <- {
    if (is.null(X)) diag(n)
    else diag(n) - X %*% solve(crossprod(X, w * X), t(w * X))
  }

  ev <- eigen(sqrt(w) * t(sqrt(w) * (M %*% V %*% t(M))), symmetric = TRUE,
              only.values = TRUE)$values
  ev[ev > 1e-12 * max(ev)]
}

fit <- lm(Math ~ poly(logBLL, 3) * Male + PIR, data = nhanes3lead)
cases <- list(
  flat_21 = test_ev(adrf(fit, treat = "logBLL", n = 21), "flat"),
  flat_51 = test_ev(adrf(fit, treat = "logBLL", n = 51), "flat"),
  linear_51 = test_ev(adrf(fit, treat = "logBLL", n = 51), "linear"),
  amef_51 = test_ev(amef(adrf(fit, treat = "logBLL", n = 51)), "zero")
)

mc_tail <- function(stat, ev, df, nsim = 1e7, chunk = 1e6) {
  set.seed(1)
  hits <- 0
  for (b in seq_len(nsim / chunk)) {
    Z2 <- matrix(rchisq(chunk * length(ev), 1), nrow = chunk)
    Q <- drop(Z2 %*% ev)
    thr <- if (is.finite(df)) stat * rchisq(chunk, df) / df else stat
    hits <- hits + sum(Q > thr)
  }
  hits / nsim
}

methods <- c("imhof", "davies", "saddlepoint", "satterthwaite", "liu")
rows <- list()

for (nm in names(cases)) {
  ev <- cases[[nm]]
  for (df in c(Inf, 30)) {
    for (target in c(.5, .1, .01, .001)) {
      # The saddlepoint approximation is used only to choose a statistic near
      # each target p-value; the Monte Carlo estimate is the reference
      stat <- uniroot(function(s) {
        log(adrftools:::.quad_form_tail(s, ev, "saddlepoint", df = df)) - log(target)
      }, lower = 1e-3 * sum(ev), upper = 1e3 * sum(ev), tol = 1e-12)$root

      mc <- mc_tail(stat, ev, df)

      for (mt in methods) {
        p <- suppressWarnings(adrftools:::.quad_form_tail(stat, ev, mt, df = df))
        rows[[length(rows) + 1L]] <- data.frame(case = nm, n_ev = length(ev),
                                                ev_ratio = min(ev) / max(ev),
                                                df = df, target = target, mc = mc,
                                                method = mt, p = p)
      }
    }
  }
}

res <- do.call("rbind", rows)
res$mc_se <- sqrt(res$mc * (1 - res$mc) / 1e7)
res$z <- (res$p - res$mc) / res$mc_se
res$rel_err <- (res$p - res$mc) / res$mc

saveRDS(res, "_dev/omnibus-pvalue-study/realistic.rds")

print(unique(res[c("case", "n_ev", "ev_ratio")]), row.names = FALSE)
print(round(with(res, tapply(abs(z), list(method, target), max)), 1))
print(round(with(res, tapply(abs(rel_err), list(method, target), max)), 4))
