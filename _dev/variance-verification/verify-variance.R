# Independent checks of the unconditional variance computed by adrf()
#
# Each section compares adrf() with an implementation that shares no code with
# it. See _dev/variance-theory.md for the derivations and references.
#
#   1. stdReg::stdGlm() (Sjolander, 2016): regression standardization with a
#      stacked sandwich variance; with a weighted glm, the prior weights are
#      used both to fit the model and to average the predictions.
#   2. geex::m_estimate() (Saul & Hudgens, 2020): general stacked M-estimation,
#      used where model weights and target weights differ and within subgroups.
#   3. A fractional weighted bootstrap of the whole procedure (fwb).
#   4. survey replicate-weight (jackknife) variances under a clustered and
#      stratified design (issue 13).
#   5. A stacked M-estimator that includes the estimation of continuous
#      treatment weights, and the corresponding check for WeightIt models
#      (issue 12, which turned out not to be a bug).
#
# Requires geex and stdReg (not package dependencies), plus survey and fwb.
# Run from the package root:
#   Rscript _dev/variance-verification/verify-variance.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

stopifnot(requireNamespace("geex", quietly = TRUE),
          requireNamespace("stdReg", quietly = TRUE),
          requireNamespace("survey", quietly = TRUE),
          requireNamespace("fwb", quietly = TRUE))

progress_r <- "~/.claude/skills/live-progress/assets/progress.R"
use_progress <- file.exists(path.expand(progress_r))

if (use_progress) {
  source(path.expand(progress_r))
  pr <- prog_init(total = 13, title = "adrftools variance verification",
                  unit = "check", kind = "verification")
}

tick <- function(label) {
  if (use_progress) {
    prog_tick(pr, label = label)
  }
}

out_dir <- "_dev/variance-verification"

data("nhanes3lead")
d <- nhanes3lead
d$Block12 <- as.numeric(d$Block >= 12)

set.seed(19)
d$wt <- rexp(nrow(d))^2  # deliberately variable target weights
d$wm <- d$MEC_wt / mean(d$MEC_wt)  # model weights

N <- nrow(d)
vals <- c(0, 1, 2)

se <- function(V) sqrt(diag(as.matrix(V)))

results <- list()
record <- function(check, method, values) {
  results[[length(results) + 1L]] <<- data.frame(check = check, method = method,
                                                 se1 = values[1], se2 = values[2],
                                                 se3 = values[3])
}

# The previous (incorrect) weighted influence function, computed by hand for
# the linear model `Math ~ logBLL + PIR + Male`, to show what was fixed
old_formula_se <- function(fit, data, w) {
  IFb <- sandwich::estfun(fit) %*% sandwich::bread(fit)
  tt <- delete.response(terms(fit))
  pieces <- lapply(vals, function(v) {
    data$logBLL <- v
    X <- model.matrix(tt, data = data)
    list(g = colSums(w * X) / sum(w), p = drop(X %*% coef(fit)))
  })
  G <- do.call("rbind", lapply(pieces, `[[`, "g"))
  P <- do.call("cbind", lapply(pieces, `[[`, "p"))
  mu <- colSums(w * P) / sum(w)
  IF <- tcrossprod(IFb, G) + sweep(P, 2L, mu) * sum(w) / sum(w)
  se(crossprod(IF) / N^2)
}

# 1. stdReg -------------------------------------------------------------------

# 1a. Unweighted logistic model
fit <- glm(Block12 ~ logBLL * Male + PIR + Race, data = d, family = binomial)
s <- stdReg::stdGlm(fit, data = d, X = "logBLL", x = vals)
a <- adrf(fit, treat = "logBLL", range = c(0, 2), n = 3)
# stdReg's middle matrix uses var() (denominator n - 1); adrf() and geex use 1/n
record("1a unweighted logistic", "stdReg::stdGlm x (n-1)/n", se(s$vcov * (N - 1) / N))
record("1a unweighted logistic", "adrf", se(attr(a, ".vcov")))
tick("1a stdReg unweighted")

# 1b. Weighted linear model, same weights for fitting and averaging
fit_w <- glm(Math ~ logBLL + PIR + Male, data = d, weights = wt)
s_w <- stdReg::stdGlm(fit_w, data = d, X = "logBLL", x = vals)
a_w <- adrf(fit_w, treat = "logBLL", range = c(0, 2), n = 3, wts = d$wt)
record("1b weighted linear", "stdReg::stdGlm x (n-1)/n", se(s_w$vcov * (N - 1) / N))
record("1b weighted linear", "adrf", se(attr(a_w, ".vcov")))
record("1b weighted linear", "adrf before fix (by hand)", old_formula_se(fit_w, d, d$wt))
tick("1b stdReg weighted")

# 1c. Weighted logistic model with clustering
d$cl <- rep(seq_len(150), length.out = N)
fit_wc <- glm(Block12 ~ logBLL * Male + PIR, data = d, weights = wt,
              family = quasibinomial)
s_wc <- stdReg::stdGlm(fit_wc, data = d, X = "logBLL", x = vals, clusterid = "cl")
a_wc <- adrf(fit_wc, treat = "logBLL", range = c(0, 2), n = 3, wts = d$wt,
             cluster = ~cl)
# With clusters, stdReg's var() over cluster totals times G / n^2 equals the
# G / (G - 1) adjustment adrf() applies, so no rescaling is needed
record("1c weighted logistic, clustered", "stdReg::stdGlm", se(s_wc$vcov))
record("1c weighted logistic, clustered", "adrf", se(attr(a_wc, ".vcov")))
tick("1c stdReg clustered")

# 2. geex ---------------------------------------------------------------------

# Stacked estimating functions for a GLM fit with model weights `wm` and an
# ADRF averaged with target weights `wt` within the units indicated by `in_b`
make_estfun <- function(formula, family, wm_var, wt_var, by_var = NULL, by_level = NULL) {
  tt <- delete.response(terms(formula))
  yname <- all.vars(formula)[1L]

  function(data) {
    x <- model.matrix(tt, data = data)
    y <- data[[yname]]
    wm <- if (is.null(wm_var)) 1 else data[[wm_var]]
    wt <- if (is.null(wt_var)) 1 else data[[wt_var]]
    in_b <- if (is.null(by_var)) 1 else as.numeric(data[[by_var]] == by_level)
    Xa <- do.call("rbind", lapply(vals, function(v) {
      data$logBLL <- v
      model.matrix(tt, data = data)
    }))
    p <- ncol(x)

    function(theta) {
      beta <- theta[seq_len(p)]
      mu <- theta[-seq_len(p)]
      eta <- drop(x %*% beta)
      m <- family$linkinv(eta)
      score <- wm * drop(x) * (y - m) * family$mu.eta(eta) / family$variance(m)
      c(score, wt * in_b * (family$linkinv(drop(Xa %*% beta)) - mu))
    }
  }
}

geex_adrf <- function(formula, family, data, wm_var, wt_var, by_var = NULL,
                      by_level = NULL, units = character(0)) {
  data$.wm <- if (is.null(wm_var)) 1 else data[[wm_var]]
  fit <- suppressWarnings(glm(formula, data = data, family = family, weights = .wm))
  sub <- if (is.null(by_var)) rep(TRUE, nrow(data)) else data[[by_var]] == by_level
  w <- if (is.null(wt_var)) rep(1, nrow(data)) else data[[wt_var]]
  mu <- vapply(vals, function(v) {
    nd <- data[sub, ]
    nd$logBLL <- v
    weighted.mean(predict(fit, newdata = nd, type = "response"), w[sub])
  }, numeric(1L))

  m <- geex::m_estimate(make_estfun(formula, family, wm_var, wt_var, by_var, by_level),
                        data = data, units = units,
                        roots = c(coef(fit), mu), compute_roots = FALSE)

  k <- length(coef(fit)) + seq_along(vals)
  list(est = mu, vcov = m@vcov[k, k])
}

# 2a. Model weights differ from target weights (logistic)
f2 <- Block12 ~ logBLL * Male + PIR
g2 <- geex_adrf(f2, binomial(), d, wm_var = "wm", wt_var = "wt")
fit2 <- glm(f2, data = d, family = quasibinomial, weights = wm)
a2 <- adrf(fit2, treat = "logBLL", range = c(0, 2), n = 3, wts = d$wt)
record("2a model weights != target weights", "geex", se(g2$vcov))
record("2a model weights != target weights", "adrf", se(attr(a2, ".vcov")))
stopifnot(isTRUE(all.equal(unname(g2$est), unname(attr(a2, ".est")), tolerance = 1e-8)))
tick("2a geex weights")

# 2b. The same within a subgroup (by = ~Male, Male == 1)
g2b <- geex_adrf(f2, binomial(), d, wm_var = "wm", wt_var = "wt",
                 by_var = "Male", by_level = 1)
a2b <- adrf(fit2, treat = "logBLL", range = c(0, 2), n = 3, wts = d$wt, by = ~Male)
record("2b subgroup Male = 1", "geex", se(g2b$vcov))
record("2b subgroup Male = 1", "adrf", se(attr(a2b, ".vcov"))[4:6])
tick("2b geex subgroup")

# 3. Bootstrap ------------------------------------------------------------------

# Fractional weighted bootstrap of the whole procedure for 2a
set.seed(99)
fb <- fwb::fwb(d, function(data, w) {
  data$.w <- data$wm * w
  f <- glm(f2, data = data, family = quasibinomial, weights = .w)
  vapply(vals, function(v) {
    nd <- data
    nd$logBLL <- v
    weighted.mean(predict(f, newdata = nd, type = "response"), data$wt * w)
  }, numeric(1L))
}, R = 2000, verbose = FALSE)
record("2a model weights != target weights", "fwb bootstrap (R = 2000)", se(vcov(fb)))
tick("3 bootstrap")

# 4. Survey designs (issue 13) -------------------------------------------------

# A stratified two-stage-style design: 12 strata, 8 PSUs per stratum
set.seed(5)
d$stratum <- as.integer(cut(rank(d$PIR, ties.method = "first"), 12))
d$psu <- interaction(d$stratum, sample(8, N, replace = TRUE), drop = TRUE)
des <- survey::svydesign(ids = ~psu, strata = ~stratum, weights = ~wm, data = d,
                         nest = TRUE)

f4 <- Math ~ logBLL * Male + PIR
sfit <- survey::svyglm(f4, design = des)
a4 <- adrf(sfit, treat = "logBLL", range = c(0, 2), n = 3)
record("4 stratified clustered design", "adrf(svyglm)", se(attr(a4, ".vcov")))

# Design-based jackknife: refit the model and recompute the ADRF in each replicate
rdes <- survey::as.svrepdesign(des, type = "JKn")
jk <- survey::withReplicates(rdes, function(w, data) {
  data$.w <- w
  f <- lm(f4, data = data, weights = .w)
  vapply(vals, function(v) {
    nd <- data
    nd$logBLL <- v
    weighted.mean(predict(f, newdata = nd), w)
  }, numeric(1L))
})
record("4 stratified clustered design", "survey JKn replicates", se(attr(jk, "var")))

# Linearization (Binder, 1983) by hand: the per-unit influence values for the
# weighted mean of predictions, totaled with the design's variance estimator
X <- model.matrix(sfit)
e <- d$Math - drop(X %*% coef(sfit))
W <- sum(d$wm)
A_inv <- solve(crossprod(X * d$wm, X) / W)
tt4 <- delete.response(terms(f4))
z <- vapply(vals, function(v) {
  nd <- d
  nd$logBLL <- v
  Xa <- model.matrix(tt4, data = nd)
  m <- drop(Xa %*% coef(sfit))
  mu <- sum(d$wm * m) / W
  gbar <- colSums(d$wm * Xa) / W
  (m - mu + drop(X %*% A_inv %*% gbar) * e) / W
}, numeric(N))
colnames(z) <- paste0("z", seq_along(vals))
des_z <- update(des, z1 = z[, 1], z2 = z[, 2], z3 = z[, 3])
lin <- survey::svytotal(~z1 + z2 + z3, des_z)
record("4 stratified clustered design", "linearization with svytotal()", se(vcov(lin)))
tick("4 survey design")

# The same three for a design with weights only (no clusters or strata), where
# the corrected adrf() variance applies directly
des0 <- survey::svydesign(ids = ~1, weights = ~wm, data = d)
sfit0 <- survey::svyglm(f4, design = des0)
a40 <- adrf(sfit0, treat = "logBLL", range = c(0, 2), n = 3)
rdes0 <- survey::as.svrepdesign(des0, type = "JK1")
jk0 <- survey::withReplicates(rdes0, function(w, data) {
  data$.w <- w
  f <- lm(f4, data = data, weights = .w)
  vapply(vals, function(v) {
    nd <- data
    nd$logBLL <- v
    weighted.mean(predict(f, newdata = nd), w)
  }, numeric(1L))
})
record("4b weights-only design", "adrf(svyglm)", se(attr(a40, ".vcov")))
record("4b weights-only design", "survey JK1 replicates", se(attr(jk0, "var")))
tick("4b survey weights only")

# 5. Estimated continuous-treatment weights (issue 12) -------------------------

# Stabilized weights f(A) / f(A | Z) with normal densities (Robins, Hernan, &
# Brumback, 2000), estimated by M-estimation together with the outcome model
# and the ADRF. Z = (1, PIR, Male, Age); outcome model Math ~ logBLL * Male + PIR.
Zf <- ~ PIR + Male + Age

check_weights <- function(f5, label) {
  tt5 <- delete.response(terms(f5))

  ps_fit <- lm(update(Zf, logBLL ~ .), data = d)
  alpha <- mean(d$logBLL)
  tau2 <- mean((d$logBLL - alpha)^2)
  sigma2 <- mean(resid(ps_fit)^2)
  d$sw5 <- dnorm(d$logBLL, alpha, sqrt(tau2)) /
    dnorm(d$logBLL, fitted(ps_fit), sqrt(sigma2))
  out_fit <- lm(f5, data = d, weights = sw5)
  mu5 <- vapply(vals, function(v) {
    nd <- d
    nd$logBLL <- v
    mean(predict(out_fit, newdata = nd))
  }, numeric(1L))

  estfun_w <- function(data) {
    Z <- model.matrix(Zf, data = data)
    A <- data$logBLL
    x <- model.matrix(tt5, data = data)
    y <- data$Math
    Xa <- do.call("rbind", lapply(vals, function(v) {
      data$logBLL <- v
      model.matrix(tt5, data = data)
    }))
    q <- ncol(Z)
    p <- ncol(x)

    function(theta) {
      al <- theta[1L]
      ta2 <- theta[2L]
      gam <- theta[2L + seq_len(q)]
      s2 <- theta[3L + q]
      beta <- theta[3L + q + seq_len(p)]
      mu <- theta[-seq_len(3L + q + p)]

      r <- A - drop(Z %*% gam)
      w <- dnorm(A, al, sqrt(ta2)) / dnorm(A, drop(Z %*% gam), sqrt(s2))

      c(A - al,
        (A - al)^2 - ta2,
        drop(Z) * r,
        r^2 - s2,
        w * drop(x) * (y - drop(x %*% beta)),
        drop(Xa %*% beta) - mu)
    }
  }

  theta5 <- c(alpha, tau2, coef(ps_fit), sigma2, coef(out_fit), mu5)
  m5 <- geex::m_estimate(estfun_w, data = d, roots = theta5, compute_roots = FALSE)
  k5 <- length(theta5) - 2:0
  record(label, "geex, weights estimated (full stack)", se(m5@vcov[k5, k5]))

  # The same stack with the weights treated as fixed and known
  estfun_fixed <- make_estfun(f5, gaussian(), wm_var = "sw5", wt_var = NULL)
  m5f <- geex::m_estimate(estfun_fixed, data = d, roots = c(coef(out_fit), mu5),
                          compute_roots = FALSE)
  k5f <- length(coef(out_fit)) + seq_along(vals)
  record(label, "geex, weights treated as fixed", se(m5f@vcov[k5f, k5f]))

  a5 <- adrf(out_fit, treat = "logBLL", range = c(0, 2), n = 3)
  record(label, "adrf(lm(weights = sw))", se(attr(a5, ".vcov")))
  tick(paste(label, "geex"))

  # Bootstrap that re-estimates the weights in every replicate
  set.seed(7)
  fb5 <- fwb::fwb(d, function(data, w) {
    data$.w <- w
    ps <- lm(update(Zf, logBLL ~ .), data = data, weights = .w)
    al <- weighted.mean(data$logBLL, w)
    ta <- sqrt(weighted.mean((data$logBLL - al)^2, w))
    sg <- sqrt(weighted.mean(resid(ps)^2, w))
    data$.sw <- w * dnorm(data$logBLL, al, ta) / dnorm(data$logBLL, fitted(ps), sg)
    f <- lm(f5, data = data, weights = .sw)
    vapply(vals, function(v) {
      nd <- data
      nd$logBLL <- v
      weighted.mean(predict(f, newdata = nd), w)
    }, numeric(1L))
  }, R = 2000, verbose = FALSE)
  record(label, "fwb bootstrap re-estimating weights (R = 2000)",
         se(vcov(fb5)))
  tick(paste(label, "bootstrap"))
}

# The outcome model also adjusts for the confounders, so the weights matter less
f5 <- Math ~ logBLL * Male + PIR
check_weights(f5, "5 estimated weights, adjusted outcome model")

# A marginal structural model in the treatment alone, where the weights carry
# all of the confounding adjustment
check_weights(Math ~ logBLL + I(logBLL^2), "5c estimated weights, MSM outcome model")

# WeightIt: glm_weightit() supplies estfun() and bread() methods whose
# estimating functions include the weight model's contribution (asympt = TRUE),
# so adrf() inherits the M-estimation correction. Compared here with the same
# model refit with the weights treated as fixed, with the delta method on
# vcov(fit), and with a bootstrap that re-estimates the weights.
if (requireNamespace("WeightIt", quietly = TRUE)) {
  f5b <- Math ~ logBLL + I(logBLL^2)
  W <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = d, method = "glm")
  gw <- WeightIt::glm_weightit(f5b, data = d, weightit = W)
  d$w5b <- W$weights
  G5 <- cbind(1, vals, vals^2)

  record("5b WeightIt glm_weightit, MSM", "adrf(glm_weightit)",
         se(attr(adrf(gw, treat = "logBLL", range = c(0, 2), n = 3), ".vcov")))
  record("5b WeightIt glm_weightit, MSM", "delta method on vcov(glm_weightit)",
         se(G5 %*% vcov(gw) %*% t(G5)))
  record("5b WeightIt glm_weightit, MSM", "adrf(lm(weights = w)), weights fixed",
         se(attr(adrf(lm(f5b, data = d, weights = w5b), treat = "logBLL",
                      range = c(0, 2), n = 3), ".vcov")))

  set.seed(3)
  fb5b <- fwb::fwb(d, function(data, w) {
    Wb <- WeightIt::weightit(logBLL ~ PIR + Male + Age, data = data, method = "glm",
                             s.weights = w)
    data$.w <- Wb$weights * w
    drop(G5 %*% coef(lm(f5b, data = data, weights = .w)))
  }, R = 2000, verbose = FALSE)
  record("5b WeightIt glm_weightit, MSM",
         "fwb bootstrap re-estimating weights (R = 2000)", se(vcov(fb5b)))
}
tick("5b WeightIt")

# Output --------------------------------------------------------------------

res <- do.call("rbind", results)
res[c("se1", "se2", "se3")] <- lapply(res[c("se1", "se2", "se3")], signif, 5)
names(res)[3:5] <- sprintf("SE at logBLL = %s", vals)

print(res, row.names = FALSE, right = FALSE)
saveRDS(res, file.path(out_dir, "results.rds"))
write.csv(res, file.path(out_dir, "results.csv"), row.names = FALSE)

if (use_progress) {
  prog_end(pr, "done")
}
