# Independent checks of the unconditional variance of shift_curve() and ade()
#
# A shift curve averages m_i(g(A_i, delta); beta), the prediction for unit i at
# its shifted treatment, which is still a function of unit i's own data, so the
# stacked estimating equations of _dev/variance-theory.md apply with the shifted
# treatment in place of the fixed treatment value. These checks confirm that
# shift_curve() computes them, using implementations that share no code with
# it:
#
#   A. A hand-computed HC2 influence function for a linear model.
#   B. geex::m_estimate() stacked M-estimation (HC0) for a weighted logistic
#      model with target weights that differ from the model weights, for an
#      additive shift and for a user-supplied proportional shift within a
#      subgroup.
#   C. The error of the ADE from the interpolated shift curve, relative to an
#      exact central difference of the predictions, and its SE relative to the
#      exact influence function.
#   D. A fractional weighted bootstrap written here (not fwb), for the shift
#      curve and the exact ADE.
#
# Run from the package root:
#   Rscript _dev/variance-verification/verify-shift.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

stopifnot(requireNamespace("geex", quietly = TRUE))

source(path.expand("~/.claude/skills/live-progress/assets/progress.R"))

out_dir <- "_dev/variance-verification"

data("nhanes3lead")
d <- nhanes3lead
d$Block12 <- as.numeric(d$Block >= 12)

set.seed(19)
d$wt <- rexp(nrow(d))^2  # deliberately variable target weights
d$wm <- d$MEC_wt / mean(d$MEC_wt)  # model weights

N <- nrow(d)
se <- function(V) sqrt(diag(as.matrix(V)))

results <- list()
record <- function(check, method, values) {
  results[[length(results) + 1L]] <<- data.frame(check = check, method = method,
                                                 se1 = values[1], se2 = values[2],
                                                 se3 = values[3])
}

# shift_curve() with the HC0 (unadjusted) unconditional variance, which is what
# geex computes
hc0 <- function(code) {
  orig <- adrftools:::.coef_influence
  testthat::with_mocked_bindings(
    code,
    .coef_influence = function(model, hc2 = TRUE, cluster = NULL) {
      orig(model, hc2 = FALSE, cluster = cluster)
    },
    .package = "adrftools"
  )
}

shifts <- c(-.4, 0, .4)

n_boot <- 10000L
pr <- prog_init(total = 4L + n_boot, title = "shift_curve() variance verification",
                unit = "check", kind = "verification")
on.exit(prog_end(pr, "failed"), add = TRUE)

# A. Hand-computed HC2 influence function --------------------------------------

fA <- Math ~ poly(logBLL, 3) * Male + PIR + Age + Race
fitA <- lm(fA, data = d)
ttA <- delete.response(terms(fitA))

X_at <- function(fit, tt, data, a) {
  data$logBLL <- a
  model.matrix(tt, model.frame(tt, data, xlev = fit$xlevels))
}

X <- model.matrix(fitA)
h <- hatvalues(fitA)
IF_beta <- (residuals(fitA) / sqrt(1 - h)) * X %*% solve(crossprod(X) / N)

hand_shift <- function(a_fun) {
  pieces <- lapply(shifts, function(v) {
    Xa <- X_at(fitA, ttA, d, a_fun(d$logBLL, v))
    m <- drop(Xa %*% coef(fitA))
    list(IF = (m - mean(m)) + drop(IF_beta %*% colMeans(Xa)), est = mean(m))
  })
  IF <- do.call("cbind", lapply(pieces, `[[`, "IF"))
  list(est = vapply(pieces, `[[`, numeric(1L), "est"), vcov = crossprod(IF) / N^2)
}

hA <- hand_shift(function(a, v) a + v)
sA <- shift_curve(fitA, treat = "logBLL", shift = c(-.4, .4), n = 3)
stopifnot(isTRUE(all.equal(hA$est, unname(attr(sA, ".est")), tolerance = 1e-10)))
record("A lm, additive, HC2", "by hand", se(hA$vcov))
record("A lm, additive, HC2", "shift_curve", se(attr(sA, ".vcov")))
prog_tick(pr, label = "A hand HC2")

# B. geex stacked M-estimation (HC0) -------------------------------------------

make_estfun <- function(formula, family, a_fun, by_var = NULL, by_level = NULL) {
  tt <- delete.response(terms(formula))
  yname <- all.vars(formula)[1L]

  function(data) {
    x <- model.matrix(tt, data = data)
    y <- data[[yname]]
    wm <- data$wm
    wt <- data$wt
    in_b <- if (is.null(by_var)) 1 else as.numeric(data[[by_var]] == by_level)
    Xa <- do.call("rbind", lapply(shifts, function(v) {
      nd <- data
      nd$logBLL <- a_fun(data$logBLL, v)
      model.matrix(tt, data = nd)
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

geex_shift <- function(formula, family, data, a_fun, by_var = NULL, by_level = NULL) {
  fit <- suppressWarnings(glm(formula, data = data, family = family, weights = wm))
  sub <- if (is.null(by_var)) rep(TRUE, nrow(data)) else data[[by_var]] == by_level
  mu <- vapply(shifts, function(v) {
    nd <- data[sub, ]
    nd$logBLL <- a_fun(nd$logBLL, v)
    weighted.mean(predict(fit, newdata = nd, type = "response"), data$wt[sub])
  }, numeric(1L))

  m <- geex::m_estimate(make_estfun(formula, family, a_fun, by_var, by_level),
                        data = data, roots = c(coef(fit), mu),
                        compute_roots = FALSE)

  k <- length(coef(fit)) + seq_along(shifts)
  list(est = mu, vcov = m@vcov[k, k])
}

fB <- Block12 ~ logBLL * Male + I(logBLL^2) + PIR
fitB <- glm(fB, data = d, family = quasibinomial, weights = wm)

# B1. Additive shift
gB1 <- geex_shift(fB, binomial(), d, function(a, v) a + v)
sB1 <- hc0(shift_curve(fitB, treat = "logBLL", shift = c(-.4, .4), n = 3, wts = d$wt))
stopifnot(isTRUE(all.equal(unname(gB1$est), unname(attr(sB1, ".est")), tolerance = 1e-8)))
record("B1 logistic, additive, model wts != target wts, HC0", "geex", se(gB1$vcov))
record("B1 logistic, additive, model wts != target wts, HC0", "shift_curve", se(attr(sB1, ".vcov")))
prog_tick(pr, label = "B1 geex additive")

# B2. A user-supplied proportional shift within the subgroup Male == 1 (logBLL
# takes negative values, so the built-in multiplicative shift, which requires a
# nonnegative treatment, is not used)
prop <- function(a, v) a * (1 + v)
gB2 <- geex_shift(fB, binomial(), d, prop, by_var = "Male", by_level = 1)
sB2 <- hc0(shift_curve(fitB, treat = "logBLL", shift = .4, shift_fun = prop,
                       n = 3, wts = d$wt, by = ~Male))
stopifnot(isTRUE(all.equal(unname(gB2$est), unname(attr(sB2, ".est"))[4:6], tolerance = 1e-8)))
record("B2 logistic, user-supplied proportional shift, subgroup Male = 1, HC0", "geex", se(gB2$vcov))
record("B2 logistic, user-supplied proportional shift, subgroup Male = 1, HC0", "shift_curve",
       se(attr(sB2, ".vcov"))[4:6])
prog_tick(pr, label = "B2 geex multiplicative")

# C. ADE from the interpolated curve vs. an exact central difference -----------

exact_ade <- function(fit, tt, data, eps = 1e-6) {
  Xp <- X_at(fit, tt, data, data$logBLL + eps)
  Xm <- X_at(fit, tt, data, data$logBLL - eps)
  linkinv <- family(fit)$linkinv
  mu.eta <- family(fit)$mu.eta
  ep <- drop(Xp %*% coef(fit))
  em <- drop(Xm %*% coef(fit))
  slope <- (linkinv(ep) - linkinv(em)) / (2 * eps)
  # gradient of the mean slope with respect to beta
  g <- colMeans(mu.eta(ep) * Xp - mu.eta(em) * Xm) / (2 * eps)
  list(est = mean(slope), slope = slope, g = g)
}

ade_accuracy <- list()
for (spec in list(list(name = "lm poly(3) * Male", fit = fitA),
                  list(name = "lm ns(5) * Male",
                       fit = lm(Math ~ splines::ns(logBLL, 5) * Male + PIR + Age, data = d)),
                  list(name = "logistic poly(3) * Male",
                       fit = glm(Block12 ~ poly(logBLL, 3) * Male + PIR, data = d,
                                 family = binomial)))) {
  fit <- spec$fit
  tt <- delete.response(terms(fit))
  ex <- exact_ade(fit, tt, d)

  # exact influence function of the ADE (HC2), for comparison of SEs
  Xf <- model.matrix(fit)
  hf <- hatvalues(fit)
  infl <- sandwich::estfun(fit) %*% sandwich::bread(fit) / sqrt(1 - hf)
  IF_ade <- (ex$slope - ex$est) + drop(infl %*% ex$g)
  se_exact <- sqrt(sum(IF_ade^2) / N^2)

  for (n in c(21L, 51L, 101L)) {
    a <- ade(shift_curve(fit, treat = "logBLL", n = n))
    ade_accuracy[[length(ade_accuracy) + 1L]] <- data.frame(
      model = spec$name, n = n,
      ade_exact = ex$est, ade = a$estimate,
      err_over_se = (a$estimate - ex$est) / se_exact,
      se_exact = se_exact, se = sqrt(diag(attr(a, ".vcov"))),
      se_ratio = sqrt(diag(attr(a, ".vcov"))) / se_exact)
  }
}
ade_accuracy <- do.call("rbind", ade_accuracy)
prog_tick(pr, label = "C ADE accuracy")

# D. Fractional weighted bootstrap, written independently ----------------------

# Each replicate refits the model with Dirichlet weights and averages the
# predictions with the same weights, as for the unconditional variance
ttA <- delete.response(terms(fitA))
Xs <- lapply(shifts, function(v) X_at(fitA, ttA, d, d$logBLL + v))
Xp <- X_at(fitA, ttA, d, d$logBLL + 1e-6)
Xm <- X_at(fitA, ttA, d, d$logBLL - 1e-6)

set.seed(20261116)
boot <- prog_lapply(pr, seq_len(n_boot), function(b) {
  w <- rexp(N)
  w <- w / mean(w)
  beta <- coef(lm.wfit(X, d$Math, w))
  c(vapply(Xs, function(Xa) weighted.mean(drop(Xa %*% beta), w), numeric(1L)),
    weighted.mean(drop((Xp - Xm) %*% beta) / 2e-6, w))
}, label = function(b) paste("boot", b))
boot <- do.call("rbind", boot)

aA <- ade(shift_curve(fitA, treat = "logBLL"))
sA_hc0 <- hc0(shift_curve(fitA, treat = "logBLL", shift = c(-.4, .4), n = 3))
aA_hc0 <- hc0(ade(shift_curve(fitA, treat = "logBLL")))

record("D lm, additive (bootstrap SEs: 3 shifts; ADE in se3 row below)", "bootstrap",
       apply(boot[, 1:3], 2, sd))
record("D lm, additive (bootstrap SEs: 3 shifts; ADE in se3 row below)", "shift_curve HC2",
       se(attr(sA, ".vcov")))
record("D lm, additive (bootstrap SEs: 3 shifts; ADE in se3 row below)", "shift_curve HC0",
       se(attr(sA_hc0, ".vcov")))
record("D lm, ADE", "bootstrap (exact ADE)", c(NA, NA, sd(boot[, 4])))
record("D lm, ADE", "ade() HC2", c(NA, NA, sqrt(attr(aA, ".vcov")[1, 1])))
record("D lm, ADE", "ade() HC0", c(NA, NA, sqrt(attr(aA_hc0, ".vcov")[1, 1])))
record("D lm, ADE", "marginaleffects::avg_slopes() (conditional)",
       c(NA, NA, marginaleffects::avg_slopes(fitA, variables = "logBLL")$std.error))

results <- do.call("rbind", results)

on.exit()
prog_end(pr, "done")

saveRDS(list(results = results, ade_accuracy = ade_accuracy), file.path(out_dir, "shift-results.rds"))
write.csv(results, file.path(out_dir, "shift-results.csv"), row.names = FALSE)

options(width = 150)
print(results, digits = 6)
print(ade_accuracy, digits = 4)
