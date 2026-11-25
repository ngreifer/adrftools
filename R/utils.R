#Note: these utils may use collapse and cli

#Strings
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (rlang::is_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (length(quotes) != 1L || !is.numeric(quotes) || (quotes %!in% c(1, 2))) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}
space <- function(n) {
  strrep(" ", n)
}
txtbar <- function(n) {
  strrep("\u2500", n)
}
str_rev <- function(x) {
  strsplit(x, NULL) |>
    lapply(rev) |>
    vapply(paste, character(1L), collapse = "")
}
center_just <- function(x, wrt = NULL) {
  if (is_null(wrt)) {
    n <- getOption("width")
  }
  else {
    n <- max(nchar(as.character(wrt)))
  }

  paste0(space(max(0, floor((n - nchar(x)) / 2))), x)
}
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }

  nas <- is.na(x)
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) {
    x <- as.character(x)
  }

  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }

  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1L))
  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  # A column whose values would all round to 0 gets enough digits to show
  # them (e.g., estimates of an outcome on a very small scale)
  col_digits <- rep.int(digits, length(df))

  for (i in which(nums)) {
    fin <- df[[i]][is.finite(df[[i]])]
    m <- if (length(fin) > 0L) max(abs(fin)) else 0

    if (m > 0 && m < .5 * 10^-digits) {
      col_digits[i] <- digits - 1 + ceiling(-log10(m))
    }
  }

  for (i in which(nums)) {
    o.negs[, i] <- !nas[, i] & df[[i]] < 0 & round(df[[i]], col_digits[i]) == 0
    df[[i]] <- round(df[[i]], digits = col_digits[i])
  }

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      s_lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[s_lengths > 1L] <- nchar(vapply(s[s_lengths > 1L], `[[`, character(1L), 2L))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(s_lengths > 1L, "", if (nzchar(as.character(pad))) "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) strrep(pad, n), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (!is_null(rn)) rownames(df) <- rn
  if (!is_null(cn)) names(df) <- cn

  df
}

#Numbers
# Tolerance for treating a treatment value as a point on the grid `v`,
# relative to the width of the grid so that it does not depend on the scale of
# the treatment
.grid_tol <- function(v) {
  width <- diff1(.range(v))

  sqrt(.Machine$double.eps) * {
    if (width > 0) width
    else max(abs(v), 1)
  }
}

check_if_zero <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  # this is the default tolerance used in all.equal
  abs(x) < tolerance
}
squish <- function(p, lo = 1e-6, hi = 1 - lo) {
  if (lo > -Inf)
    p[p < lo] <- lo

  if (hi < Inf)
    p[p > hi] <- hi

  p
}

#Statistics
# Weighted standard deviation as computed by cobalt::col_w_sd(): the weights
# are normalized to sum to 1 and the variance uses the denominator
# 1 - sum(w^2), so it does not depend on the scale of the weights and equals
# sd() when the weights are equal. (collapse::fvar() treats the weights as
# frequencies, so its result depends on their scale.)
.w_sd <- function(x, w = NULL) {
  if (is_null(w)) {
    return(stats::sd(x))
  }

  w <- w / sum(w)

  xc <- x - sum(w * x)

  sqrt(sum(w * xc^2) / (1 - sum(w^2)))
}

bw.nrd <- function(x) {
  #R's bw.nrd doesn't always work, but bw.nrd0 does
  stats::bw.nrd0(x) * 1.06 / .9
}

#Uniqueness
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na_rm(x)

    if (!na.rm) {
      return(is_null(x))
    }
  }

  if (length(x) == 1L) {
    return(TRUE)
  }

  if (is.numeric(x)) {
    return(check_if_zero(diff1(.range(x))))
  }

  allv(x, x[1L])
}

#R Processing
make_list <- function(n) {
  if (is_null(n)) {
    vector("list", 0L)
  }
  else if (length(n) == 1L && is.numeric(n) && n >= 0) {
    vector("list", as.integer(n))
  }
  else if (is.atomic(n)) {
    setNames(vector("list", length(n)),
             as.character(n))
  }
  else {
    stop("'n' must be an integer(ish) scalar or an atomic variable.")
  }
}
make_df <- function(ncol, nrow = 0L, types = "numeric") {
  if (is_null(ncol)) {
    ncol <- 0L
  }

  if (length(ncol) == 1L && is.numeric(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }

  if (is_null(nrow)) {
    nrow <- 0L
  }

  if (length(nrow) == 1L && is.numeric(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }

  df <- matrix(NA_real_, nrow = nrow, ncol = ncol) |>
    as.data.frame.matrix()

  names(df) <- col_names
  rownames(df) <- row_names

  if (is_null(types)) {
    return(df)
  }

  if (!any(c(1L, ncol) == length(types))) {
    stop("'types' must be equal to the number of columns.")
  }

  if (!is.character(types) ||
      !all(types %in% c("numeric", "integer", "logical", "character", NA))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }

  if (length(types) == 1L) {
    types <- rep.int(types, ncol)
  }

  for (i in which(!is.na(types))) {
    if (types[i] != "numeric") {
      df[[i]] <- get(types[i])(nrow)
    }
  }

  df
}
rep_with <- function(x, y) {
  #Helper function to fill named vectors with x and given names of y
  setNames(rep.int(x, length(y)), names(y))
}
is_null <- function(x) {length(x) == 0L}
is_number <- function(x) {
  length(x) == 1L && is.numeric(x) && !anyNA(x)
}
is_count <- function(x) {
  is_number(x) && x >= 0 && check_if_zero(x - trunc(x))
}
`%or%` <- function(x, y) {
  # like `%||%` but works for non-NULL length 0 objects
  if (is_null(x)) y else x
}
null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}
.attr <- function(x, which, exact = TRUE) {
  attr(x, which, exact = exact)
}
block_diag <- function(m, n = 1L) {
  arg::arg_count(n)

  m <- as.matrix(m)

  if (n == 1L || is_null(m)) {
    return(m)
  }

  nr <- nrow(m)
  nc <- ncol(m)

  out <- matrix(0, nrow = nr * n, ncol = nc * n)

  for (i in seq_len(n)) {
    out[(i - 1L) * nr + seq_len(nr), (i - 1L) * nc + seq_len(nc)] <- m
  }

  out
}
quad_mult <- function(A, B) {
  tcrossprod(A, tcrossprod(A, B))
}
get_varnames <- function(expr) {
  recurse <- function(e) {
    if (is.symbol(e)) {
      # bare variable like age
      return(as.character(e))
    }

    if (!is.call(e)) {
      return(NULL)
    }

    # keep as-is for $, [[, and [
    fn <- e[[1L]]
    if (fn == as.name("$") || fn == as.name("[[") || fn == as.name("[")) {
      return(deparse1(e))
    }

    # strip outer function, recurse into arguments
    lapply(as.list(e)[-1L], recurse) |>
      unlist()
  }

  recurse(expr)
}

diff1 <- function(x) {
  x[-1L] - x[-length(x)]
}

#Extract variables from ..., similar to ...elt(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      ...elt(.m1[1L]) %or% .ifnotfound
    }
  })

  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}
...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))

  not_found <- is.na(found)

  if (all(not_found)) {
    return(list())
  }

  lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }) |>
    setNames(x[!not_found])
}

any_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (isTRUE(FUN(x, ...))) {
      return(TRUE)
    }
  }

  FALSE
}
all_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (!isTRUE(FUN(x, ...))) {
      return(FALSE)
    }
  }

  TRUE
}

# cli utilities
.it <- function(...) cli::style_italic(...)

# mvtnorm::pmvt() and mvtnorm::qmvt() require integer degrees of freedom (and
# treat 0 as infinite), so finite degrees of freedom are rounded down, which is
# conservative; the simulated critical values use the same rounding so that
# simultaneous intervals and p-values agree
.mvt_df <- function(df) {
  if (is.finite(df)) max(1, floor(df))
  else Inf
}

# Settings of the Genz-Bretz algorithm used by mvtnorm::pmvt() for simultaneous
# p-values and by mvtnorm::qmvt() for critical values of a few estimates. With
# 51 highly correlated estimates, these give p-values within about 1e-4 of
# their exact values.
.mvt_maxpts <- 1e5
.mvt_abseps <- 1e-5

# Settings of the simulation used for simultaneous critical values: directions
# are drawn until the Monte Carlo standard error of the critical value is at
# most `.sim_crit_se` or `.sim_max_draws` have been drawn, starting with
# `.sim_start`, in chunks of `.sim_chunk` so that little memory is used at once
.sim_crit_se <- 1e-3
.sim_start <- 2e4
# Settings of the simulation used for simultaneous p-values (see .sim_tail()):
# the number of directions drawn first, and the target Monte Carlo standard
# error, max(.sim_p_se_abs, .sim_p_se_rel * p) for a p-value p
.sim_start_p <- 1e5
.sim_p_se_abs <- 1e-4
.sim_p_se_rel <- .005

# The largest number of estimates for which mvtnorm::qmvt() is used when their
# correlation matrix has full rank (see .mvt_crit())
.qmvt_max_dim <- 10
.sim_max_draws <- 1e6
.sim_chunk <- 1e3

# The simultaneous critical value for Wald inference, the `conf_level`
# quantile of max_j |T_j| for T multivariate t (or normal) with correlation
# matrix `corr`.
#
# The estimates along effect curves are linear in a few model coefficients, so
# `corr` usually has a rank r far below its dimension k (e.g., 12 for 306
# estimates), and T can be written as L u / sqrt(W / df) with L a k x r factor
# of `corr`, u standard normal in r dimensions, and W chi-squared with df
# degrees of freedom. Writing u = R theta, with theta uniform on the unit
# sphere and R^2 chi-squared with r degrees of freedom independent of it,
# max_j |T_j| = R m(theta) / sqrt(W / df) with m(theta) = max_j |l_j' theta|,
# so given theta, the probability that it exceeds c is that of an F(r, df)
# variable exceeding c^2 / (r m(theta)^2) (or of a chi-squared variable with r
# degrees of freedom exceeding c^2 / m(theta)^2 when df is infinite). The tail
# probability is the mean of these over random directions theta, and the
# critical value is the root of its difference from 1 - conf_level. This is
# exact apart from Monte Carlo error, which is much smaller than that of
# simulating T directly because R and W are integrated analytically. It is
# also much faster than numerical integration over all k estimates: for 306
# estimates, mvtnorm::qmvt() took 18 s with a standard error of about .005,
# and this takes under 1 s with a standard error of .001.
#
# For a few estimates whose correlation matrix has full rank (e.g., when
# summary() is applied to an effect curve evaluated at three points),
# mvtnorm::qmvt() is as fast and more accurate (a standard error of about
# 2e-5 for three estimates and 1e-4 for 10 equicorrelated ones), so it is used
# instead. For rank-deficient matrices, such as those of the estimates at the
# grid points of a curve, qmvt() is much less accurate (a standard error of
# .002 to .006 even for five estimates), and for full-rank matrices of more
# than 10 estimates, it is no more accurate than the simulation and slower.
#
# Eigenvalues below 1e-6 of the largest are dropped, which removes the ridge
# .to_psd() adds to rank-deficient matrices. With a single dimension, every T_j
# is the same t variable up to sign, so the pointwise critical value is exact.
.mvt_crit <- function(conf_level, corr, df = Inf) {
  df <- .mvt_df(df)
  alpha <- 1 - conf_level

  .q <- function(p) {
    if (is.finite(df)) stats::qt(p, df)
    else stats::qnorm(p)
  }

  L <- .corr_factor(corr)
  r <- fncol(L)

  # The simultaneous critical value lies between the pointwise and Bonferroni
  # critical values
  lower <- .q(1 - alpha / 2)

  if (r == 1L) {
    return(lower)
  }

  if (r == fnrow(L) && r <= .qmvt_max_dim) {
    return(.qmvt_crit(conf_level, corr, df))
  }

  upper <- .q(1 - alpha / (2 * fnrow(L)))

  m <- .sim_sphere_max(L, .sim_start)

  repeat {
    crit <- stats::uniroot(function(cr) mean(.tail_given_dir(cr, m, r, df)) - alpha,
                           interval = c(lower, upper), extendInt = "downX",
                           tol = 1e-8)$root

    se <- .crit_se(crit, m, r, df)

    if (se <= .sim_crit_se || length(m) >= .sim_max_draws) {
      break
    }

    # Enough further directions to reach the target, with a margin
    n_total <- min(.sim_max_draws, ceiling(1.1 * length(m) * (se / .sim_crit_se)^2))

    m <- c(m, .sim_sphere_max(L, n_total - length(m)))
  }

  crit
}

.qmvt_crit <- function(conf_level, corr, df = Inf) {
  t_crit <- try(mvtnorm::qmvt(conf_level,
                              tail = "both.tails",
                              df = if (is.finite(df)) df else 0,
                              corr = corr,
                              maxpts = .mvt_maxpts,
                              abseps = .mvt_abseps)$quantile,
                silent = TRUE)

  if (null_or_error(t_crit)) {
    arg::err(c("There was an error computing simultaneous confidence intervals.",
               x = if (inherits(t_crit, "try-error")) conditionMessage(.attr(t_crit, "condition"))))
  }

  t_crit
}

# A factor L of a correlation matrix with tcrossprod(L) equal to it after
# dropping eigenvalues below `rel_tol` of the largest, with its rows rescaled
# to unit length so that each implied variance is 1
.corr_factor <- function(corr, rel_tol = 1e-6) {
  es <- eigen(corr, symmetric = TRUE)

  keep <- es$values > rel_tol * es$values[1L]

  L <- es$vectors[, keep, drop = FALSE] %*%
    diag(sqrt(es$values[keep]), nrow = sum(keep))

  L / sqrt(rowSums(L^2))
}

# `n` draws of m(theta) = max_j |l_j' theta| for directions theta uniform on
# the unit sphere, made in chunks so that only one chunk of L theta is held at
# once
.sim_sphere_max <- function(L, n) {
  r <- fncol(L)
  out <- numeric(n)
  done <- 0L

  while (done < n) {
    b <- min(.sim_chunk, n - done)

    U <- matrix(stats::rnorm(r * b), nrow = r, ncol = b)
    U <- U / rep(sqrt(colSums(U^2)), each = r)

    out[done + seq_len(b)] <- (L %*% U) |>
      abs() |>
      fmax()

    done <- done + b
  }

  out
}

# P(max_j |T_j| > crit) given each direction, with `m` the values of m(theta)
.tail_given_dir <- function(crit, m, r, df = Inf) {
  if (is.finite(df)) stats::pf(crit^2 / (r * m^2), r, df, lower.tail = FALSE)
  else stats::pchisq(crit^2 / m^2, r, lower.tail = FALSE)
}

# The Monte Carlo standard error of the critical value: that of the estimated
# tail probability divided by the derivative of the tail probability
.crit_se <- function(crit, m, r, df = Inf) {
  se_p <- stats::sd(.tail_given_dir(crit, m, r, df)) / sqrt(length(m))

  d_p <- {
    if (is.finite(df)) mean(stats::df(crit^2 / (r * m^2), r, df) * 2 * crit / (r * m^2))
    else mean(stats::dchisq(crit^2 / m^2, r) * 2 * crit / m^2)
  }

  se_p / d_p
}

# P(max_j |T_j| >= |z|) for each value of `z`. The p-value is a decreasing
# function of |z|, so it is computed once per distinct value and made
# monotone. It lies between the pointwise p-value and its Bonferroni bound, so
# when these are within 1e-4, the bound is used; this skips the estimates far
# from the null. When the correlation matrix is singular (e.g., for estimates
# along curves, which are smooth functions of a few model coefficients), the
# others are computed by simulation from a factor of the correlation matrix
# (see .sim_tail()), which was 2 to 10 times more accurate than
# mvtnorm::pmvt() for p-values below .05 and 10 to 40 times faster (e.g., 0.3 s
# rather than 10 s for six p-values of 306 estimates). When it has full rank,
# pmvt() is used: it is nearly exact for a few estimates (an error of about
# 1e-6 for three or four) and was more accurate than the simulation for
# p-values above .1 with 50 equicorrelated estimates.
.mvt_p_values <- function(z, corr, df = Inf) {
  k <- fnrow(corr)
  df <- .mvt_df(df)

  az <- abs(z)
  u <- funique(az, sort = TRUE)

  p_point <- {
    if (is.finite(df)) 2 * stats::pt(-u, df)
    else 2 * stats::pnorm(-u)
  }

  p_bonf <- pmin(1, k * p_point)

  p <- p_bonf

  todo <- which(p_bonf - p_point > 1e-4)

  if (!is_null(todo)) {
    L <- .corr_factor(corr)

    p[todo] <- {
      if (fncol(L) == k) .pmvt_tail(u[todo], corr, df)
      else .sim_tail(u[todo], L, df)
    }

    p[todo] <- pmin(pmax(p[todo], p_point[todo]), p_bonf[todo])
  }

  p <- rev(cummax(rev(p)))

  p[fmatch(az, u)]
}

# P(max_j |T_j| >= u) for each value of `u` using mvtnorm::pmvt()
.pmvt_tail <- function(u, corr, df = Inf) {
  k <- fnrow(corr)

  vapply(u, function(.u) {
    pr <- try(1 - mvtnorm::pmvt(lower = rep.int(-.u, k),
                                upper = rep.int(.u, k),
                                df = if (is.finite(df)) df else 0,
                                corr = corr,
                                maxpts = .mvt_maxpts,
                                abseps = .mvt_abseps),
              silent = TRUE)

    if (null_or_error(pr)) {
      arg::err("there was an error computing simultaneous p-values")
    }

    pr
  }, numeric(1L))
}

# P(max_j |T_j| >= u) for each value of `u` by simulation, with L a factor of
# the correlation matrix from .corr_factor(): the mean over random directions
# of the probability given each direction, as in .mvt_crit(). The same
# directions are used for every value of `u`. Directions are added until the
# Monte Carlo standard error of every p-value is at most
# max(.sim_p_se_abs, .sim_p_se_rel * p), so small p-values are estimated with
# small absolute error and large ones (whose exact value matters little) are
# not estimated more precisely than needed, or until `.sim_max_draws` have
# been drawn.
.sim_tail <- function(u, L, df = Inf) {
  r <- fncol(L)

  m <- .sim_sphere_max(L, .sim_start_p)

  repeat {
    n <- length(m)

    stats_u <- vapply(u, function(.u) {
      g <- .tail_given_dir(.u, m, r, df)
      c(mean(g), stats::sd(g) / sqrt(n))
    }, numeric(2L))

    p <- stats_u[1L, ]
    se <- stats_u[2L, ]
    target <- pmax(.sim_p_se_abs, .sim_p_se_rel * p)

    if (all(se <= target) || n >= .sim_max_draws) {
      break
    }

    # Enough further directions for the p-value farthest from its target,
    # with a margin
    n_total <- min(.sim_max_draws, ceiling(1.1 * n * max((se / target)^2)))

    m <- c(m, .sim_sphere_max(L, n_total - n))
  }

  p
}

