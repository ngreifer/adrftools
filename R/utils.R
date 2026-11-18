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
# conservative
.mvt_df <- function(df) {
  if (is.finite(df)) max(1, floor(df))
  else Inf
}

# Simultaneous critical values and p-values for Wald inference, from the
# multivariate t (or normal) distribution with correlation matrix `corr`. Both
# use the same settings of the Genz-Bretz algorithm, so intervals from
# summary() and bands from plot() agree to within its Monte Carlo error. With
# 51 highly correlated estimates, these settings give critical values within
# about 1e-3 and p-values within about 1e-4 of their exact values.
.mvt_maxpts <- 1e5
.mvt_abseps <- 1e-5

.mvt_crit <- function(conf_level, corr, df = Inf) {
  t_crit <- try(mvtnorm::qmvt(conf_level,
                              tail = "both.tails",
                              df = .mvt_df(df), corr = corr,
                              maxpts = .mvt_maxpts,
                              abseps = .mvt_abseps)$quantile,
                silent = TRUE)

  if (null_or_error(t_crit)) {
    arg::err(c("There was an error computing simultaneous confidence intervals.",
               x = if (inherits(t_crit, "try-error")) conditionMessage(.attr(t_crit, "condition"))))
  }

  t_crit
}

# P(max_j |T_j| >= |z|) for each value of `z`. The p-value is a decreasing
# function of |z|, so it is computed once per distinct value and made
# monotone (taking the larger, i.e., conservative, value when Monte Carlo error
# reverses the order). It lies between the pointwise p-value and its
# Bonferroni bound, so when these are within the precision of the algorithm,
# the bound is used without calling pmvt(); this skips the estimates far from
# the null, which are the slowest to compute accurately.
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

  for (i in todo) {
    pr <- try(1 - mvtnorm::pmvt(lower = rep.int(-u[i], k),
                                upper = rep.int(u[i], k),
                                df = df, corr = corr,
                                maxpts = .mvt_maxpts,
                                abseps = .mvt_abseps),
              silent = TRUE)

    if (null_or_error(pr)) {
      arg::err("there was an error computing simultaneous p-values")
    }

    p[i] <- min(max(pr, p_point[i]), p_bonf[i])
  }

  p <- rev(cummax(rev(p)))

  p[fmatch(az, u)]
}
