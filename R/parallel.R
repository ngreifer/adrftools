# Parallel processing through the future framework, following the r-parallel
# skill's assets/parallel.R.
#
# The interface to the user is `future::plan()` and nothing else. A function
# that parallelizes calls `par_lapply()` where it would have called `lapply()`,
# with the random number streams drawn by `parallel_streams()` first.
#
# None of the units draw random numbers, so the estimates do not depend on the
# streams. Some functions they call start the session's random number stream
# when there is none, as on a new worker, without drawing from it (e.g.,
# collapse::qG() and the predictions of a *bartisan* fit), which *future* takes
# for drawing random numbers when no seeds are given. Giving each unit its own
# stream, as for units that do draw random numbers, keeps the check quiet and
# keeps the results independent of the plan if a unit ever does draw one.

# Whether to send work to workers: future and future.apply are installed and
# the plan has more than one worker. Without either package, or with a
# sequential plan, the same units run one after another in this session, which
# is a question of speed rather than of whether the function works.
use_future <- function() {
  rlang::is_installed(c("future", "future.apply")) &&
    future::nbrOfWorkers() > 1L
}

# One L'Ecuyer-CMRG stream per unit, advanced from the session's current seed,
# so that one `set.seed()` before the call reproduces the whole run whatever
# RNG kind the session uses: switching the kind to L'Ecuyer seeds it from a draw
# of the current generator, which is a deterministic function of the seed.
#
# The streams belong to the units, not to the workers, which is what makes the
# result independent of how many workers there are and of which worker ran
# which unit.
#
# The session's own stream is restored and then moved on by one draw. Restoring
# it exactly would hand the next call the same streams, so two calls in a row
# without a `set.seed()` between them would return identical results.
parallel_streams <- function(n) {
  old <- restore_stream()

  on.exit({
    old()
    stats::runif(1L)
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")

  seed <- get(".Random.seed", envir = globalenv())
  out <- vector("list", n)

  for (i in seq_len(n)) {
    out[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }

  out
}

# Captures the session's RNG state and returns the function that puts it back,
# including the kind, so that switching to L'Ecuyer does not leave the session
# on a generator it did not choose.
restore_stream <- function() {
  kind <- RNGkind()

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    return(function() {
      RNGkind(kind[1L], kind[2L], kind[3L])
    })
  }

  seed <- get(".Random.seed", envir = globalenv())

  function() {
    RNGkind(kind[1L], kind[2L], kind[3L])
    assign(".Random.seed", seed, envir = globalenv())
  }
}

# `lapply(X, FUN, ...)`, on the workers of the current plan when there are
# several and in this session otherwise, with the same result either way.
#
# `seeds` is `NULL` when `FUN` draws no random numbers, and otherwise the output
# of `parallel_streams(length(X))`, drawn by the caller before this call. Each
# unit runs on its own stream in both branches, and the caller's stream is put
# back afterwards in both: `future_lapply()` given a list of seeds moves the
# caller's `.Random.seed` itself, and the sequential branch leaves it on the
# last unit's stream, so without the restore the code after the call would draw
# different numbers depending on the plan.
#
# `chunk_size = 1L` is dynamic scheduling: each unit is its own future, and a
# worker is handed the next one as soon as it finishes. `chunk_size = NULL`
# cuts `X` into one piece per worker up front, which sends the globals once per
# worker rather than once per unit.
#
# `FUN` should be a function defined at the top level of the package, with the
# per-unit data passed as the elements of `X`, so that a worker receives its
# own unit and no more. `packages` are attached on the workers in addition to
# this package; see .worker_packages(). `local = TRUE` runs the units in this
# session whatever the plan (e.g., for a model that cannot be sent to the
# workers; see .prepare_for_workers()), on the same streams.
par_lapply <- function(X, FUN, ..., seeds = NULL, chunk_size = 1L,
                       packages = NULL, local = FALSE) {
  if (!is_null(seeds)) {
    restore <- restore_stream()
    on.exit(restore(), add = TRUE)
  }

  if (!local && use_future()) {
    return(future.apply::future_lapply(X, FUN, ...,
                                       future.seed = seeds %or% FALSE,
                                       future.chunk.size = chunk_size,
                                       future.packages = unique(c("adrftools", packages))))
  }

  if (is_null(seeds)) {
    return(lapply(X, FUN, ...))
  }

  # Indexed rather than over `X` so each unit can be given its stream, and
  # named afterwards so the result matches `lapply()`'s and `future_lapply()`'s.
  lapply(seq_along(X), function(i) {
    assign(".Random.seed", seeds[[i]], envir = globalenv())
    FUN(X[[i]], ...)
  }) |>
    stats::setNames(names(X))
}

# The packages a worker needs beyond this one to predict from a model `x`. A
# model formula is evaluated again to predict and can call a function from an
# attached package (e.g., `ns()` after `library(splines)`), so the packages
# attached in this session are attached on the workers too, in the same order.
# The packages that define S3 methods for the classes of the model are loaded
# as well, because a worker that has not loaded them would not find
# predict() or marginaleffects::get_predict() methods for the model; *future*
# finds the packages of the functions it sends, not those of the objects. A
# registered method is stored with its generic rather than with the package
# that defines it, so the methods are looked for in the S3 tables of every
# loaded namespace, and their packages read from their environments.
.worker_packages <- function(x) {
  attached <- search() |>
    grep(pattern = "^package:", value = TRUE) |>
    sub(pattern = "^package:", replacement = "") |>
    setdiff(c("base", "stats", "graphics", "grDevices", "utils",
              "datasets", "methods")) |>
    rev()

  suffixes <- paste0(".", class(x))

  defining <- unlist(lapply(loadedNamespaces(), function(ns) {
    tbl <- asNamespace(ns)[[".__S3MethodsTable__."]]

    if (is_null(tbl)) {
      return(NULL)
    }

    methods <- names(tbl)
    methods <- methods[Reduce("|", lapply(suffixes, endsWith, x = methods))]

    vapply(methods, function(m) {
      e <- environment(get(m, envir = tbl))

      if (isNamespace(e)) getNamespaceName(e)
      else NA_character_
    }, character(1L), USE.NAMES = FALSE)
  }), use.names = FALSE)

  # The attached packages come last so that they are attached last and so
  # mask the others, as they do here
  defining <- defining |>
    na_rm() |>
    unique() |>
    setdiff(c("base", attached))

  c(defining, attached)
}

# Whether a model can be sent to the workers of the plan and predict there as
# it does here, readying it to be sent if it needs that. A *stan4bart* fit predicts through a pointer to trees
# held outside R, which is lost when the fit is copied, so its predictions
# stay in this session (unless the workers are forked, and so share its
# memory). A *dbarts* fit holds its trees outside R as well, but stores them in
# the R object once its state is asked for, which *dbarts* requires before a
# fit is saved; without that, a copy predicts from different trees, with no
# error.
.prepare_for_workers <- function(x) {
  if (!use_future()) {
    return(TRUE)
  }

  if (inherits(x, "stan4bartFit")) {
    return(inherits(future::plan(), "multicore"))
  }

  if (inherits(x, "bart") && !is_null(x[["fit"]])) {
    invisible(x[["fit"]]$state)
  }

  TRUE
}
