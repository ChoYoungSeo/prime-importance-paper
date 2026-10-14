source("experiments/03_pme/00_common.R")

.libPaths(unique(c(
  normalizePath("r-library", mustWork = TRUE),
  file.path(Sys.getenv("HOME"), "R", "library"),
  .libPaths()
)))

benchmark_exact_dimension_budget <- 12L

first_selector <- function(X, y, seed) {
  set.seed(seed)
  importance <- first::first(
    X = as.data.frame(X),
    y = as.numeric(y),
    importance = "total",
    verbose = FALSE
  )
  which(is.finite(importance) & abs(importance) > 0)
}

lasso_pipeline_selector <- function(X, y, seed) {
  set.seed(seed)
  fit <- glmnet::cv.glmnet(as.matrix(X), as.numeric(y), alpha = 1)
  coefficient <- as.numeric(stats::coef(fit, s = "lambda.1se"))[-1L]
  which(coefficient != 0)
}

foci_pipeline_selector <- function(X, y, seed) {
  set.seed(seed)
  fit <- FOCI::foci(
    Y = as.numeric(y),
    X = as.matrix(X),
    num_features = NULL,
    stop = TRUE,
    numCores = 1,
    printIntermed = FALSE
  )
  if (is.null(fit$selectedVar$index)) integer() else
    sort(unique(as.integer(fit$selectedVar$index)))
}

ford_pipeline_selector <- function(X, y, seed) {
  set.seed(seed)
  fit <- FORD::ford(
    Y = as.numeric(y),
    X = as.matrix(X),
    stop = TRUE,
    standardize = "scale",
    numCores = 1,
    printIntermed = FALSE
  )
  if (is.null(fit$selectedVar$index)) integer() else
    sort(unique(as.integer(fit$selectedVar$index)))
}

ncmd_screen_value <- function(X, y, neighbors = 10L) {
  X <- standardize_predictors(X)
  y <- as.numeric(y)
  cache <- new.env(parent = emptyenv())
  value <- function(subset) {
    subset <- sort(unique(as.integer(subset)))
    if (!length(subset)) return(-Inf)
    key <- paste(subset, collapse = ",")
    if (!exists(key, cache, inherits = FALSE)) {
      nearest <- prime:::.prime_query_neighbors(
        X[, subset, drop = FALSE], seq_len(nrow(X)), neighbors
      )
      estimate <- abs(mean(y * rowMeans(matrix(
        y[nearest],
        nrow = nrow(X),
        ncol = neighbors
      ))))
      assign(key, estimate, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
  attr(value, "cache") <- cache
  value
}

ncmd_selector <- function(X, y, seed, neighbors = 10L) {
  set.seed(seed)
  p <- ncol(X)
  value <- ncmd_screen_value(X, y, neighbors)
  selected <- integer()
  current <- -Inf
  trace <- list()

  while (length(selected) < p) {
    candidates <- setdiff(seq_len(p), selected)
    candidate_values <- vapply(candidates, function(variable) {
      value(c(selected, variable))
    }, numeric(1))
    best_position <- which.max(candidate_values)
    best <- candidates[best_position]
    best_value <- candidate_values[best_position]
    if (!is.finite(best_value) || best_value < current) break
    trace[[length(trace) + 1L]] <- data.frame(
      step = length(trace) + 1L,
      added = best,
      value = best_value,
      gain = if (is.finite(current)) best_value - current else NA_real_,
      stringsAsFactors = FALSE
    )
    selected <- c(selected, best)
    current <- best_value
  }

  list(
    selected = as.integer(selected),
    trace = if (length(trace)) do.call(rbind, trace) else data.frame(),
    value_calls = length(ls(attr(value, "cache")))
  )
}

selected_pme_pipeline <- function(X, y, selected, p, estimator) {
  selected <- sort(unique(as.integer(selected)))
  if (!length(selected)) {
    return(list(
      effects = rep(NA_real_, p),
      selected = integer(),
      game = numeric(),
      value_calls = NA_integer_,
      feasible = FALSE,
      reason = "empty selected set"
    ))
  }
  if (length(selected) > benchmark_exact_dimension_budget) {
    return(list(
      effects = rep(NA_real_, p),
      selected = selected,
      game = numeric(),
      value_calls = NA_integer_,
      feasible = FALSE,
      reason = paste0(
        "selected dimension ", length(selected),
        " exceeds the benchmark exact-enumeration budget of ",
        benchmark_exact_dimension_budget
      )
    ))
  }

  value <- switch(
    estimator,
    local_linear = make_llr_value(X, y),
    pme_knn = make_local_constant_value(X, y, neighbors = 2L),
    ncmd = make_ncmd_value(X, y, neighbors = 5L),
    stop("Unknown PME estimator: ", estimator, call. = FALSE)
  )
  result <- estimate_pme_with_value(value, selected, p)
  result$feasible <- TRUE
  result$reason <- ""
  result
}

complete_pipeline_methods <- function() {
  c(
    "PRIME",
    "FIRST screen + PME-kNN",
    "NNVS screen + PME-kNN",
    "LASSO screen + PME-kNN",
    "FOCI screen + PME-kNN",
    "FORD screen + PME-kNN"
  )
}

run_complete_pipeline <- function(method, X, y, seed) {
  p <- ncol(X)
  selection_start <- proc.time()[["elapsed"]]
  selector_calls <- NA_integer_
  selected <- switch(
    method,
    "PRIME" = {
      prime_screen(X, y)$selected
    },
    "FIRST screen + PME-kNN" = first_selector(X, y, seed),
    "NNVS screen + PME-kNN" = {
      fit <- ncmd_selector(X, y, seed, neighbors = 10L)
      selector_calls <- fit$value_calls
      fit$selected
    },
    "LASSO screen + PME-kNN" = lasso_pipeline_selector(X, y, seed),
    "FOCI screen + PME-kNN" = foci_pipeline_selector(X, y, seed),
    "FORD screen + PME-kNN" = ford_pipeline_selector(X, y, seed),
    stop("Unknown pipeline: ", method, call. = FALSE)
  )
  selection_elapsed <- proc.time()[["elapsed"]] - selection_start

  allocation_start <- proc.time()[["elapsed"]]
  estimator <- if (method == "PRIME") {
    "local_linear"
  } else {
    "pme_knn"
  }
  error <- ""
  result <- tryCatch(
    selected_pme_pipeline(X, y, selected, p, estimator),
    error = function(condition) {
      error <<- conditionMessage(condition)
      list(
        effects = rep(NA_real_, p),
        selected = selected,
        value_calls = NA_integer_,
        feasible = FALSE,
        reason = error
      )
    }
  )
  allocation_elapsed <- proc.time()[["elapsed"]] - allocation_start
  if (!nzchar(error) && !result$feasible) error <- result$reason

  list(
    effects = result$effects,
    selected = sort(unique(as.integer(selected))),
    selector_calls = selector_calls,
    value_calls = result$value_calls,
    selection_elapsed = selection_elapsed,
    allocation_elapsed = allocation_elapsed,
    elapsed_sec = selection_elapsed + allocation_elapsed,
    feasible = result$feasible,
    error = error
  )
}

rmvnorm_eigen <- function(n, mean, covariance) {
  mean <- as.numeric(mean)
  dimension <- length(mean)
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  tolerance <- max(1, max(abs(eig$values))) * 1e-12
  if (min(eig$values) < -tolerance) {
    stop("The covariance matrix is not positive semidefinite.", call. = FALSE)
  }
  root <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0)), nrow = dimension)
  draws <- matrix(stats::rnorm(n * dimension), nrow = n) %*% t(root)
  sweep(draws, 2L, mean, "+")
}

simulate_correlated_pipeline <- function(
  benchmark,
  n,
  p,
  rho,
  seed,
  noise_ratio = 0.10
) {
  specification <- proxy_registry()[[benchmark]]
  active <- specification$active
  proxies <- max(active) + seq_along(active)
  if (p < max(proxies)) stop("p is too small for the proxy design.")
  covariance <- diag(p)
  for (index in seq_along(active)) {
    covariance[active[index], proxies[index]] <- rho
    covariance[proxies[index], active[index]] <- rho
  }
  set.seed(seed)
  X <- stats::pnorm(rmvnorm_eigen(n, numeric(p), covariance))
  signal <- specification$fun(X[, active, drop = FALSE])
  noise_sd <- noise_ratio * stats::sd(signal)
  y <- signal + stats::rnorm(n, sd = noise_sd)
  colnames(X) <- paste0("X", seq_len(p))
  list(
    X = X,
    y = y,
    signal = signal,
    active = active,
    proxies = proxies,
    p = p,
    noise_sd = noise_sd
  )
}

correlated_pipeline_seed <- function(benchmark, n, replication) {
  benchmark_id <- match(benchmark, names(proxy_registry()))
  202608040L + 100000L * benchmark_id + 10000000L + n + replication
}

scaling_seed <- function(benchmark, p, replication) {
  benchmark_id <- match(
    benchmark,
    c("branin", "friedman", "otl", "piston", "ishigami", "borehole")
  )
  if (is.na(benchmark_id)) stop("Unknown scaling benchmark: ", benchmark)
  20260803L + 1000000L * benchmark_id + 1000L * p + replication
}
