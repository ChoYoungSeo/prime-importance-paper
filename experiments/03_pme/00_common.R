source("experiments/00_reference/benchmarks.R")

standardize_predictors <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X)
  column_scale <- apply(X, 2L, stats::sd)
  column_scale[!is.finite(column_scale) | column_scale == 0] <- 1
  sweep(sweep(X, 2L, center, "-"), 2L, column_scale, "/")
}

coalition_subset <- function(mask, variables, removed = FALSE) {
  present <- vapply(seq_along(variables), function(position) {
    bitwAnd(mask, bitwShiftL(1L, position - 1L)) != 0L
  }, logical(1))
  variables[if (removed) !present else present]
}

make_local_constant_value <- function(X, y, neighbors = 2L) {
  X <- standardize_predictors(X)
  y <- as.numeric(y)
  variance_y <- stats::var(y)
  cache <- new.env(parent = emptyenv())
  value <- function(subset) {
    subset <- sort(unique(as.integer(subset)))
    if (!length(subset)) return(0)
    key <- paste(subset, collapse = ",")
    if (!exists(key, cache, inherits = FALSE)) {
      projected <- X[, subset, drop = FALSE]
      if (neighbors == 1L) {
        indices <- matrix(seq_len(nrow(X)), ncol = 1L)
      } else {
        nearest <- prime:::.prime_query_neighbors(
          projected, seq_len(nrow(X)), neighbors - 1L
        )
        indices <- cbind(seq_len(nrow(X)), nearest)
      }
      local_variance <- apply(
        matrix(y[indices], nrow = nrow(indices)),
        1L,
        stats::var
      )
      estimate <- 1 - mean(local_variance) / variance_y
      assign(key, estimate, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
  attr(value, "cache") <- cache
  value
}

make_ncmd_value <- function(X, y, neighbors = 5L) {
  X <- standardize_predictors(X)
  y <- as.numeric(y)
  n <- length(y)
  denominator <- stats::var(y)
  cross_mean <- (sum(y)^2 - sum(y^2)) / (n * (n - 1))
  cache <- new.env(parent = emptyenv())
  value <- function(subset) {
    subset <- sort(unique(as.integer(subset)))
    if (!length(subset)) return(0)
    key <- paste(subset, collapse = ",")
    if (!exists(key, cache, inherits = FALSE)) {
      projected <- X[, subset, drop = FALSE]
      nearest <- prime:::.prime_query_neighbors(
        projected, seq_len(n), neighbors
      )
      neighbor_term <- mean(y * rowMeans(matrix(
        y[nearest],
        nrow = n,
        ncol = neighbors
      )))
      estimate <- (neighbor_term - cross_mean) / denominator
      assign(key, estimate, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
  attr(value, "cache") <- cache
  value
}

make_llr_value <- function(X, y) {
  prepared <- prime:::.prime_prepare(
    X, y,
    rescale = TRUE,
    twin_r = NULL,
    lambda = 1e-6
  )
  value <- function(subset) {
    subset <- sort(unique(as.integer(subset)))
    if (!length(subset)) return(0)
    k <- prime:::.prime_neighborhood_size(
      prepared$n, length(subset)
    )
    prime:::.prime_value(prepared, subset, k)
  }
  attr(value, "cache") <- prepared$value_cache
  # Local-linear allocations use the public noise-adjusted risk constructor.
  # Keep the subset-value interface for diagnostics and the other estimators.
  attr(value, "local_linear_data") <- prepared
  value
}

estimate_pme_with_value <- function(value, variables, p) {
  variables <- sort(unique(as.integer(variables)))
  s <- length(variables)
  if (!s) {
    stop("empty_selected_set: allocation unavailable", call. = FALSE)
  }
  full_mask <- 2^s - 1L
  q_full <- value(variables)
  if (!is.finite(q_full) || q_full <= sqrt(.Machine$double.eps)) {
    stop("The full retained set has no positive estimated signal.")
  }
  local_linear_data <- attr(value, "local_linear_data", exact = TRUE)
  game <- if (!is.null(local_linear_data)) {
    prime:::.prime_estimate_subset_values(local_linear_data, variables)$raw
  } else {
    vapply(0:full_mask, function(mask) {
      keep <- coalition_subset(mask, variables, removed = TRUE)
      (q_full - value(keep)) / q_full
    }, numeric(1))
  }
  game[1L] <- 0
  game[length(game)] <- 1
  local_effects <- prime_allocate(game)
  effects <- numeric(p)
  effects[variables] <- as.numeric(local_effects)
  list(
    effects = effects,
    selected = variables,
    game = game,
    projected_values = attr(local_effects, "projected_values"),
    value_calls = length(ls(attr(value, "cache"))) + 1L
  )
}

proxy_registry <- function() {
  list(
    branin = list(
      label = "Branin",
      active = 1:2,
      p = 8L,
      fun = branin_function
    ),
    bratley = list(
      label = "Bratley",
      active = 1:5,
      p = 10L,
      fun = bratley_function
    )
  )
}

simulate_proxy_data <- function(
  benchmark,
  n,
  rho,
  seed,
  noise_ratio = 0.10
) {
  specification <- proxy_registry()[[benchmark]]
  active <- specification$active
  p <- specification$p
  proxies <- max(active) + seq_along(active)
  set.seed(seed)
  latent <- matrix(stats::rnorm(n * p), nrow = n)
  for (index in seq_along(active)) {
    latent[, proxies[index]] <- rho * latent[, active[index]] +
      sqrt(1 - rho^2) * stats::rnorm(n)
  }
  X <- stats::pnorm(latent)
  colnames(X) <- paste0("X", seq_len(p))
  signal <- specification$fun(X[, active, drop = FALSE])
  noise_sd <- noise_ratio * stats::sd(signal)
  y <- signal + stats::rnorm(n, sd = noise_sd)
  list(
    X = X,
    y = y,
    active = active,
    proxies = proxies,
    p = p,
    noise_sd = noise_sd
  )
}

reference_pme <- function(benchmark, p) {
  references <- readRDS("results/provenance/pme_reference_objects.rds")
  key <- paste0(benchmark, "_rho0.0")
  reference <- references[[key]]
  if (is.null(reference)) stop("Missing PME reference for ", benchmark)
  effects <- numeric(p)
  effects[reference$active] <- reference$pme
  effects
}

pme_metrics <- function(estimate, truth, active, proxies) {
  inactive <- setdiff(seq_along(estimate), active)
  data.frame(
    half_l1 = 0.5 * sum(abs(estimate - truth)),
    active_mae = mean(abs(estimate[active] - truth[active])),
    inactive_mass = sum(estimate[inactive]),
    proxy_mass = sum(estimate[proxies]),
    stringsAsFactors = FALSE
  )
}

run_pme_method <- function(method, data, seed) {
  set.seed(seed)
  start <- proc.time()[["elapsed"]]
  error <- ""
  result <- tryCatch(
    switch(
      method,
      "PRIME" = {
        fit <- prime(data$X, data$y)
        if (!length(fit$selected)) {
          stop("empty_screen: allocation unavailable", call. = FALSE)
        }
        list(
          effects = fit$effects,
          selected = fit$selected,
          value_calls = if (is.null(fit$subset_values)) 0L else nrow(fit$subset_values)
        )
      },
      "Full local-linear PME" = estimate_pme_with_value(
        make_llr_value(data$X, data$y),
        seq_len(data$p),
        data$p
      ),
      "PME-kNN" = estimate_pme_with_value(
        make_local_constant_value(data$X, data$y, neighbors = 2L),
        seq_len(data$p),
        data$p
      ),
      "NCMD-PME" = estimate_pme_with_value(
        make_ncmd_value(data$X, data$y, neighbors = 5L),
        seq_len(data$p),
        data$p
      ),
      "Oracle local-linear PME" = estimate_pme_with_value(
        make_llr_value(data$X, data$y),
        data$active,
        data$p
      ),
      stop("Unknown method: ", method)
    ),
    error = function(condition) {
      error <<- conditionMessage(condition)
      list(
        effects = rep(NA_real_, data$p),
        selected = integer(),
        value_calls = NA_integer_
      )
    }
  )
  result$elapsed_sec <- proc.time()[["elapsed"]] - start
  result$error <- error
  result
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save checkpoint: ", path, call. = FALSE)
  }
}
