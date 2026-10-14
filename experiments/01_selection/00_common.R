source("experiments/00_reference/benchmarks.R")

selection_profile <- function(profile) {
  switch(
    profile,
    smoke = list(
      benchmarks = c("branin", "ackley3", "friedman"),
      n = 100L,
      designs = c("independent", "balanced_proxy"),
      replications = 1L
    ),
    pilot = list(
      benchmarks = c("branin", "ackley3", "friedman"),
      n = c(100L, 500L),
      designs = c("independent", "spaced_ar1", "balanced_proxy"),
      replications = 5L
    ),
    main = list(
      benchmarks = c("branin", "ackley3", "friedman"),
      n = c(100L, 200L, 500L, 1000L, 2000L),
      designs = c("independent", "spaced_ar1", "balanced_proxy"),
      replications = 50L
    ),
    stop("Unknown selection profile: ", profile, call. = FALSE)
  )
}

selection_registry <- function() {
  list(
    branin = list(
      active_dimension = 2L,
      noise_sd = 2,
      fun = branin_function
    ),
    ackley3 = list(
      active_dimension = 3L,
      noise_ratio = 0.10,
      fun = ackley3_function
    ),
    friedman = list(
      active_dimension = 5L,
      noise_sd = 1,
      fun = friedman_function
    )
  )
}

selection_active_positions <- function(p, dimension, design) {
  if (design %in% c("independent", "balanced_proxy")) {
    return(seq_len(dimension))
  }
  if (design == "spaced_ar1") {
    positions <- as.integer(round((seq_len(dimension) - 0.5) * p / dimension))
    if (anyDuplicated(positions)) {
      stop("Spaced active positions are not distinct.", call. = FALSE)
    }
    return(positions)
  }
  stop("Unknown selection design: ", design, call. = FALSE)
}

selection_covariance <- function(p, active, design, rho = 0.5) {
  if (design == "independent") return(diag(p))
  if (design == "spaced_ar1") {
    index <- seq_len(p)
    return(rho^abs(outer(index, index, "-")))
  }
  if (design == "balanced_proxy") {
    proxies <- max(active) + seq_along(active)
    if (max(proxies) > p) {
      stop("The ambient dimension is too small for the proxy design.",
           call. = FALSE)
    }
    covariance <- diag(p)
    for (position in seq_along(active)) {
      covariance[active[position], proxies[position]] <- rho
      covariance[proxies[position], active[position]] <- rho
    }
    return(covariance)
  }
  stop("Unknown selection design: ", design, call. = FALSE)
}

simulate_selection_data <- function(
  benchmark,
  n,
  p,
  design,
  rho,
  seed
) {
  specification <- selection_registry()[[benchmark]]
  if (is.null(specification)) {
    stop("Unknown selection benchmark: ", benchmark, call. = FALSE)
  }
  active <- selection_active_positions(
    p,
    specification$active_dimension,
    design
  )
  covariance <- selection_covariance(p, active, design, rho)
  set.seed(seed)
  latent <- MASS::mvrnorm(n, mu = numeric(p), Sigma = covariance)
  X <- stats::pnorm(latent)
  colnames(X) <- paste0("X", seq_len(p))
  signal <- specification$fun(X[, active, drop = FALSE])
  noise_sd <- if (!is.null(specification$noise_sd)) {
    specification$noise_sd
  } else {
    specification$noise_ratio * stats::sd(signal)
  }
  y <- signal + stats::rnorm(n, sd = noise_sd)
  proxies <- if (design == "balanced_proxy") {
    max(active) + seq_along(active)
  } else {
    integer()
  }
  list(X = X, y = y, active = active, proxies = proxies)
}

selection_metrics <- function(selected, active, p) {
  selected <- sort(unique(as.integer(selected)))
  active <- sort(unique(as.integer(active)))
  true_positive <- length(intersect(selected, active))
  false_positive <- length(setdiff(selected, active))
  false_negative <- length(setdiff(active, selected))
  data.frame(
    selected = paste(selected, collapse = ","),
    selected_size = length(selected),
    active_selected = true_positive,
    tpr = true_positive / length(active),
    fpr = false_positive / max(1L, p - length(active)),
    false_positives = false_positive,
    false_negatives = false_negative,
    inclusion = false_negative == 0L,
    exact = identical(selected, active),
    stringsAsFactors = FALSE
  )
}

lasso_selector <- function(X, y) {
  fit <- glmnet::cv.glmnet(as.matrix(X), as.numeric(y), alpha = 1)
  coefficient <- as.numeric(stats::coef(fit, s = "lambda.1se"))[-1L]
  which(coefficient != 0)
}

foci_selector <- function(X, y) {
  X <- as.matrix(X)
  fit <- FOCI::foci(
    Y = as.numeric(y),
    X = X,
    num_features = NULL,
    stop = TRUE,
    numCores = 1,
    printIntermed = FALSE
  )
  if (is.null(fit$selectedVar$index)) integer() else
    sort(unique(as.integer(fit$selectedVar$index)))
}

ford_selector <- function(X, y) {
  X <- as.matrix(X)
  fit <- FORD::ford(
    Y = as.numeric(y),
    X = X,
    stop = TRUE,
    standardize = "scale",
    numCores = 1,
    printIntermed = FALSE
  )
  if (is.null(fit$selectedVar$index)) integer() else
    sort(unique(as.integer(fit$selectedVar$index)))
}

run_selection_method <- function(method, X, y, seed) {
  set.seed(seed)
  error <- ""
  start <- proc.time()[["elapsed"]]
  selected <- tryCatch(
    switch(
      method,
      "PRIME screen" = prime::prime_screen(X, y)$selected,
      "FIRST" = {
        importance <- first::first(
          X = as.data.frame(X),
          y = as.numeric(y),
          importance = "total",
          verbose = FALSE
        )
        which(is.finite(importance) & abs(importance) > 0)
      },
      "LASSO" = lasso_selector(X, y),
      "FOCI" = foci_selector(X, y),
      "FORD" = ford_selector(X, y),
      stop("Unknown selection method: ", method, call. = FALSE)
    ),
    error = function(condition) {
      error <<- conditionMessage(condition)
      integer()
    }
  )
  list(
    selected = sort(unique(as.integer(selected))),
    elapsed_sec = proc.time()[["elapsed"]] - start,
    error = error
  )
}

selection_atomic_save <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save checkpoint: ", path, call. = FALSE)
  }
}

selection_output_path <- function(profile) {
  if (profile == "main") {
    "results/raw/selection_v052_runs.csv"
  } else {
    file.path("results/raw", paste0("selection_v052_", profile, ".csv"))
  }
}
