source("experiments/03_pme/00_common.R")

selection_registry <- function() {
  list(
    branin = list(
      active = 1:2,
      noise_sd = 2,
      fun = branin_function
    ),
    ackley3 = list(
      active = 1:3,
      noise_ratio = 0.10,
      fun = ackley3_function
    ),
    friedman = list(
      active = 1:5,
      noise_sd = 1,
      fun = friedman_function
    )
  )
}

selection_design_covariance <- function(p, active, design, rho) {
  if (design == "independent") return(diag(p))
  if (design == "balanced_proxy") {
    proxies <- max(active) + seq_along(active)
    covariance <- diag(p)
    for (position in seq_along(active)) {
      covariance[active[position], proxies[position]] <- rho
      covariance[proxies[position], active[position]] <- rho
    }
    return(covariance)
  }
  stop("Unsupported selector-to-PME design: ", design, call. = FALSE)
}

simulate_selection_task <- function(benchmark, n, p, design, rho, seed) {
  specification <- selection_registry()[[benchmark]]
  if (is.null(specification)) stop("Unknown benchmark: ", benchmark)
  active <- specification$active
  covariance <- selection_design_covariance(p, active, design, rho)
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
  list(X = X, y = y, active = active, proxies = proxies, p = p)
}

parse_selected <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(integer())
  sort(unique(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])))
}

selection_pme_reference <- function(benchmark, p) {
  references <- readRDS("results/provenance/pme_reference_objects.rds")
  reference <- references[[paste0(benchmark, "_rho0.0")]]
  if (is.null(reference)) stop("Missing PME truth for ", benchmark)
  active_effects <- reference$pme
  active <- reference$active
  effects <- numeric(p)
  effects[active] <- active_effects
  effects
}

selected_set_pme <- function(X, y, selected, p) {
  selected <- sort(unique(as.integer(selected)))
  if (!length(selected)) {
    return(list(
      effects = rep(NA_real_, p),
      selected = integer(),
      value_calls = NA_integer_,
      feasible = FALSE,
      reason = "empty selected set"
    ))
  }
  if (length(selected) > 12L) {
    return(list(
      effects = rep(NA_real_, p),
      selected = selected,
      value_calls = NA_integer_,
      feasible = FALSE,
      reason = "selected dimension exceeds the experiment budget"
    ))
  }
  result <- estimate_pme_with_value(
    make_llr_value(X, y), selected, p
  )
  result$feasible <- TRUE
  result$reason <- ""
  result
}
