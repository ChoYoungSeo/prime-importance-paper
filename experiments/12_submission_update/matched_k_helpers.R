
shapley_exact <- function(game, p) {
  masks <- 0:(2^p - 1L)
  vapply(seq_len(p), function(j) {
    bit <- bitwShiftL(1L, j - 1L)
    without <- masks[bitwAnd(masks, bit) == 0L]
    sizes <- vapply(without, function(mask) {
      sum(vapply(seq_len(p), function(k) {
        bitwAnd(mask, bitwShiftL(1L, k - 1L)) != 0L
      }, logical(1)))
    }, integer(1))
    weights <- 1 / (p * choose(p - 1L, sizes))
    sum(weights * (game[without + bit + 1L] - game[without + 1L]))
  }, numeric(1))
}

matched_k_local_constant_game <- function(X, y) {
  X <- standardize_predictors(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  response_variance <- stats::var(y)
  cache <- new.env(parent = emptyenv())

  risk <- function(subset) {
    subset <- sort(unique(as.integer(subset)))
    if (!length(subset)) return(response_variance)
    key <- paste(subset, collapse = ",")
    if (!exists(key, cache, inherits = FALSE)) {
      k <- prime:::.prime_neighborhood_size(n, length(subset))
      stopifnot(length(k) == 1L, is.finite(k), k > 1L)
      projected <- X[, subset, drop = FALSE]
      neighbors <- prime:::.prime_query_neighbors(projected, seq_len(n), k)
      stopifnot(
        nrow(neighbors) == n, ncol(neighbors) == k,
        !any(neighbors == matrix(seq_len(n), nrow = n, ncol = k))
      )
      local_variance <- apply(
        matrix(y[neighbors], nrow = n, ncol = k),
        1L,
        stats::var
      )
      value <- mean(local_variance)
      stopifnot(is.finite(value), value >= 0)
      assign(key, value, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }

  full <- seq_len(p)
  noise_floor <- risk(full)
  signal <- max(response_variance - noise_floor, 0)
  guard <- sqrt(.Machine$double.eps) * response_variance
  if (!is.finite(signal) || signal <= guard) {
    stop("Matched-K local-constant allocation is unavailable.", call. = FALSE)
  }

  full_mask <- 2^p - 1L
  raw <- numeric(full_mask + 1L)
  complement_risk <- numeric(full_mask + 1L)
  for (mask in 0:full_mask) {
    keep <- which(!vapply(seq_len(p), function(position) {
      bitwAnd(mask, bitwShiftL(1L, position - 1L)) != 0L
    }, logical(1)))
    risk_keep <- risk(keep)
    complement_risk[mask + 1L] <- risk_keep
    raw[mask + 1L] <- max(risk_keep - noise_floor, 0) / signal
  }
  raw[1L] <- 0
  raw[length(raw)] <- 1

  pme <- prime::prime_allocate(raw)
  projected <- as.numeric(attr(pme, "projected_values"))
  shapley <- shapley_exact(projected, p)
  stopifnot(
    length(raw) == 2^p, length(projected) == 2^p,
    all(is.finite(raw)), all(is.finite(projected)),
    abs(sum(pme) - 1) < 1e-10,
    abs(sum(shapley) - 1) < 1e-10
  )
  list(
    pme = as.numeric(pme),
    shapley = shapley,
    raw_game = raw,
    projected_game = projected,
    complement_risk = complement_risk,
    noise_floor = noise_floor,
    signal = signal,
    response_variance = response_variance,
    value_calls = length(ls(cache, all.names = TRUE)),
    k_by_dimension = vapply(seq_len(p), function(d) {
      prime:::.prime_neighborhood_size(n, d)
    }, integer(1))
  )
