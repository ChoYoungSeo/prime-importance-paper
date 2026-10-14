ar1_uniform_design <- function(n, p, rho = 0) {
  if (rho == 0) return(matrix(stats::runif(n * p), nrow = n, ncol = p))
  index <- seq_len(p)
  sigma <- rho^abs(outer(index, index, "-"))
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p) %*% chol(sigma)
  stats::pnorm(z)
}

scale_range <- function(u, lower, upper) lower + (upper - lower) * u

branin_function <- function(U) {
  x1 <- scale_range(U[, 1], -5, 10)
  x2 <- scale_range(U[, 2], 0, 15)
  (x2 - 5.1 * x1^2 / (4 * pi^2) + 5 * x1 / pi - 6)^2 +
    10 * (1 - 1 / (8 * pi)) * cos(x1) + 10
}

friedman_function <- function(U) {
  10 * sin(pi * U[, 1] * U[, 2]) +
    20 * (U[, 3] - 0.5)^2 + 10 * U[, 4] + 5 * U[, 5]
}

moon_function <- function(U) {
  0.5 * (U[, 1] + U[, 2])^2 + U[, 3]^2
}

ishigami_function <- function(U, a = 7, b = 0.1) {
  x <- -pi + 2 * pi * U[, 1:3, drop = FALSE]
  sin(x[, 1]) + a * sin(x[, 2])^2 + b * x[, 3]^4 * sin(x[, 1])
}

otl_function <- function(U) {
  rb1 <- scale_range(U[, 1], 50, 150)
  rb2 <- scale_range(U[, 2], 25, 70)
  rf <- scale_range(U[, 3], 0.5, 3)
  rc1 <- scale_range(U[, 4], 1.2, 2.5)
  rc2 <- scale_range(U[, 5], 0.25, 1.2)
  beta <- scale_range(U[, 6], 50, 300)
  vb1 <- 12 * rb2 / (rb1 + rb2)
  denominator <- beta * (rc2 + 9) + rf
  (vb1 + 0.74) * beta * (rc2 + 9) / denominator +
    11.35 * rf / denominator +
    0.74 * rf * beta * (rc2 + 9) / (denominator * rc1)
}

piston_function <- function(U) {
  m <- scale_range(U[, 1], 30, 60)
  s <- scale_range(U[, 2], 0.005, 0.020)
  v0 <- scale_range(U[, 3], 0.002, 0.010)
  k <- scale_range(U[, 4], 1000, 5000)
  p0 <- scale_range(U[, 5], 90000, 110000)
  ta <- scale_range(U[, 6], 290, 296)
  t0 <- scale_range(U[, 7], 340, 360)
  a <- p0 * s + 19.62 * m - k * v0 / s
  v <- s * (sqrt(a^2 + 4 * k * p0 * v0 * ta / t0) - a) / (2 * k)
  2 * pi * sqrt(m / (k + s^2 * p0 * v0 * ta / (t0 * v^2)))
}

borehole_function <- function(U) {
  rw <- scale_range(U[, 1], 0.05, 0.15)
  r <- scale_range(U[, 2], 100, 50000)
  tu <- scale_range(U[, 3], 63070, 115600)
  hu <- scale_range(U[, 4], 990, 1110)
  tl <- scale_range(U[, 5], 63.1, 116)
  hl <- scale_range(U[, 6], 700, 820)
  l <- scale_range(U[, 7], 1120, 1680)
  kw <- scale_range(U[, 8], 9855, 12045)
  log_ratio <- log(r / rw)
  2 * pi * tu * (hu - hl) /
    (log_ratio * (1 + 2 * l * tu / (log_ratio * rw^2 * kw) + tu / tl))
}

bratley_function <- function(U) {
  if (ncol(U) < 5L) stop("Bratley requires five columns.")
  rowSums(vapply(seq_len(5L), function(i) {
    (-1)^i * apply(U[, seq_len(i), drop = FALSE], 1L, prod)
  }, numeric(nrow(U))))
}

sobol_g_function <- function(U, a = c(0, 0.5, 3, 9, 20, 50)) {
  if (ncol(U) < length(a)) {
    stop("Sobol g requires ", length(a), " columns.")
  }
  factors <- sweep(
    abs(4 * U[, seq_along(a), drop = FALSE] - 2),
    2L,
    a,
    "+"
  )
  factors <- sweep(factors, 2L, 1 + a, "/")
  apply(factors, 1L, prod)
}

sobol_g3_function <- function(U) {
  sobol_g_function(U, a = c(0, 0.5, 1))
}

ackley3_function <- function(U) {
  if (ncol(U) < 3L) stop("Ackley-3 requires three columns.")
  X <- -5 + 10 * U[, 1:3, drop = FALSE]
  -20 * exp(-0.2 * sqrt(rowMeans(X^2))) -
    exp(rowMeans(cos(2 * pi * X))) + 20 + exp(1)
}

hartmann3_function <- function(U) {
  if (ncol(U) < 3L) stop("Hartmann-3 requires three columns.")
  alpha <- c(1, 1.2, 3, 3.2)
  A <- rbind(
    c(3, 10, 30),
    c(0.1, 10, 35),
    c(3, 10, 30),
    c(0.1, 10, 35)
  )
  P <- 1e-4 * rbind(
    c(3689, 1170, 2673),
    c(4699, 4387, 7470),
    c(1091, 8732, 5547),
    c(381, 5743, 8828)
  )
  terms <- vapply(seq_len(4L), function(i) {
    alpha[i] * exp(-rowSums(
      sweep(U[, 1:3, drop = FALSE], 2L, P[i, ], "-")^2 *
        matrix(A[i, ], nrow = nrow(U), ncol = 3L, byrow = TRUE)
    ))
  }, numeric(nrow(U)))
  -rowSums(terms)
}

benchmark_registry <- function() {
  list(
    branin = list(active = 1:2, noise_sd = 2, fun = branin_function),
    moon = list(active = 1:3, noise_sd = 0.5, fun = moon_function),
    friedman = list(active = 1:5, noise_sd = 1, fun = friedman_function),
    ishigami = list(active = 1:3, noise_ratio = 0.10, fun = ishigami_function),
    otl = list(active = 1:6, noise_ratio = 0.10, fun = otl_function),
    piston = list(active = 1:7, noise_ratio = 0.10, fun = piston_function),
    borehole = list(active = 1:8, noise_ratio = 0.10, fun = borehole_function),
    bratley = list(active = 1:5, noise_ratio = 0.10, fun = bratley_function),
    ackley3 = list(active = 1:3, noise_ratio = 0.10, fun = ackley3_function),
    hartmann3 = list(active = 1:3, noise_ratio = 0.10,
                     fun = hartmann3_function),
    sobol_g3 = list(active = 1:3, noise_ratio = 0.10, fun = sobol_g3_function),
    sobol_g = list(active = 1:6, noise_ratio = 0.10, fun = sobol_g_function)
  )
}

simulate_benchmark <- function(name, n, p, rho = 0, seed = 1L) {
  benchmark <- benchmark_registry()[[name]]
  if (is.null(benchmark)) stop("Unknown benchmark: ", name)
  if (p < max(benchmark$active)) stop("p is smaller than the active dimension.")
  set.seed(seed)
  X <- ar1_uniform_design(n, p, rho)
  signal <- benchmark$fun(X[, benchmark$active, drop = FALSE])
  sigma <- if (!is.null(benchmark$noise_sd)) {
    benchmark$noise_sd
  } else {
    benchmark$noise_ratio * stats::sd(signal)
  }
  list(
    X = X,
    y = signal + stats::rnorm(n, sd = sigma),
    signal = signal,
    active = benchmark$active,
    noise_sd = sigma
  )
}

shapley_from_coalition_values <- function(values) {
  q <- log2(length(values))
  if (!isTRUE(all.equal(q, round(q))) || q > 20L) {
    stop("Coalition values must have length 2^q with q <= 20.")
  }
  q <- as.integer(round(q))
  phi <- numeric(q)
  denominator <- factorial(q)
  for (j in seq_len(q)) {
    bit <- bitwShiftL(1L, j - 1L)
    for (mask in 0:(2^q - 1L)) {
      if (bitwAnd(mask, bit) == 0L) {
        size <- sum(as.logical(intToBits(mask)[seq_len(q)]))
        weight <- factorial(size) * factorial(q - size - 1L) / denominator
        phi[j] <- phi[j] + weight *
          (values[bitwOr(mask, bit) + 1L] - values[mask + 1L])
      }
    }
  }
  phi <- pmax(phi, 0)
  if (sum(phi) > 0) phi <- phi / sum(phi)
  phi
}

conditional_gaussian_copula_pair <- function(n, p, rho, S) {
  S <- sort(unique(as.integer(S)))
  complement <- setdiff(seq_len(p), S)
  sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))

  if (length(S) == 0L) {
    z1 <- matrix(stats::rnorm(n * p), nrow = n) %*% chol(sigma)
    z2 <- matrix(stats::rnorm(n * p), nrow = n) %*% chol(sigma)
    return(list(X1 = stats::pnorm(z1), X2 = stats::pnorm(z2)))
  }

  sigma_ss <- sigma[S, S, drop = FALSE]
  z_s <- matrix(stats::rnorm(n * length(S)), nrow = n) %*% chol(sigma_ss)
  if (length(complement) == 0L) {
    z <- matrix(0, nrow = n, ncol = p)
    z[, S] <- z_s
    X <- stats::pnorm(z)
    return(list(X1 = X, X2 = X))
  }

  conditional_mean_map <- solve(
    sigma_ss,
    sigma[S, complement, drop = FALSE]
  )
  conditional_mean <- z_s %*% conditional_mean_map
  conditional_covariance <-
    sigma[complement, complement, drop = FALSE] -
    sigma[complement, S, drop = FALSE] %*% conditional_mean_map
  conditional_covariance <-
    (conditional_covariance + t(conditional_covariance)) / 2
  eigen_decomposition <- eigen(conditional_covariance, symmetric = TRUE)
  eigen_decomposition$values <- pmax(eigen_decomposition$values, 0)
  conditional_factor <-
    eigen_decomposition$vectors %*%
    diag(sqrt(eigen_decomposition$values), nrow = length(complement))

  z1 <- matrix(0, nrow = n, ncol = p)
  z2 <- matrix(0, nrow = n, ncol = p)
  z1[, S] <- z_s
  z2[, S] <- z_s
  z1[, complement] <- conditional_mean +
    matrix(stats::rnorm(n * length(complement)), nrow = n) %*%
    t(conditional_factor)
  z2[, complement] <- conditional_mean +
    matrix(stats::rnorm(n * length(complement)), nrow = n) %*%
    t(conditional_factor)
  list(X1 = stats::pnorm(z1), X2 = stats::pnorm(z2))
}

conditional_pickfreeze_reference <- function(name, rho, n_ref = 100000L,
                                             seed = 1L) {
  benchmark <- benchmark_registry()[[name]]
  if (is.null(benchmark)) stop("Unknown benchmark: ", name)
  active <- benchmark$active
  q <- length(active)
  if (q > 20L) stop("Reference calculation is limited to 20 active variables.")
  if (n_ref < 100L || n_ref %% 2L != 0L) {
    stop("n_ref must be an even integer of at least 100.")
  }

  set.seed(seed)
  denominator_design <- ar1_uniform_design(n_ref, q, rho)
  denominator_signal <- benchmark$fun(denominator_design)
  halves <- rep(1:2, each = n_ref / 2L)
  variance_all <- stats::var(denominator_signal)
  variance_half <- vapply(1:2, function(h) {
    stats::var(denominator_signal[halves == h])
  }, numeric(1))

  number_values <- 2^q
  values <- numeric(number_values)
  values_half <- matrix(0, nrow = number_values, ncol = 2L)
  values[number_values] <- 1
  values_half[number_values, ] <- 1

  if (number_values > 2L) {
    for (mask in 1:(number_values - 2L)) {
      S <- which(as.logical(intToBits(mask)[seq_len(q)]))
      set.seed(seed + 1009L * mask)
      pair <- conditional_gaussian_copula_pair(n_ref, q, rho, S)
      y1 <- benchmark$fun(pair$X1)
      y2 <- benchmark$fun(pair$X2)
      values[mask + 1L] <- stats::cov(y1, y2) / variance_all
      values_half[mask + 1L, ] <- vapply(1:2, function(h) {
        index <- halves == h
        stats::cov(y1[index], y2[index]) / variance_half[h]
      }, numeric(1))
    }
  }

  values <- pmin(1, pmax(0, values))
  values_half[] <- pmin(1, pmax(0, values_half))
  phi <- shapley_from_coalition_values(values)
  phi_half <- apply(values_half, 2L, shapley_from_coalition_values)
  list(
    phi = phi,
    coalition_values = values,
    coalition_values_half = values_half,
    phi_half = phi_half,
    max_half_phi_difference = max(abs(phi_half[, 1L] - phi_half[, 2L])),
    benchmark = name,
    active = active,
    rho = rho,
    n_ref = n_ref,
    seed = seed,
    method = "conditional Gaussian-copula pick-freeze covariance"
  )
}

reference_shapley <- function(name, p, rho, cache_dir) {
  benchmark <- benchmark_registry()[[name]]
  if (is.null(benchmark)) stop("Unknown benchmark: ", name)
  active_dimension <- length(benchmark$active)
  rho_key <- gsub("\\.", "p", sprintf("%.4f", rho))
  pattern <- sprintf("^reference_%s_s%s_rho%s_cpfr_.*\\.rds$",
                     name, active_dimension, rho_key)
  paths <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
  if (length(paths) != 1L) {
    stop("Expected one cached reference for ", name, " at rho=", rho,
         "; found ", length(paths), ".")
  }
  reference <- readRDS(paths)
  if (!is.list(reference) || is.null(reference$phi)) {
    stop("Reference cache does not contain conditional pick-freeze metadata: ", paths)
  }
  active_phi <- as.numeric(reference$phi)
  out <- numeric(p)
  out[benchmark$active] <- active_phi
  out
}
