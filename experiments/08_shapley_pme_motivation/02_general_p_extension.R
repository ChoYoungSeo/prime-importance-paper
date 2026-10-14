source("software/prime/R/constants.R")
source("software/prime/R/pme.R")

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

total_game_state <- function(a, b, ell, p, beta, rho) {
  weak_count <- p - 2L
  retained_weak <- weak_count - ell
  retained_x1 <- 1L - a
  retained_x2 <- 1L - b
  total_variance <-
    (1 + weak_count * beta * rho)^2 +
    1 + weak_count * beta^2 * (1 - rho^2)

  residual_variance <- (1 - retained_x2) +
    ell * beta^2 * (1 - rho^2)
  if (!retained_x1) {
    conditional_z0_variance <- if (retained_weak == 0L) {
      1
    } else {
      (1 - rho^2) / (1 + (retained_weak - 1L) * rho^2)
    }
    residual_variance <- residual_variance +
      (1 + ell * beta * rho)^2 * conditional_z0_variance
  }
  residual_variance / total_variance
}

general_p_pme <- function(p, beta, rho) {
  stopifnot(p >= 3L, beta > 0, abs(rho) < 1)
  weak_count <- p - 2L
  log_potential <- array(
    NA_real_,
    dim = c(2L, 2L, weak_count + 1L),
    dimnames = list(a = 0:1, b = 0:1, ell = 0:weak_count)
  )
  log_potential[1L, 1L, 1L] <- 0

  for (order in seq_len(p)) {
    for (a in 0:1) {
      for (b in 0:1) {
        ell <- order - a - b
        if (ell < 0L || ell > weak_count) next
        predecessors <- numeric()
        if (a == 1L) {
          predecessors <- c(
            predecessors,
            -log_potential[1L, b + 1L, ell + 1L]
          )
        }
        if (b == 1L) {
          predecessors <- c(
            predecessors,
            -log_potential[a + 1L, 1L, ell + 1L]
          )
        }
        if (ell > 0L) {
          predecessors <- c(
            predecessors,
            log(ell) - log_potential[a + 1L, b + 1L, ell]
          )
        }
        game_value <- total_game_state(a, b, ell, p, beta, rho)
        log_potential[a + 1L, b + 1L, ell + 1L] <-
          log(game_value) - log_sum_exp(predecessors)
      }
    }
  }

  log_full <- log_potential[2L, 2L, weak_count + 1L]
  effects <- c(
    x1 = exp(log_full - log_potential[1L, 2L, weak_count + 1L]),
    x2 = exp(log_full - log_potential[2L, 1L, weak_count + 1L]),
    weak_each = exp(log_full - log_potential[2L, 2L, weak_count])
  )
  c(effects, weak_total = weak_count * effects[["weak_each"]])
}

general_p_shapley <- function(p, beta, rho) {
  weak_count <- p - 2L
  x1 <- 0
  x2 <- 0
  weak <- 0

  for (b in 0:1) {
    for (ell in 0:weak_count) {
      coalition_size <- b + ell
      multiplicity_weight <- exp(
        lchoose(weak_count, ell) - log(p) -
          lchoose(p - 1L, coalition_size)
      )
      x1 <- x1 + multiplicity_weight * (
        total_game_state(1L, b, ell, p, beta, rho) -
          total_game_state(0L, b, ell, p, beta, rho)
      )
    }
  }

  for (a in 0:1) {
    for (ell in 0:weak_count) {
      coalition_size <- a + ell
      multiplicity_weight <- exp(
        lchoose(weak_count, ell) - log(p) -
          lchoose(p - 1L, coalition_size)
      )
      x2 <- x2 + multiplicity_weight * (
        total_game_state(a, 1L, ell, p, beta, rho) -
          total_game_state(a, 0L, ell, p, beta, rho)
      )
    }
  }

  for (a in 0:1) {
    for (b in 0:1) {
      for (ell in 0:(weak_count - 1L)) {
        coalition_size <- a + b + ell
        multiplicity_weight <- exp(
          lchoose(weak_count - 1L, ell) - log(p) -
            lchoose(p - 1L, coalition_size)
        )
        weak <- weak + multiplicity_weight * (
          total_game_state(a, b, ell + 1L, p, beta, rho) -
            total_game_state(a, b, ell, p, beta, rho)
        )
      }
    }
  }

  c(
    x1 = x1,
    x2 = x2,
    weak_each = weak,
    weak_total = weak_count * weak
  )
}

exact_shapley <- function(game) {
  p <- as.integer(round(log2(length(game))))
  effects <- numeric(p)
  for (coordinate in seq_len(p)) {
    bit <- bitwShiftL(1L, coordinate - 1L)
    for (mask in 0:(2^p - 1L)) {
      if (bitwAnd(mask, bit) != 0L) next
      size <- sum(as.logical(intToBits(mask)[seq_len(p)]))
      weight <- 1 / (p * choose(p - 1L, size))
      effects[coordinate] <- effects[coordinate] + weight *
        (game[bitwOr(mask, bit) + 1L] - game[mask + 1L])
    }
  }
  effects
}

exact_game <- function(p, beta, rho) {
  vapply(0:(2^p - 1L), function(mask) {
    deleted <- as.logical(intToBits(mask)[seq_len(p)])
    total_game_state(
      a = as.integer(deleted[1L]),
      b = as.integer(deleted[2L]),
      ell = sum(deleted[-c(1L, 2L)]),
      p = p,
      beta = beta,
      rho = rho
    )
  }, numeric(1))
}

beta <- 0.001
rho <- 0.9

for (p in 3:9) {
  game <- exact_game(p, beta, rho)
  state_pme <- general_p_pme(p, beta, rho)
  state_shapley <- general_p_shapley(p, beta, rho)
  exact_pme <- as.numeric(prime_allocate(game))
  exact_shapley_value <- exact_shapley(game)
  stopifnot(
    max(abs(exact_pme[1:2] - state_pme[1:2])) < 1e-10,
    max(abs(exact_pme[-c(1L, 2L)] - state_pme[["weak_each"]])) < 1e-10,
    max(abs(exact_shapley_value[1:2] - state_shapley[1:2])) < 1e-10,
    max(abs(
      exact_shapley_value[-c(1L, 2L)] - state_shapley[["weak_each"]]
    )) < 1e-10
  )
}

p_grid <- c(3L, 5L, 10L, 20L, 50L, 100L, 500L, 1000L, 5000L)
summary <- do.call(rbind, lapply(p_grid, function(p) {
  shapley <- general_p_shapley(p, beta, rho)
  pme <- general_p_pme(p, beta, rho)
  data.frame(
    p = p,
    beta = beta,
    rho = rho,
    shapley_x1 = shapley[["x1"]],
    shapley_x2 = shapley[["x2"]],
    shapley_weak_each = shapley[["weak_each"]],
    shapley_weak_total = shapley[["weak_total"]],
    pme_x1 = pme[["x1"]],
    pme_x2 = pme[["x2"]],
    pme_weak_each = pme[["weak_each"]],
    pme_weak_total = pme[["weak_total"]]
  )
}))

stopifnot(
  max(abs(rowSums(summary[c(
    "shapley_x1", "shapley_x2", "shapley_weak_total"
  )]) - 1)) < 1e-10,
  max(abs(rowSums(summary[c(
    "pme_x1", "pme_x2", "pme_weak_total"
  )]) - 1)) < 1e-10
)

write.csv(
  summary,
  "results/summary/shapley_pme_proxy_general_p.csv",
  row.names = FALSE
)

print(summary, digits = 7)
cat("Validated the exchangeable general-p calculation against exact enumeration.\n")
