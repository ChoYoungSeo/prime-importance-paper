library(prime)

set.seed(5)
n <- 500L
p <- 10L
X <- matrix(0, nrow = n, ncol = p)
X[, 1] <- runif(n, -1, 1)
X[, 2] <- runif(n, -1, 1)
X[, 3] <- 0.9 * X[, 2] + 0.1 * rnorm(n)
X[, 4] <- 0.9 * X[, 1] + 0.1 * rnorm(n)
X[, 5:p] <- matrix(runif(n * (p - 4L), -1, 1), nrow = n)
y <- X[, 1] + X[, 2] - X[, 3] + X[, 5]^2 +
  0.1 * rnorm(n)

fit <- prime(X, y)
stopifnot(
  identical(fit$screen$forward$ridge_initial, c(1L, 2L, 3L)),
  identical(fit$selected, c(1L, 2L, 3L, 5L)),
  abs(sum(fit$effects) - 1) < 1e-8,
  all(fit$effects[-fit$selected] == 0)
)

rho <- 0.9
proxy_game <- c(
  0,
  (1 - rho^2) / 2,
  1 / 2,
  1 - rho^2 / 2,
  0,
  1 / 2,
  1 / 2,
  1
)
proxy_pme <- prime_allocate(proxy_game)
stopifnot(max(abs(proxy_pme - c(0.5, 0.5, 0))) < 1e-10)

output <- data.frame(
  variable = names(fit$effects),
  selected = seq_along(fit$effects) %in% fit$selected,
  pme = as.numeric(fit$effects)
)
dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)
write.csv(
  output,
  "results/summary/smoke_correlated_suppressor.csv",
  row.names = FALSE
)

cat("Smoke validation passed.\n")
cat("Selected:", paste(fit$selected, collapse = ","), "\n")
cat("PME:", paste(round(fit$effects, 4), collapse = ","), "\n")
