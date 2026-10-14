library(prime)
library(firstllr)

validate_case <- function(seed, n, p, correlated = FALSE) {
  set.seed(seed)
  X <- matrix(runif(n * p, -1, 1), nrow = n)
  if (correlated && p >= 5L) {
    X[, 3] <- 0.9 * X[, 2] + 0.1 * rnorm(n)
    X[, 4] <- 0.9 * X[, 1] + 0.1 * rnorm(n)
    y <- X[, 1] + X[, 2] - X[, 3] + X[, 5]^2 +
      0.1 * rnorm(n)
  } else {
    y <- 2 * X[, 1] - X[, 2] + X[, 3]^2 +
      sin(pi * X[, 4] * X[, 5]) + 0.2 * rnorm(n)
  }

  reference <- firstllr::first_llr(X, y)
  candidate <- prime::prime_screen(X, y)
  stopifnot(
    identical(candidate$candidates, reference$candidates),
    identical(candidate$selected, reference$selected),
    isTRUE(all.equal(
      candidate$forward$trace$candidate_value,
      reference$forward$trace$candidate_value,
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      candidate$backward$trace$statistic + 1,
      reference$backward$trace$statistic,
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      candidate$backward$trace$cutoff + 1,
      reference$backward$trace$f_critical,
      tolerance = 1e-12
    )),
    identical(
      candidate$backward$trace$removed,
      reference$backward$trace$removed
    )
  )
  data.frame(
    seed = seed,
    n = n,
    p = p,
    correlated = correlated,
    candidates = paste(candidate$candidates, collapse = ","),
    selected = paste(candidate$selected, collapse = ","),
    stringsAsFactors = FALSE
  )
}

grid <- expand.grid(
  seed = c(3L, 5L, 11L),
  n = c(100L, 500L),
  p = c(10L, 50L),
  correlated = c(FALSE, TRUE),
  KEEP.OUT.ATTRS = FALSE
)
results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(row) {
  do.call(validate_case, as.list(grid[row, ]))
}))
write.csv(
  results,
  "results/provenance/selection_implementation_parity.csv",
  row.names = FALSE
)
cat("Selection parity passed for", nrow(results), "cases.\n")

