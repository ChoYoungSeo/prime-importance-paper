library(prime)
source("experiments/03_pme/00_common.R")

direct_knn_value <- function(X, y, subset, neighbors = 2L) {
  projected <- X[, subset, drop = FALSE]
  distance <- as.matrix(stats::dist(projected))
  index <- t(apply(distance, 2L, order)[seq_len(neighbors), , drop = FALSE])
  local_variance <- apply(
    matrix(y[index], ncol = neighbors),
    1L,
    stats::var
  )
  1 - mean(local_variance) / stats::var(y)
}

set.seed(20260816)
X <- standardize_predictors(matrix(stats::runif(1200L), ncol = 3L))
y <- X[, 1] + 0.7 * X[, 2]^2 + stats::rnorm(nrow(X), sd = 0.05)
subsets <- list(1L, 2L, 3L, c(1L, 2L), c(1L, 3L), c(2L, 3L), 1:3)
value <- make_local_constant_value(X, y, neighbors = 2L)
rows <- do.call(rbind, lapply(subsets, function(subset) {
  direct <- direct_knn_value(X, y, subset, neighbors = 2L)
  cached <- value(subset)
  reference <- getFromNamespace(
    "estim.VE.knn", "sensitivity"
  )(
    dataX = X,
    y = y,
    subset = subset,
    n = nrow(X),
    n.knn = 2L,
    n.limit = 2000L
  )
  data.frame(
    subset = paste(subset, collapse = ","),
    direct = direct,
    cached = cached,
    sensitivity_reference = reference,
    direct_difference = abs(direct - cached),
    reference_difference = abs(reference - cached),
    stringsAsFactors = FALSE
  )
}))
stopifnot(
  max(rows$direct_difference) < 1e-12,
  max(rows$reference_difference) < 1e-12
)
write.csv(
  rows,
  "results/provenance/pme_knn_subset_value_parity.csv",
  row.names = FALSE
)
cat("PME reduction and kNN subset-value parity validated.\n")
