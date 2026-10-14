source("experiments/03_pme/00_common.R")
source("experiments/09_observed_data/00_common.R")

concrete_representations <- function() {
  data <- load_observed_data("concrete")
  raw <- data$X
  engineering <- raw
  binder <- rowSums(raw[, c("cement", "slag", "ash"), drop = FALSE])
  if (any(!is.finite(binder) | binder <= 0)) {
    stop("Concrete binder mass must be positive.", call. = FALSE)
  }
  engineering[, "water"] <- raw[, "water"] / binder
  colnames(engineering)[colnames(engineering) == "water"] <-
    "water_binder_ratio"
  list(y = data$y, raw = raw, engineering = engineering)
}

concrete_method <- function(method, X, y) {
  switch(
    method,
    "PRIME" = {
      fit <- prime::prime(X, y)
      if (!length(fit$selected)) {
        stop("empty_screen: allocation unavailable", call. = FALSE)
      }
      list(effects = fit$effects, selected = fit$selected,
           value_calls = if (is.null(fit$subset_values)) 0L else nrow(fit$subset_values))
    },
    "Full LL-PME" = estimate_pme_with_value(
      make_llr_value(X, y), seq_len(ncol(X)), ncol(X)
    ),
    "PME-kNN" = estimate_pme_with_value(
      make_local_constant_value(X, y, neighbors = 2L),
      seq_len(ncol(X)), ncol(X)
    ),
    stop("Unknown concrete method: ", method, call. = FALSE)
  )
}

concrete_standardize <- function(train, test) {
  center <- colMeans(train)
  scale <- apply(train, 2L, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(
    train = sweep(sweep(train, 2L, center, "-"), 2L, scale, "/"),
    test = sweep(sweep(test, 2L, center, "-"), 2L, scale, "/")
  )
}

concrete_knn_predict <- function(train_X, train_y, test_X, k) {
  FNN::knn.reg(
    train = as.matrix(train_X),
    test = as.matrix(test_X),
    y = as.numeric(train_y),
    k = as.integer(k)
  )$pred
}

concrete_tune_k <- function(X, y, variables, fold, k_grid) {
  errors <- vapply(k_grid, function(k) {
    fold_error <- vapply(sort(unique(fold)), function(level) {
      validation <- which(fold == level)
      training <- which(fold != level)
      scaled <- concrete_standardize(
        X[training, variables, drop = FALSE],
        X[validation, variables, drop = FALSE]
      )
      prediction <- concrete_knn_predict(
        scaled$train, y[training], scaled$test, k
      )
      mean((y[validation] - prediction)^2)
    }, numeric(1))
    sqrt(mean(fold_error))
  }, numeric(1))
  as.integer(k_grid[which.min(errors)])
}

concrete_atomic_save <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save concrete checkpoint: ", path, call. = FALSE)
  }
}
