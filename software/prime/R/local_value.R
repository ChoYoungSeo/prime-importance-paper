.prime_query_neighbors <- function(X, query_index, k) {
  algorithm <- if (nrow(X) <= 2000L && ncol(X) >= 5L) {
    "brute"
  } else {
    "kd_tree"
  }
  all_queries <- length(query_index) == nrow(X) &&
    identical(query_index, seq_len(nrow(X)))
  if (all_queries) {
    neighbors <- FNN::get.knn(
      X,
      k = as.integer(k),
      algorithm = algorithm
    )$nn.index
    self <- neighbors == matrix(query_index, nrow(neighbors), ncol(neighbors))
    affected <- which(rowSums(self) > 0L)
    if (length(affected)) {
      # Preserve every already-valid neighbor entry. FNN can discard another
      # zero-distance row instead of self when projected predictors tie.
      candidates <- FNN::get.knnx(
        X, X[query_index[affected], , drop = FALSE],
        k = as.integer(k) + 1L, algorithm = algorithm
      )$nn.index
      for (position in seq_along(affected)) {
        row <- affected[position]
        eligible <- setdiff(
          candidates[position, ], c(query_index[row], neighbors[row, ])
        )
        if (!length(eligible)) {
          stop("No eligible replacement for a query's self neighbor.", call. = FALSE)
        }
        neighbors[row, which(self[row, ])] <- eligible[1L]
      }
    }
    if (any(neighbors == matrix(query_index, nrow(neighbors), ncol(neighbors)))) {
      stop("A query observation remained in its own neighborhood.", call. = FALSE)
    }
    return(neighbors)
  }

  candidates <- FNN::get.knnx(
    X,
    X[query_index, , drop = FALSE],
    k = as.integer(k) + 1L,
    algorithm = algorithm
  )$nn.index
  query_matrix <- matrix(
    query_index,
    nrow = length(query_index),
    ncol = ncol(candidates)
  )
  self_match <- candidates == query_matrix
  self_position <- max.col(self_match, ties.method = "first")
  self_position[rowSums(self_match) == 0L] <- ncol(candidates) + 1L

  target_column <- matrix(
    rep(seq_len(k), each = length(query_index)),
    nrow = length(query_index),
    ncol = k
  )
  source_column <- target_column + (target_column >= self_position)
  row_index <- matrix(
    rep(seq_along(query_index), times = k),
    nrow = length(query_index),
    ncol = k
  )
  neighbors <- matrix(
    candidates[cbind(as.vector(row_index), as.vector(source_column))],
    nrow = length(query_index),
    ncol = k
  )
  if (any(neighbors == matrix(
    query_index,
    nrow = length(query_index),
    ncol = k
  ))) {
    stop("A query observation remained in its own neighborhood.",
         call. = FALSE)
  }
  neighbors
}

.prime_variance_details <- function(data, subset, k) {
  subset <- sort(unique(as.integer(subset)))
  if (!length(subset)) {
    return(list(
      risk = data$response_variance,
      mean_sse = NA_real_,
      residual_df = NA_real_,
      contribution = rep(
        data$response_variance, length(data$query_index)
      ),
      query_count = length(data$query_index)
    ))
  }
  key <- paste(k, paste(subset, collapse = ","), sep = ":")
  if (exists(key, envir = data$value_cache, inherits = FALSE)) {
    return(get(key, envir = data$value_cache, inherits = FALSE))
  }
  Xs <- data$X[, subset, drop = FALSE]
  neighbors <- .prime_query_neighbors(
    Xs, data$query_index, k
  ) - 1L
  details <- query_variance_cpp(
    Xs,
    data$y,
    data$query_index - 1L,
    neighbors,
    data$lambda
  )
  details$contribution <- as.numeric(details$contribution)
  details$query_count <- length(data$query_index)
  assign(key, details, envir = data$value_cache)
  details
}

.prime_risk <- function(data, subset, k) {
  .prime_variance_details(data, subset, k)$risk
}

.prime_value <- function(data, subset, k) {
  if (!length(subset)) return(0)
  risk <- .prime_risk(data, subset, k)
  max(0, data$response_variance - risk) / data$response_variance
}

.prime_forward_score <- function(data, subset, k) {
  if (!length(subset)) return(0)
  risk <- .prime_risk(data, subset, k)
  (data$response_variance - risk) / data$response_variance
}

.prime_deletion_diagnostic <- function(
  data, selected, variable, k, full_details = NULL
) {
  selected <- sort(unique(as.integer(selected)))
  if (!variable %in% selected) {
    stop("The tested coordinate is not in the selected set.", call. = FALSE)
  }
  if (is.null(full_details)) {
    full_details <- .prime_variance_details(data, selected, k)
  }
  reduced <- setdiff(selected, variable)
  reduced_details <- .prime_variance_details(data, reduced, k)
  paired_loss <- reduced_details$contribution - full_details$contribution
  raw_loss <- mean(paired_loss)
  query_scale <- if (length(paired_loss) > 1L) {
    stats::sd(paired_loss) / sqrt(length(paired_loss))
  } else {
    0
  }
  threshold <- query_scale * sqrt(
    2 * log(max(2L, length(selected)))
  )
  list(
    raw_loss = raw_loss,
    clipped_loss = max(raw_loss, 0),
    query_scale = query_scale,
    threshold = threshold,
    calibrated_loss = max(raw_loss - threshold, 0),
    full_risk = full_details$risk,
    reduced_risk = reduced_details$risk,
    paired_loss = paired_loss
  )
}

.prime_forward_diagnostic <- function(
  data, selected, variable, k
) {
  selected <- sort(unique(as.integer(selected)))
  current <- .prime_variance_details(data, selected, k)
  augmented <- .prime_variance_details(
    data, sort(c(selected, variable)), k
  )
  paired_gain <- current$contribution - augmented$contribution
  raw_gain <- mean(paired_gain)
  query_scale <- if (length(paired_gain) > 1L) {
    stats::sd(paired_gain) / sqrt(length(paired_gain))
  } else {
    0
  }
  augmented_dimension <- length(selected) + 1L
  threshold <- query_scale * sqrt(
    2 * log(max(2L, augmented_dimension))
  )
  list(
    raw_gain = raw_gain,
    query_scale = query_scale,
    threshold = threshold,
    calibrated_gain = max(raw_gain - threshold, 0),
    paired_gain = paired_gain
  )
}

.prime_forward_values <- function(data, selected, candidates, k) {
  batch <- forward_values_cpp(
    data$X,
    data$y,
    data$query_index - 1L,
    as.integer(selected) - 1L,
    as.integer(candidates) - 1L,
    as.integer(k),
    data$lambda,
    data$response_variance,
    .prime_worker_count(length(candidates))
  )
  values <- as.numeric(batch$values)
  boundary_ties <- which(as.logical(batch$boundary_ties))
  if (length(boundary_ties)) {
    values[boundary_ties] <- vapply(boundary_ties, function(index) {
      .prime_forward_score(data, c(selected, candidates[index]), k)
    }, numeric(1))
  }
  values
}

.prime_variance_details_r <- function(data, subset, k) {
  subset <- sort(unique(as.integer(subset)))
  if (!length(subset)) {
    return(list(
      risk = data$response_variance,
      mean_sse = NA_real_,
      residual_df = NA_real_,
      contribution = rep(
        data$response_variance, length(data$query_index)
      )
    ))
  }
  Xs <- data$X[, subset, drop = FALSE]
  neighbors <- .prime_query_neighbors(Xs, data$query_index, k)
  local_sse <- vapply(seq_along(data$query_index), function(i) {
    index <- neighbors[i, ]
    centered <- sweep(
      Xs[index, , drop = FALSE],
      2L,
      Xs[data$query_index[i], ],
      "-"
    )
    design <- cbind(1, centered)
    cross <- crossprod(design)
    diag(cross)[-1L] <- diag(cross)[-1L] + data$lambda
    coefficient <- tryCatch(
      solve(cross, crossprod(design, data$y[index])),
      error = function(e) {
        qr.solve(cross, crossprod(design, data$y[index]))
      }
    )
    sum((data$y[index] - design %*% coefficient)^2)
  }, numeric(1))
  residual_df <- k - length(subset) - 1
  mean_sse <- mean(local_sse)
  list(
    risk = mean_sse / residual_df,
    mean_sse = mean_sse,
    residual_df = residual_df,
    contribution = local_sse / residual_df
  )
}
