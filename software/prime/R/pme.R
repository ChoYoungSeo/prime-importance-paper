.prime_mask_size <- function(mask, s) {
  sum(bitwAnd(mask, bitwShiftL(1L, 0:(s - 1L))) != 0L)
}

.prime_mask_contains <- function(mask, position) {
  bitwAnd(mask, bitwShiftL(1L, position - 1L)) != 0L
}

.prime_hasse_edges <- function(s) {
  full <- bitwShiftL(1L, s) - 1L
  lower <- integer(s * bitwShiftL(1L, s - 1L))
  upper <- integer(length(lower))
  position <- 0L
  for (mask in 0:full) {
    for (variable in seq_len(s)) {
      bit <- bitwShiftL(1L, variable - 1L)
      if (bitwAnd(mask, bit) == 0L) {
        position <- position + 1L
        lower[position] <- mask
        upper[position] <- bitwOr(mask, bit)
      }
    }
  }
  list(lower = lower, upper = upper)
}

.prime_project_monotone <- function(
  game,
  s,
  tolerance = 1e-10,
  max_iterations = 20000L
) {
  full <- bitwShiftL(1L, s) - 1L
  full_value <- game[full + 1L]
  if (!is.finite(full_value) || full_value < 0) {
    stop("The full-set value must be finite and nonnegative.",
         call. = FALSE)
  }
  if (s == 1L) {
    projected <- c(0, full_value)
    return(list(
      game = projected,
      iterations = 0L,
      converged = TRUE,
      max_violation = 0
    ))
  }

  interior_masks <- 1:(full - 1L)
  x <- as.numeric(game[interior_masks + 1L])
  variable_index <- integer(full + 1L)
  variable_index[interior_masks + 1L] <- seq_along(interior_masks)
  edges <- .prime_hasse_edges(s)
  edge_count <- length(edges$lower)
  correction <- numeric(edge_count)
  norm_squared <- numeric(edge_count)
  rhs <- numeric(edge_count)
  lower_index <- variable_index[edges$lower + 1L]
  upper_index <- variable_index[edges$upper + 1L]

  for (edge in seq_len(edge_count)) {
    norm_squared[edge] <- (lower_index[edge] > 0L) +
      (upper_index[edge] > 0L)
    fixed_lower <- if (edges$lower[edge] == 0L) 0 else 0
    fixed_upper <- if (edges$upper[edge] == full) full_value else 0
    rhs[edge] <- fixed_lower - fixed_upper
  }

  converged <- FALSE
  max_violation <- Inf
  for (iteration in seq_len(max_iterations)) {
    max_move <- 0
    for (edge in seq_len(edge_count)) {
      li <- lower_index[edge]
      ui <- upper_index[edge]
      dot_current <- 0
      if (ui > 0L) dot_current <- dot_current + x[ui]
      if (li > 0L) dot_current <- dot_current - x[li]
      dot_shifted <- dot_current +
        correction[edge] * norm_squared[edge]
      delta <- if (dot_shifted < rhs[edge]) {
        (rhs[edge] - dot_shifted) / norm_squared[edge]
      } else {
        0
      }
      move <- correction[edge] + delta
      if (ui > 0L) x[ui] <- x[ui] + move
      if (li > 0L) x[li] <- x[li] - move
      correction[edge] <- -delta
      max_move <- max(max_move, abs(move))
    }

    projected <- numeric(full + 1L)
    projected[1L] <- 0
    projected[interior_masks + 1L] <- x
    projected[full + 1L] <- full_value
    violations <- projected[edges$lower + 1L] -
      projected[edges$upper + 1L]
    max_violation <- max(0, max(violations))
    if (max_move <= tolerance && max_violation <= tolerance) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    stop("Monotone subset-value projection did not converge.", call. = FALSE)
  }
  projected[abs(projected) <= tolerance] <- 0
  projected[abs(projected - full_value) <= tolerance] <- full_value
  projected <- pmin(full_value, pmax(0, projected))
  list(
    game = projected,
    iterations = iteration,
    converged = converged,
    max_violation = max_violation
  )
}

.prime_ratio_potential_shifted <- function(game, s, zero_mask) {
  full <- bitwShiftL(1L, s) - 1L
  remaining <- bitwXor(full, zero_mask)
  potential <- numeric(full + 1L)
  potential[1L] <- 1
  admissible <- which(
    vapply(0:full, function(mask) {
      bitwAnd(mask, zero_mask) == 0L
    }, logical(1))
  ) - 1L
  sizes <- vapply(
    admissible,
    .prime_mask_size,
    integer(1),
    s = s
  )
  for (order in seq_len(.prime_mask_size(remaining, s))) {
    masks <- admissible[sizes == order]
    for (mask in masks) {
      predecessors <- integer()
      for (position in seq_len(s)) {
        bit <- bitwShiftL(1L, position - 1L)
        if (bitwAnd(mask, bit) != 0L) {
          predecessors <- c(
            predecessors,
            bitwXor(mask, bit)
          )
        }
      }
      denominator <- sum(1 / potential[predecessors + 1L])
      shifted_value <- game[bitwOr(zero_mask, mask) + 1L]
      potential[mask + 1L] <- shifted_value / denominator
    }
  }
  list(potential = potential, remaining = remaining)
}

.prime_extended_proportional_value <- function(game, s, zero_tolerance) {
  full <- bitwShiftL(1L, s) - 1L
  if (game[full + 1L] <= zero_tolerance) return(numeric(s))
  masks <- 0:full
  zero_masks <- masks[game <= zero_tolerance]
  sizes <- vapply(
    zero_masks,
    .prime_mask_size,
    integer(1),
    s = s
  )
  largest_zero <- zero_masks[sizes == max(sizes)]

  denominator <- 0
  numerator <- numeric(s)
  for (zero_mask in largest_zero) {
    ratio <- .prime_ratio_potential_shifted(game, s, zero_mask)
    full_potential <- ratio$potential[ratio$remaining + 1L]
    denominator <- denominator + 1 / full_potential
    for (position in seq_len(s)) {
      bit <- bitwShiftL(1L, position - 1L)
      if (bitwAnd(zero_mask, bit) == 0L) {
        without <- bitwXor(ratio$remaining, bit)
        numerator[position] <- numerator[position] +
          1 / ratio$potential[without + 1L]
      }
    }
  }
  effects <- numerator / denominator
  effects[abs(effects) <= zero_tolerance] <- 0
  effects
}

#' Allocate precomputed subset values
#'
#' This lower-level interface accepts values ordered by binary subset masks
#' from zero through eqn{2^d-1}. The returned allocation is the canonical
#' proportional-value extension used to define proportional marginal effects.
#' Computation and storage scale exponentially with the retained dimension.
#'
#' @param subset_values Numeric vector of subset values.
#' @param zero_tolerance Numerical threshold for null subsets.
#' @return A numeric PRIME allocation with the projected values and
#'   projection diagnostics attached as attributes.
#' @export
prime_allocate <- function(
  subset_values,
  zero_tolerance = sqrt(.Machine$double.eps)
) {
  game <- as.numeric(subset_values)
  if (!length(game) || any(!is.finite(game))) {
    stop(
      "subset_values must be a nonempty finite numeric vector.",
      call. = FALSE
    )
  }
  s <- log2(length(game))
  if (!isTRUE(all.equal(s, round(s))) || s < 1 ||
      s > .prime_max_mask_dimension) {
    stop(
      "The subset_values length must be 2^d with 1 <= d <= ",
      .prime_max_mask_dimension,
      " because subset masks use 32-bit integers.",
      call. = FALSE
    )
  }
  s <- as.integer(round(s))
  game[1L] <- 0
  projection <- .prime_project_monotone(game, s)
  projected <- projection$game
  threshold <- zero_tolerance * max(1, projected[length(projected)])
  projected[projected <= threshold] <- 0
  effects <- .prime_extended_proportional_value(
    projected, s, threshold
  )
  names(effects) <- paste0("X", seq_len(s))
  attr(effects, "projected_values") <- projected
  attr(effects, "projection") <- projection[
    c("iterations", "converged", "max_violation")
  ]
  effects
}

.prime_subset_label <- function(mask, selected, feature_names) {
  positions <- which(vapply(
    seq_along(selected),
    function(position) .prime_mask_contains(mask, position),
    logical(1)
  ))
  if (!length(positions)) return("{}")
  paste0("{", paste(feature_names[selected[positions]], collapse = ","), "}")
}

.prime_estimate_subset_values <- function(data, selected) {
  selected <- sort(unique(as.integer(selected)))
  s <- length(selected)
  if (!s) stop("At least one selected coordinate is required.",
               call. = FALSE)
  if (s > .prime_max_mask_dimension) {
    stop(
      "The selected dimension exceeds the 32-bit subset-mask limit of ",
      .prime_max_mask_dimension, "; the screen retained ", s, ".",
      call. = FALSE
    )
  }
  full_mask <- bitwShiftL(1L, s) - 1L
  k_full <- .prime_neighborhood_size(data$n, s)
  noise_floor <- .prime_risk(data, selected, k_full)
  signal <- max(data$response_variance - noise_floor, 0)
  signal_guard <- sqrt(.Machine$double.eps) * data$response_variance
  if (!is.finite(signal) || signal <= signal_guard) {
    stop("The retained coordinates have no estimable explained variance.",
         call. = FALSE)
  }

  raw <- numeric(full_mask + 1L)
  complement_risk <- numeric(full_mask + 1L)
  adjusted_risk <- numeric(full_mask + 1L)
  for (mask in 0:full_mask) {
    keep_positions <- which(!vapply(
      seq_len(s),
      function(position) .prime_mask_contains(mask, position),
      logical(1)
    ))
    keep <- selected[keep_positions]
    if (!length(keep)) {
      risk_keep <- data$response_variance
    } else {
      k <- .prime_neighborhood_size(data$n, length(keep))
      risk_keep <- .prime_risk(data, keep, k)
    }
    complement_risk[mask + 1L] <- risk_keep
    adjusted_risk[mask + 1L] <- max(risk_keep - noise_floor, 0)
    raw[mask + 1L] <- adjusted_risk[mask + 1L] / signal
  }
  raw[1L] <- 0
  raw[full_mask + 1L] <- 1
  labels <- vapply(
    0:full_mask,
    .prime_subset_label,
    character(1),
    selected = selected,
    feature_names = data$feature_names
  )
  list(
    raw = raw,
    subset = labels,
    complement_risk = complement_risk,
    adjusted_risk = adjusted_risk,
    signal = signal,
    noise_floor = noise_floor,
    risk_full = noise_floor
  )
}
