library(prime)

closed_to_dual <- function(closed) {
  closed <- as.numeric(closed)
  s <- as.integer(round(log2(length(closed))))
  full <- 2^s - 1L
  dual <- vapply(0:full, function(mask) {
    closed[full + 1L] - closed[bitwXor(full, mask) + 1L]
  }, numeric(1))
  dual / closed[full + 1L]
}

bratley_value <- function(subset) {
  subset <- sort(unique(as.integer(subset)))
  terms <- lapply(seq_len(5L), seq_len)
  coefficients <- (-1)^seq_len(5L)
  retained_terms <- lapply(terms, intersect, y = subset)
  integrated_count <- lengths(terms) - lengths(retained_terms)
  retained_coefficients <- coefficients * 0.5^integrated_count
  mean_value <- sum(
    retained_coefficients * 0.5^lengths(retained_terms)
  )
  second_moment <- 0
  for (left in seq_len(5L)) {
    for (right in seq_len(5L)) {
      left_set <- retained_terms[[left]]
      right_set <- retained_terms[[right]]
      shared <- length(intersect(left_set, right_set))
      single <- length(union(left_set, right_set)) - shared
      moment <- (1 / 3)^shared * (1 / 2)^single
      second_moment <- second_moment +
        retained_coefficients[left] *
        retained_coefficients[right] * moment
    }
  }
  second_moment - mean_value^2
}

reference_directory <- "results/provenance/reference_cache"
paths <- sort(list.files(
  reference_directory,
  pattern = "_cpfr_.*\\.rds$",
  full.names = TRUE
))
rows <- list()
objects <- list()
for (path in paths) {
  reference <- readRDS(path)
  if (!reference$benchmark %in% c(
    "ackley3", "branin", "friedman", "hartmann3", "ishigami", "otl", "piston",
    "borehole"
  )) {
    next
  }
  dual <- closed_to_dual(reference$coalition_values)
  effects <- prime_allocate(dual)
  key <- paste(
    reference$benchmark,
    paste0("rho", format(reference$rho, nsmall = 1)),
    sep = "_"
  )
  objects[[key]] <- list(
    benchmark = reference$benchmark,
    rho = reference$rho,
    active = reference$active,
    closed_values = reference$coalition_values,
    dual_game = dual,
    pme = as.numeric(effects),
    source = basename(path),
    n_ref = reference$n_ref,
    seed = reference$seed
  )
  rows[[length(rows) + 1L]] <- data.frame(
    benchmark = reference$benchmark,
    rho = reference$rho,
    variable = seq_along(effects),
    pme = as.numeric(effects),
    source = basename(path),
    n_ref = reference$n_ref,
    seed = reference$seed,
    stringsAsFactors = FALSE
  )
}

if (!is.null(objects$ackley3_rho0.0)) {
  objects$ackley3_rho0.0$pme <- rep(1 / 3, 3L)
  objects$ackley3_rho0.0$reference_note <-
    "exact PME by permutation symmetry"
  keep <- vapply(
    rows,
    function(row) !all(row$benchmark == "ackley3"),
    logical(1)
  )
  rows <- rows[keep]
  rows[[length(rows) + 1L]] <- data.frame(
    benchmark = "ackley3",
    rho = 0,
    variable = 1:3,
    pme = rep(1 / 3, 3L),
    source = paste0(
      objects$ackley3_rho0.0$source,
      "; exact PME by permutation symmetry"
    ),
    n_ref = objects$ackley3_rho0.0$n_ref,
    seed = objects$ackley3_rho0.0$seed,
    stringsAsFactors = FALSE
  )
}

closed_bratley <- vapply(0:(2^5 - 1L), function(mask) {
  subset <- which(as.logical(intToBits(mask)[seq_len(5L)]))
  bratley_value(subset)
}, numeric(1))
closed_bratley <- closed_bratley / closed_bratley[length(closed_bratley)]
dual_bratley <- closed_to_dual(closed_bratley)
pme_bratley <- prime_allocate(dual_bratley)
objects$bratley_rho0.0 <- list(
  benchmark = "bratley",
  rho = 0,
  active = 1:5,
  closed_values = closed_bratley,
  dual_game = dual_bratley,
  pme = as.numeric(pme_bratley),
  source = "analytic nested-product moments",
  n_ref = NA_integer_,
  seed = NA_integer_
)
rows[[length(rows) + 1L]] <- data.frame(
  benchmark = "bratley",
  rho = 0,
  variable = seq_along(pme_bratley),
  pme = as.numeric(pme_bratley),
  source = "analytic nested-product moments",
  n_ref = NA_integer_,
  seed = NA_integer_,
  stringsAsFactors = FALSE
)

sobol_a <- c(0, 0.5, 3, 9, 20, 50)
sobol_component_variance <- 1 / (3 * (1 + sobol_a)^2)
closed_sobol <- vapply(0:(2^length(sobol_a) - 1L), function(mask) {
  subset <- which(as.logical(intToBits(mask)[seq_along(sobol_a)]))
  if (!length(subset)) return(0)
  prod(1 + sobol_component_variance[subset]) - 1
}, numeric(1))
closed_sobol <- closed_sobol / closed_sobol[length(closed_sobol)]
dual_sobol <- closed_to_dual(closed_sobol)
pme_sobol <- prime_allocate(dual_sobol)
objects$sobol_g_rho0.0 <- list(
  benchmark = "sobol_g",
  rho = 0,
  active = seq_along(sobol_a),
  closed_values = closed_sobol,
  dual_game = dual_sobol,
  pme = as.numeric(pme_sobol),
  source = "analytic product moments",
  n_ref = NA_integer_,
  seed = NA_integer_
)
rows[[length(rows) + 1L]] <- data.frame(
  benchmark = "sobol_g",
  rho = 0,
  variable = seq_along(pme_sobol),
  pme = as.numeric(pme_sobol),
  source = "analytic product moments",
  n_ref = NA_integer_,
  seed = NA_integer_,
  stringsAsFactors = FALSE
)

sobol_g3_a <- c(0, 0.5, 1)
sobol_g3_component_variance <- 1 / (3 * (1 + sobol_g3_a)^2)
closed_sobol_g3 <- vapply(0:(2^length(sobol_g3_a) - 1L), function(mask) {
  subset <- which(as.logical(intToBits(mask)[seq_along(sobol_g3_a)]))
  if (!length(subset)) return(0)
  prod(1 + sobol_g3_component_variance[subset]) - 1
}, numeric(1))
closed_sobol_g3 <- closed_sobol_g3 / closed_sobol_g3[length(closed_sobol_g3)]
dual_sobol_g3 <- closed_to_dual(closed_sobol_g3)
pme_sobol_g3 <- prime_allocate(dual_sobol_g3)
objects$sobol_g3_rho0.0 <- list(
  benchmark = "sobol_g3",
  rho = 0,
  active = seq_along(sobol_g3_a),
  closed_values = closed_sobol_g3,
  dual_game = dual_sobol_g3,
  pme = as.numeric(pme_sobol_g3),
  source = "analytic balanced product moments",
  n_ref = NA_integer_,
  seed = NA_integer_
)
rows[[length(rows) + 1L]] <- data.frame(
  benchmark = "sobol_g3",
  rho = 0,
  variable = seq_along(pme_sobol_g3),
  pme = as.numeric(pme_sobol_g3),
  source = "analytic balanced product moments",
  n_ref = NA_integer_,
  seed = NA_integer_,
  stringsAsFactors = FALSE
)

saveRDS(objects, "results/provenance/pme_reference_objects.rds")
write.csv(
  do.call(rbind, rows),
  "results/summary/pme_reference_allocations.csv",
  row.names = FALSE
)
cat("Wrote", length(objects), "PME reference objects.\n")
