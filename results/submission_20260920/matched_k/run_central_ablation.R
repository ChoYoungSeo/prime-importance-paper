args <- commandArgs(trailingOnly = TRUE)
limit <- if (length(args)) as.integer(args[1]) else Inf

root <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
plan <- normalizePath(file.path(root, ".."))
kit <- file.path(plan, "execution_20260919")
input_dir <- file.path(kit, "results/4f29d1c54aff9649/checkpoints")
checkpoint_dir <- file.path(root, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)
.libPaths(c(file.path(kit, "library"), file.path(Sys.getenv("HOME"), "R/library"), .libPaths()))
stopifnot(
  normalizePath(find.package("prime")) == normalizePath(file.path(kit, "library/prime")),
  as.character(utils::packageVersion("prime")) == "0.5.2.9101"
)
suppressPackageStartupMessages(library(prime))
setwd(kit)
source("experiments/03_pme/00_common.R")
setwd(root)

atomic_rds <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(value, temporary, version = 3L)
  stopifnot(file.rename(temporary, path))
}

atomic_json <- function(value, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  jsonlite::write_json(
    value, temporary, auto_unbox = TRUE, pretty = TRUE,
    digits = 16, na = "null", null = "null"
  )
  stopifnot(file.rename(temporary, path))
}

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
}

population <- utils::read.csv(
  file.path(plan, "allocation_rule_comparison_20260920/population/population_allocations.csv"),
  stringsAsFactors = FALSE
)
population <- population[population$rho == 0.5, ]
stopifnot(nrow(population) == 10L)
pme_truth <- population$PME
shapley_truth <- population$Shapley_full10

inputs <- sort(list.files(
  input_dir,
  pattern = "^friedman_n200_rho0\\.5_r[0-9]{3}\\.input\\.rds$",
  full.names = TRUE
))
stopifnot(length(inputs) == 30L)
if (is.finite(limit)) inputs <- head(inputs, limit)

records <- lapply(inputs, function(input_path) {
  id <- sub("\\.input\\.rds$", "", basename(input_path))
  output_path <- file.path(checkpoint_dir, paste0(id, ".matched_k_lc.rds"))
  json_path <- sub("\\.rds$", ".json", output_path)
  if (file.exists(output_path) && file.exists(json_path)) {
    return(readRDS(output_path)$record)
  }
  input <- readRDS(input_path)
  stopifnot(
    input$task$id == id,
    input$task$n == 200L,
    input$task$p == 10L,
    identical(input$task$rho, 0.5)
  )
  gc()
  started <- proc.time()[["elapsed"]]
  result <- matched_k_local_constant_game(input$X, input$y)
  elapsed <- proc.time()[["elapsed"]] - started
  record <- list(
    task_id = id,
    replication = as.integer(input$task$replication),
    seed = as.integer(input$task$seed),
    n = 200L,
    p = 10L,
    rho = 0.5,
    status = "complete",
    method = "Matched-K local-constant",
    neighborhood_rule = "K_n(d)=ceil(max(1.5*sqrt(n),3*d)), capped at floor(n/2)",
    k_by_dimension = result$k_by_dimension,
    value_calls = result$value_calls,
    elapsed_sec = elapsed,
    pme = result$pme,
    shapley = result$shapley,
    pme_half_l1 = 0.5 * sum(abs(result$pme - pme_truth)),
    shapley_half_l1 = 0.5 * sum(abs(result$shapley - shapley_truth)),
    pme_proxy_mass = sum(result$pme[6:10]),
    shapley_proxy_mass = sum(result$shapley[6:10]),
    input_md5 = unname(tools::md5sum(input_path)),
    prime_version = as.character(utils::packageVersion("prime"))
  )
  atomic_rds(list(record = record, details = result), output_path)
  atomic_json(record, json_path)
  record
})

summary <- do.call(rbind, lapply(records, function(record) {
  data.frame(
    task_id = record$task_id,
    replication = record$replication,
    pme_half_l1 = record$pme_half_l1,
    shapley_half_l1 = record$shapley_half_l1,
    pme_proxy_mass = record$pme_proxy_mass,
    shapley_proxy_mass = record$shapley_proxy_mass,
    elapsed_sec = record$elapsed_sec,
    stringsAsFactors = FALSE
  )
}))
summary <- summary[order(summary$replication), ]
utils::write.csv(summary, file.path(root, "central_summary.csv"), row.names = FALSE)

receipt <- list(
  status = if (nrow(summary) == length(inputs) && all(is.finite(as.matrix(summary[, -c(1, 2)])))) "PASS" else "FAIL",
  scope = "Friedman n=200, p=10, rho=0.5; frozen 30 paired datasets",
  datasets = nrow(summary),
  estimator = "full-coordinate local-constant PME and Shapley with the same dimension-dependent K_n(d) rule as PRIME",
  neighbor_identity = "K nearest other observations; query self-excluded",
  k_by_dimension = records[[1]]$k_by_dimension,
  mean = as.list(colMeans(summary[, c("pme_half_l1", "shapley_half_l1", "pme_proxy_mass", "shapley_proxy_mass", "elapsed_sec")])),
  median_elapsed_sec = stats::median(summary$elapsed_sec),
  input_md5_unique = length(unique(vapply(records, `[[`, character(1), "input_md5"))) == nrow(summary),
  prime_version = as.character(utils::packageVersion("prime")),
  script_md5 = unname(tools::md5sum(file.path(root, "run_central_ablation.R")))
)
atomic_json(receipt, file.path(root, "RECEIPT.json"))
print(receipt)
stopifnot(identical(receipt$status, "PASS"))
