source("experiments/lib/verified_io.R")
verified_require_context()
suppressPackageStartupMessages({
  library(prime)
  library(MASS)
})
source("experiments/01_selection/00_common.R")

if (utils::packageVersion("prime") < "0.5.2") {
  stop("Install the current PRIME package before running this comparison.",
       call. = FALSE)
}

prime_internal <- function(name) getFromNamespace(name, "prime")

clipped_value <- function(data, subset, k) {
  if (!length(subset)) return(0)
  risk <- prime_internal(".prime_risk")(data, subset, k)
  max(data$response_variance - risk, 0)
}

first_level_forward <- function(data) {
  ridge <- prime_internal(".prime_ridge_initial")(data)
  selected <- sort(unique(as.integer(ridge$selected)))
  dimension_limit <- prime_internal(".prime_dimension_limit")(
    data$n, data$p
  )

  repeat {
    if (length(selected) >= data$p ||
        length(selected) >= dimension_limit) break
    candidates <- setdiff(seq_len(data$p), selected)
    dimension <- length(selected) + 1L
    k <- prime_internal(".prime_neighborhood_size")(data$n, dimension)
    if (!is.finite(k)) break

    current <- clipped_value(data, selected, k)
    candidate_scores <- prime_internal(".prime_forward_values")(
      data, selected, candidates, k
    )
    candidate_scores <- pmax(candidate_scores * data$response_variance, 0)
    position <- which.max(candidate_scores)
    if (!length(position) || !is.finite(candidate_scores[position]) ||
        candidate_scores[position] <= current) break
    selected <- sort(c(selected, candidates[position]))
  }

  as.integer(selected)
}

first_level_backward <- function(data, candidates) {
  selected <- sort(unique(as.integer(candidates)))
  if (!length(selected)) return(integer())

  repeat {
    d <- length(selected)
    k <- prime_internal(".prime_neighborhood_size")(data$n, d)
    if (!is.finite(k)) {
      stop("Backward reduction is infeasible at dimension ", d, ".",
           call. = FALSE)
    }
    current <- clipped_value(data, selected, k)
    reduced_scores <- vapply(selected, function(variable) {
      clipped_value(data, setdiff(selected, variable), k)
    }, numeric(1))
    position <- which.max(reduced_scores)
    if (!length(position) || !is.finite(reduced_scores[position]) ||
        reduced_scores[position] <= current) break
    selected <- selected[-position]
    if (!length(selected)) break
  }

  as.integer(selected)
}

profile <- Sys.getenv("PRIME_LEVEL_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(n = c(100L, 500L), replications = 3L),
  pilot = list(n = c(100L, 200L, 500L, 1000L), replications = 20L),
  main = list(n = c(100L, 200L, 500L, 1000L, 2000L), replications = 50L),
  stop("Unknown PRIME_LEVEL_PROFILE: ", profile, call. = FALSE)
)

benchmarks <- c("branin", "ackley3", "friedman")
designs <- c("independent", "spaced_ar1", "balanced_proxy")
p <- 100L
rho <- 0.5
tasks <- expand.grid(
  benchmark = benchmarks,
  n = configuration$n,
  design = designs,
  replication = seq_len(configuration$replications),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$task_id <- seq_len(nrow(tasks))
tasks$seed <- 202609090L +
  1000000L * match(tasks$benchmark, benchmarks) +
  100000L * match(tasks$design, designs) +
  1000L * match(tasks$n, configuration$n) + tasks$replication

run_task <- function(task_id) {
  task <- tasks[task_id, , drop = FALSE]
  generated <- simulate_selection_data(
    task$benchmark, task$n, p, task$design, rho, task$seed
  )
  data <- prime_internal(".prime_prepare")(
    generated$X, generated$y, TRUE, NULL, 1e-6
  )

  start <- proc.time()[["elapsed"]]
  current <- prime_internal(".prime_forward_core")(data)
  current_selected <- prime_internal(".prime_backward_core")(
    data, current$candidates
  )$selected
  current_time <- proc.time()[["elapsed"]] - start

  start <- proc.time()[["elapsed"]]
  level_candidates <- first_level_forward(data)
  level_selected <- first_level_backward(data, level_candidates)
  level_time <- proc.time()[["elapsed"]] - start

  selections <- list(
    "buffered paired contrast" = list(
      candidates = current$candidates,
      selected = current_selected,
      elapsed = current_time
    ),
    "FIRST-style clipped level" = list(
      candidates = level_candidates,
      selected = level_selected,
      elapsed = level_time
    )
  )
  rows <- lapply(names(selections), function(method) {
    result <- selections[[method]]
    data.frame(
      task,
      p = p,
      rho = if (task$design == "independent") 0 else rho,
      method = method,
      candidate_size = length(result$candidates),
      elapsed_sec = result$elapsed,
      selection_metrics(result$selected, generated$active, p),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

cores <- min(
  8L,
  max(1L, suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "8"))))
)
Sys.setenv(PRIME_NUM_THREADS = "1")
results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    seq_len(nrow(tasks)), run_task,
    mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(seq_len(nrow(tasks)), run_task)
}
raw <- do.call(rbind, results)

output_directory <- file.path(
  "results", "raw", paste0("selection_calibration_v052_", profile)
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_directory, "clipped_level.csv")
write.csv(raw, raw_path, row.names = FALSE)

summary <- aggregate(
  cbind(
    candidate_size, selected_size, tpr, fpr, false_positives,
    inclusion, exact, elapsed_sec
  ) ~ method,
  raw,
  mean
)
write.csv(
  summary,
  file.path(output_directory, "clipped_level_summary.csv"),
  row.names = FALSE
)
print(summary, row.names = FALSE)
cat("Completed", nrow(tasks), "paired tasks.\n")
