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

raw_gain_forward <- function(data) {
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

    values <- prime_internal(".prime_forward_values")(
      data, selected, candidates, k
    )
    position <- which.max(values)
    if (!length(position) || !is.finite(values[position])) {
      stop("No finite forward candidate value was available.", call. = FALSE)
    }
    variable <- candidates[position]
    diagnostic <- prime_internal(".prime_forward_diagnostic")(
      data, selected, variable, k
    )
    if (!is.finite(diagnostic$raw_gain) || diagnostic$raw_gain <= 0) break
    selected <- sort(c(selected, variable))
  }

  as.integer(selected)
}

positive_part_backward <- function(data, candidates) {
  selected <- sort(unique(as.integer(candidates)))
  if (!length(selected)) return(integer())

  repeat {
    d <- length(selected)
    k <- prime_internal(".prime_neighborhood_size")(data$n, d)
    if (!is.finite(k)) {
      stop("Backward reduction is infeasible at dimension ", d, ".",
           call. = FALSE)
    }
    full <- prime_internal(".prime_variance_details")(data, selected, k)
    raw_loss <- vapply(selected, function(variable) {
      prime_internal(".prime_deletion_diagnostic")(
        data, selected, variable, k, full_details = full
      )$raw_loss
    }, numeric(1))
    clipped <- pmax(raw_loss, 0)
    position <- which.min(clipped)
    if (!length(position) || !is.finite(clipped[position]) ||
        clipped[position] > 0) break
    selected <- selected[-position]
    if (!length(selected)) break
  }

  as.integer(selected)
}

profile <- Sys.getenv("PRIME_RULE_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(n = c(100L, 500L), replications = 3L),
  pilot = list(n = c(100L, 200L, 500L, 1000L), replications = 20L),
  main = list(n = c(100L, 200L, 500L, 1000L, 2000L), replications = 50L),
  stop("Unknown PRIME_RULE_PROFILE: ", profile, call. = FALSE)
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

  buffered_candidates <- prime_internal(".prime_forward_core")(data)$candidates
  positive_candidates <- raw_gain_forward(data)
  candidate_sets <- list(
    "buffered forward" = buffered_candidates,
    "positive-part forward" = positive_candidates
  )

  rows <- list()
  position <- 0L
  for (forward_name in names(candidate_sets)) {
    candidates <- candidate_sets[[forward_name]]
    selected_sets <- list(
      "buffered backward" = prime_internal(".prime_backward_core")(
        data, candidates
      )$selected,
      "positive-part backward" = positive_part_backward(data, candidates)
    )
    for (backward_name in names(selected_sets)) {
      selected <- selected_sets[[backward_name]]
      metrics <- selection_metrics(selected, generated$active, p)
      position <- position + 1L
      rows[[position]] <- data.frame(
        task,
        p = p,
        rho = if (task$design == "independent") 0 else rho,
        active = paste(generated$active, collapse = ","),
        forward_rule = forward_name,
        backward_rule = backward_name,
        method = paste(forward_name, backward_name, sep = " + "),
        candidate_size = length(candidates),
        metrics,
        stringsAsFactors = FALSE
      )
    }
  }
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
raw_path <- file.path(
  output_directory, "buffer_ablation.csv"
)
write.csv(raw, raw_path, row.names = FALSE)

summary <- aggregate(
  cbind(
    candidate_size, selected_size, tpr, fpr, false_positives,
    inclusion, exact
  ) ~ method,
  raw,
  mean
)
summary_path <- file.path(
  output_directory,
  "buffer_ablation_summary.csv"
)
write.csv(summary, summary_path, row.names = FALSE)

cell_summary <- aggregate(
  cbind(selected_size, tpr, fpr, false_positives, inclusion, exact) ~
    benchmark + n + design + method,
  raw,
  mean
)
cell_path <- file.path(
  output_directory,
  "buffer_ablation_cells.csv"
)
write.csv(cell_summary, cell_path, row.names = FALSE)

print(summary, row.names = FALSE)
cat("Completed", nrow(tasks), "paired tasks.\n")
cat("Raw:", raw_path, "\n")
cat("Summary:", summary_path, "\n")
