source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/01_selection/00_common.R")

profile <- Sys.getenv("PRIME_TWINNING_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(
    benchmarks = "branin", p = 20L, rho = 0,
    replications = 1L, twin_r = c(NA_integer_, 3L)
  ),
  main = list(
    benchmarks = c("branin", "ackley3", "friedman"),
    p = c(20L, 100L), rho = c(0, 0.5),
    replications = 5L, twin_r = c(NA_integer_, 2L, 3L, 5L)
  ),
  stop("Unknown Data-Twinning profile: ", profile, call. = FALSE)
)

tasks <- expand.grid(
  benchmark = configuration$benchmarks,
  p = configuration$p,
  rho = configuration$rho,
  replication = seq_len(configuration$replications),
  twin_r = configuration$twin_r,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$n <- 500L
tasks$task_id <- seq_len(nrow(tasks))
tasks$seed <- 202608170L +
  1000000L * match(tasks$benchmark, configuration$benchmarks) +
  10000L * match(tasks$p, configuration$p) +
  1000L * match(tasks$rho, configuration$rho) +
  tasks$replication

output_directory <- file.path(
  "results/raw", paste0("twinning_v052_", profile)
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

collapse_set <- function(x) paste(sort(unique(as.integer(x))), collapse = ",")

jaccard <- function(left, right) {
  left <- sort(unique(as.integer(left)))
  right <- sort(unique(as.integer(right)))
  denominator <- length(union(left, right))
  if (!denominator) return(1)
  length(intersect(left, right)) / denominator
}

atomic_save <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save checkpoint: ", path, call. = FALSE)
  }
}

run_task <- function(task_id) {
  task <- tasks[task_id, , drop = FALSE]
  checkpoint <- file.path(
    output_directory, sprintf("task_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    definition_valid <- task$benchmark != "ackley3" ||
      identical(attr(existing, "benchmark_definition"),
                "ackley3_v1")
    if (identical(as.integer(existing$task_id), as.integer(task_id)) &&
        identical(as.character(existing$benchmark), task$benchmark) &&
        definition_valid) {
      return(existing)
    }
  }
  design <- if (task$rho == 0) "independent" else "spaced_ar1"
  data <- simulate_selection_data(
    task$benchmark, task$n, task$p, design, task$rho, task$seed
  )
  ratio <- if (is.na(task$twin_r)) NULL else as.integer(task$twin_r)
  error <- ""
  start <- Sys.time()
  fit <- tryCatch(
    prime_screen(data$X, data$y, twin_r = ratio),
    error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  row <- data.frame(
    task_id = task_id,
    benchmark = task$benchmark,
    p = task$p,
    n = task$n,
    rho = task$rho,
    design = design,
    replication = task$replication,
    seed = task$seed,
    twin_r = if (is.null(ratio)) "all" else as.character(ratio),
    query_count = if (is.null(fit)) NA_integer_ else
      length(fit$data$query_index),
    candidate = if (is.null(fit)) "" else collapse_set(fit$candidates),
    selected = if (is.null(fit)) "" else collapse_set(fit$selected),
    selected_size = if (is.null(fit)) NA_integer_ else length(fit$selected),
    active_recall = if (is.null(fit)) NA_real_ else
      mean(data$active %in% fit$selected),
    exact_recovery = if (is.null(fit)) NA else
      identical(sort(fit$selected), sort(data$active)),
    elapsed_sec = elapsed,
    error = error,
    stringsAsFactors = FALSE
  )
  if (task$benchmark == "ackley3") {
    attr(row, "benchmark_definition") <- "ackley3_v1"
  }
  atomic_save(row, checkpoint)
  row
}

cores <- min(
  6L,
  suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
)
if (!is.finite(cores) || cores < 1L) cores <- 1L
Sys.setenv(PRIME_NUM_THREADS = "1")
rows <- if (.Platform$OS.type == "unix" && cores > 1L) {
  do.call(rbind, parallel::mclapply(
    tasks$task_id, run_task,
    mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
  ))
} else {
  do.call(rbind, lapply(tasks$task_id, run_task))
}

key <- c("benchmark", "p", "rho", "replication", "seed")
all_query <- rows[rows$twin_r == "all", c(
  key, "candidate", "selected", "elapsed_sec"
)]
names(all_query)[-(seq_along(key))] <- c(
  "candidate_all", "selected_all", "elapsed_all"
)
rows <- merge(rows, all_query, by = key, all.x = TRUE, sort = FALSE)
parse_set <- function(value) {
  if (is.na(value) || !nzchar(value)) return(integer())
  as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
}
rows$candidate_jaccard <- mapply(
  function(x, y) jaccard(parse_set(x), parse_set(y)),
  rows$candidate, rows$candidate_all
)
rows$selected_jaccard <- mapply(
  function(x, y) jaccard(parse_set(x), parse_set(y)),
  rows$selected, rows$selected_all
)
rows$speedup <- rows$elapsed_all / rows$elapsed_sec
rows <- rows[order(
  rows$benchmark, rows$p, rows$rho, rows$replication,
  match(rows$twin_r, c("all", "2", "3", "5"))
), , drop = FALSE]
write.csv(
  rows,
  file.path(output_directory, "twinning_runs.csv"),
  row.names = FALSE
)
cat("Completed", nrow(rows), "Data-Twinning runs.\n")
