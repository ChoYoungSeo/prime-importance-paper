source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/lib/pme_pipelines.R")

methods <- complete_pipeline_methods()
new_methods <- c(
  "LASSO screen + PME-kNN",
  "FOCI screen + PME-kNN",
  "FORD screen + PME-kNN"
)
benchmarks <- c("branin", "friedman", "otl", "piston", "borehole")
tasks <- expand.grid(
  benchmark = benchmarks,
  p = c(20L, 50L, 100L, 250L, 500L),
  replication = seq_len(30L),
  method = methods,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$n <- 500L
tasks$rho <- 0
tasks$task_id <- seq_len(nrow(tasks))
tasks <- tasks[tasks$method %in% new_methods, , drop = FALSE]

output_directory <- "results/raw/scaling_v052_main"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(index) {
  task <- tasks[index, , drop = FALSE]
  checkpoint <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    if (identical(as.character(existing$row$method), task$method) &&
        identical(as.integer(existing$row$task_id),
                  as.integer(task$task_id))) {
      return(existing)
    }
  }

  seed <- scaling_seed(task$benchmark, task$p, task$replication)
  data <- simulate_benchmark(
    task$benchmark, task$n, task$p, rho = 0, seed = seed
  )
  truth <- reference_pme(task$benchmark, task$p)
  result <- run_complete_pipeline(
    task$method, data$X, data$y, seed + 10000L
  )
  valid <- result$feasible && all(is.finite(result$effects))
  metrics <- if (valid) {
    pme_metrics(result$effects, truth, data$active, integer())
  } else {
    data.frame(
      half_l1 = NA_real_, active_mae = NA_real_,
      inactive_mass = NA_real_, proxy_mass = NA_real_
    )
  }
  selected <- result$selected
  row <- data.frame(
    task_id = task$task_id,
    benchmark = task$benchmark,
    n = task$n,
    p = task$p,
    replication = task$replication,
    seed = seed,
    method = task$method,
    selected = paste(selected, collapse = ","),
    selected_size = length(selected),
    active_inclusion = all(data$active %in% selected),
    exact_active_set = identical(selected, data$active),
    false_positives = length(setdiff(selected, data$active)),
    feasible = result$feasible,
    selector_calls = result$selector_calls,
    value_calls = result$value_calls,
    selection_elapsed = result$selection_elapsed,
    allocation_elapsed = result$allocation_elapsed,
    elapsed_sec = result$elapsed_sec,
    metrics,
    error = result$error,
    result_origin = "paired baseline rerun",
    stringsAsFactors = FALSE
  )
  phi <- data.frame(
    task_id = task$task_id,
    benchmark = task$benchmark,
    p = task$p,
    replication = task$replication,
    method = task$method,
    variable = seq_len(task$p),
    active = seq_len(task$p) %in% data$active,
    estimate = result$effects,
    truth = truth,
    result_origin = "paired baseline rerun",
    stringsAsFactors = FALSE
  )
  object <- list(row = row, phi = phi)
  atomic_save_rds(object, checkpoint)
  object
}

cores <- min(
  6L,
  suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
)
if (!is.finite(cores) || cores < 1L) cores <- 1L
Sys.setenv(PRIME_NUM_THREADS = "1")
objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    seq_len(nrow(tasks)), run_task,
    mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(seq_len(nrow(tasks)), run_task)
}
new_runs <- do.call(rbind, lapply(objects, `[[`, "row"))
new_phi <- do.call(rbind, lapply(objects, `[[`, "phi"))

atomic_write_csv <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  verified_write_csv(object, temporary, compressed = endsWith(path, ".gz"), row.names = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not replace ", path, call. = FALSE)
  }
}

run_path <- file.path(output_directory, "scaling_runs.csv")
phi_path <- file.path(output_directory, "scaling_phi.csv")
current_runs <- read.csv(run_path, stringsAsFactors = FALSE)
current_phi <- read.csv(phi_path, stringsAsFactors = FALSE)
all_runs <- rbind(
  current_runs[!current_runs$method %in% new_methods, ],
  new_runs[names(current_runs)]
)
all_phi <- rbind(
  current_phi[!current_phi$method %in% new_methods, ],
  new_phi[names(current_phi)]
)
all_runs <- all_runs[order(
  all_runs$benchmark, all_runs$p, all_runs$replication,
  match(all_runs$method, methods)
), ]
all_phi <- all_phi[order(
  all_phi$benchmark, all_phi$p, all_phi$replication,
  match(all_phi$method, methods), all_phi$variable
), ]
stopifnot(
  nrow(all_runs) == length(benchmarks) * 5L * 30L * length(methods),
  !anyDuplicated(all_runs[c("benchmark", "p", "replication", "method")]),
  !anyDuplicated(all_phi[c(
    "benchmark", "p", "replication", "method", "variable"
  )])
)
atomic_write_csv(all_runs, run_path)
atomic_write_csv(all_phi, phi_path)
atomic_write_csv(
  all_runs[all_runs$method != methods[1L], ],
  "results/provenance/paired_baselines/scaling_competitors.csv.gz"
)
atomic_write_csv(
  all_phi[all_phi$method != methods[1L], ],
  "results/provenance/paired_baselines/scaling_competitors_phi.csv.gz"
)
cat("Added", nrow(new_runs), "scaling competitor runs.\n")
