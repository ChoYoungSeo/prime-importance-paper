source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/lib/pme_pipelines.R")

profile <- Sys.getenv("PRIME_SCALING_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(benchmarks = "branin", p = 20L, replications = 1L),
  pilot = list(
    benchmarks = c("branin", "friedman", "otl", "piston", "borehole"),
    p = c(20L, 500L),
    replications = 1L
  ),
  main = list(
    benchmarks = c("branin", "friedman", "otl", "piston", "borehole"),
    p = c(20L, 50L, 100L, 250L, 500L),
    replications = 30L
  ),
  stop("Unknown scaling profile: ", profile, call. = FALSE)
)
methods <- complete_pipeline_methods()
methods_to_run <- if (profile == "main") methods[1L] else methods
tasks <- expand.grid(
  benchmark = configuration$benchmarks,
  p = configuration$p,
  replication = seq_len(configuration$replications),
  method = methods,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$n <- 500L
tasks$rho <- 0
tasks$task_id <- seq_len(nrow(tasks))
all_tasks <- tasks
run_tasks <- tasks[tasks$method %in% methods_to_run, , drop = FALSE]

output_directory <- file.path("results/raw", paste0("scaling_v052_", profile))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(run_index) {
  task <- run_tasks[run_index, , drop = FALSE]
  task_id <- task$task_id
  checkpoint <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    valid <- identical(as.character(existing$row$method), task$method) &&
      identical(as.integer(existing$row$task_id), as.integer(task_id))
    if (valid) return(existing)
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
      half_l1 = NA_real_,
      active_mae = NA_real_,
      inactive_mass = NA_real_,
      proxy_mass = NA_real_
    )
  }
  selected <- result$selected
  row <- data.frame(
    task_id = task_id,
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
    stringsAsFactors = FALSE
  )
  phi <- data.frame(
    task_id = task_id,
    benchmark = task$benchmark,
    p = task$p,
    replication = task$replication,
    method = task$method,
    variable = seq_len(task$p),
    active = seq_len(task$p) %in% data$active,
    estimate = result$effects,
    truth = truth,
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
    seq_len(nrow(run_tasks)),
    run_task,
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  lapply(seq_len(nrow(run_tasks)), run_task)
}
runs <- do.call(rbind, lapply(objects, `[[`, "row"))
phi <- do.call(rbind, lapply(objects, `[[`, "phi"))
if (profile == "main") {
  baseline_runs <- read.csv(
    "results/provenance/paired_baselines/scaling_competitors.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_phi <- read.csv(
    "results/provenance/paired_baselines/scaling_competitors_phi.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_runs <- baseline_runs[
    baseline_runs$method != methods[1L], , drop = FALSE
  ]
  baseline_phi <- baseline_phi[
    baseline_phi$method != methods[1L], , drop = FALSE
  ]
  expected <- all_tasks[all_tasks$method != methods[1L], , drop = FALSE]
  key <- c("benchmark", "p", "replication", "method")
  observed_keys <- do.call(paste, c(baseline_runs[key], sep = "|"))
  expected_keys <- do.call(paste, c(expected[key], sep = "|"))
  if (!setequal(observed_keys, expected_keys)) {
    stop("Reusable scaling baseline keys do not match.", call. = FALSE)
  }
  expected$seed <- mapply(
    scaling_seed, expected$benchmark, expected$p, expected$replication
  )
  seed_check <- merge(
    baseline_runs[c(key, "seed")], expected[c(key, "seed")],
    by = key, suffixes = c("_old", "_expected"), sort = FALSE
  )
  if (nrow(seed_check) != nrow(expected) || !isTRUE(all.equal(
    seed_check$seed_old, seed_check$seed_expected,
    check.attributes = FALSE
  ))) {
    stop("Reusable scaling baseline seeds do not match.", call. = FALSE)
  }
  runs$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  phi$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  baseline_runs$result_origin[
    is.na(baseline_runs$result_origin) |
      !nzchar(baseline_runs$result_origin)
  ] <- "reused paired baseline"
  baseline_phi$result_origin[
    is.na(baseline_phi$result_origin) |
      !nzchar(baseline_phi$result_origin)
  ] <- "reused paired baseline"
  runs <- rbind(runs[names(baseline_runs)], baseline_runs)
  phi <- rbind(phi[names(baseline_phi)], baseline_phi)
}
write.csv(
  runs[order(runs$benchmark, runs$p, runs$replication,
             match(runs$method, methods)), ],
  file.path(output_directory, "scaling_runs.csv"),
  row.names = FALSE
)
write.csv(
  phi[order(phi$benchmark, phi$p, phi$replication,
            match(phi$method, methods), phi$variable), ],
  file.path(output_directory, "scaling_phi.csv"),
  row.names = FALSE
)
cat("Completed", nrow(runs), "scaling runs.\n")
