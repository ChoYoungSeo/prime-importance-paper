source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/02_selection_to_pme/00_common.R")

methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
selection <- read.csv(
  "results/raw/selection_v052_runs.csv",
  stringsAsFactors = FALSE
)
tasks <- selection[
  selection$n == 100L &
    selection$design %in% c("independent", "balanced_proxy") &
    selection$replication <= 30L &
    selection$method %in% methods,
  ,
  drop = FALSE
]
tasks <- tasks[order(
  tasks$benchmark,
  tasks$design,
  tasks$replication,
  match(tasks$method, methods)
), , drop = FALSE]
tasks$link_id <- seq_len(nrow(tasks))
tasks <- tasks[tasks$benchmark == "ackley3", , drop = FALSE]

output_directory <- "results/raw/selection_to_pme_v052_main"
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(index) {
  task <- tasks[index, , drop = FALSE]
  checkpoint <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task$link_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    if (identical(attr(existing, "benchmark_definition"),
                  "ackley3_v1") &&
        identical(as.character(existing$row$method), task$method)) {
      return(existing)
    }
  }

  data <- simulate_selection_task(
    task$benchmark,
    task$n,
    task$p,
    task$design,
    task$rho,
    task$seed
  )
  truth <- selection_pme_reference(task$benchmark, task$p)
  selected <- parse_selected(task$selected)
  error <- ""
  start <- proc.time()[["elapsed"]]
  result <- tryCatch(
    selected_set_pme(data$X, data$y, selected, task$p),
    error = function(condition) {
      error <<- conditionMessage(condition)
      list(
        effects = rep(NA_real_, task$p),
        value_calls = NA_integer_,
        feasible = FALSE
      )
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start
  if (isFALSE(result$feasible) && !nzchar(error) &&
      is.character(result$reason) && length(result$reason) == 1L &&
      !is.na(result$reason) && nzchar(result$reason)) {
    error <- result$reason
  }
  valid <- result$feasible && all(is.finite(result$effects))
  origin <- if (task$method == "PRIME screen") {
    paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  } else {
    "paired Ackley-3 rerun"
  }
  row <- data.frame(
    link_id = task$link_id,
    task_id = task$task_id,
    benchmark = task$benchmark,
    n = task$n,
    p = task$p,
    design = task$design,
    rho = task$rho,
    replication = task$replication,
    seed = task$seed,
    method = task$method,
    selected = task$selected,
    selected_size = length(selected),
    active_inclusion = all(data$active %in% selected),
    feasible = result$feasible,
    value_calls = result$value_calls,
    half_l1 = if (valid) {
      0.5 * sum(abs(result$effects - truth))
    } else {
      NA_real_
    },
    inactive_mass = if (valid) {
      sum(result$effects[-data$active])
    } else {
      NA_real_
    },
    elapsed_sec = elapsed,
    error = error,
    result_origin = origin,
    stringsAsFactors = FALSE
  )
  phi <- data.frame(
    link_id = task$link_id,
    benchmark = task$benchmark,
    design = task$design,
    replication = task$replication,
    method = task$method,
    variable = seq_len(task$p),
    active = seq_len(task$p) %in% data$active,
    proxy = seq_len(task$p) %in% data$proxies,
    estimate = result$effects,
    truth = truth,
    result_origin = origin,
    stringsAsFactors = FALSE
  )
  object <- list(row = row, phi = phi)
  attr(object, "benchmark_definition") <- "ackley3_v1"
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
    seq_len(nrow(tasks)),
    run_task,
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  lapply(seq_len(nrow(tasks)), run_task)
}
ackley_rows <- do.call(rbind, lapply(objects, `[[`, "row"))
ackley_phi <- do.call(rbind, lapply(objects, `[[`, "phi"))

replace_rows <- function(path, replacement, key_order, keep_methods = NULL) {
  current <- read.csv(path, stringsAsFactors = FALSE)
  if (!is.null(keep_methods)) {
    replacement <- replacement[
      replacement$method %in% keep_methods, , drop = FALSE
    ]
  }
  missing_current <- setdiff(names(replacement), names(current))
  missing_replacement <- setdiff(names(current), names(replacement))
  current[missing_current] <- NA
  replacement[missing_replacement] <- NA
  updated <- rbind(
    current[
      !current$benchmark %in% c(
        "moon", "sobol_g3", "hartmann3", "ackley3", "currin2"
      ),
      names(current),
      drop = FALSE
    ],
    replacement[names(current)]
  )
  ordering <- do.call(
    order,
    c(unname(updated[key_order]), list(na.last = TRUE))
  )
  updated <- updated[ordering, , drop = FALSE]
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  verified_write_csv(updated, temporary, compressed = endsWith(path, ".gz"), row.names = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not replace ", path, call. = FALSE)
  }
}

replace_rows(
  file.path(output_directory, "selection_to_pme.csv"),
  ackley_rows,
  c("benchmark", "design", "replication", "method")
)
replace_rows(
  file.path(output_directory, "selection_to_pme_phi.csv"),
  ackley_phi,
  c("benchmark", "design", "replication", "method", "variable")
)
replace_rows(
  "results/provenance/paired_baselines/selection_to_pme_competitors.csv.gz",
  ackley_rows,
  c("benchmark", "design", "replication", "method"),
  setdiff(methods, "PRIME screen")
)
replace_rows(
  "results/provenance/paired_baselines/selection_to_pme_competitors_phi.csv.gz",
  ackley_phi,
  c("benchmark", "design", "replication", "method", "variable"),
  setdiff(methods, "PRIME screen")
)

mark_empty_infeasible <- function(run_path, phi_path) {
  runs <- read.csv(run_path, stringsAsFactors = FALSE)
  phi <- read.csv(phi_path, stringsAsFactors = FALSE)
  empty <- runs$selected_size == 0L
  if (any(empty)) {
    keys <- do.call(paste, c(runs[
      empty,
      c("benchmark", "design", "replication", "method")
    ], sep = "|"))
    phi_key <- do.call(paste, c(phi[
      c("benchmark", "design", "replication", "method")
    ], sep = "|"))
    runs$feasible[empty] <- FALSE
    runs$value_calls[empty] <- NA_integer_
    runs$half_l1[empty] <- NA_real_
    runs$inactive_mass[empty] <- NA_real_
    runs$error[empty] <- "empty selected set"
    phi$estimate[phi_key %in% keys] <- NA_real_
  }
  temporary_runs <- paste0(run_path, ".tmp-", Sys.getpid())
  temporary_phi <- paste0(phi_path, ".tmp-", Sys.getpid())
  verified_write_csv(runs, temporary_runs, compressed = endsWith(run_path, ".gz"), row.names = FALSE)
  verified_write_csv(phi, temporary_phi, compressed = endsWith(phi_path, ".gz"), row.names = FALSE)
  if (!file.rename(temporary_runs, run_path) ||
      !file.rename(temporary_phi, phi_path)) {
    unlink(c(temporary_runs, temporary_phi))
    stop("Could not normalize empty selected sets.", call. = FALSE)
  }
}

mark_empty_infeasible(
  file.path(output_directory, "selection_to_pme.csv"),
  file.path(output_directory, "selection_to_pme_phi.csv")
)
mark_empty_infeasible(
  "results/provenance/paired_baselines/selection_to_pme_competitors.csv.gz",
  "results/provenance/paired_baselines/selection_to_pme_competitors_phi.csv.gz"
)

canonical_link_map <- function(runs) {
  benchmark_order <- c("ackley3", "branin", "friedman")
  design_order <- c("balanced_proxy", "independent")
  ordering <- order(
    match(runs$benchmark, benchmark_order),
    match(runs$design, design_order),
    runs$replication,
    match(runs$method, methods)
  )
  ordered <- runs[ordering, c(
    "benchmark", "design", "replication", "method"
  )]
  if (anyDuplicated(ordered)) {
    stop("Duplicate selector-to-PME keys before link normalization.",
         call. = FALSE)
  }
  ordered$link_id <- seq_len(nrow(ordered))
  ordered
}

apply_link_map <- function(run_path, phi_path, mapping) {
  runs <- read.csv(run_path, stringsAsFactors = FALSE)
  phi <- read.csv(phi_path, stringsAsFactors = FALSE)
  key <- c("benchmark", "design", "replication", "method")
  run_key <- do.call(paste, c(runs[key], sep = "|"))
  phi_key <- do.call(paste, c(phi[key], sep = "|"))
  map_key <- do.call(paste, c(mapping[key], sep = "|"))
  runs$link_id <- mapping$link_id[match(run_key, map_key)]
  phi$link_id <- mapping$link_id[match(phi_key, map_key)]
  if (anyNA(runs$link_id) || anyNA(phi$link_id) ||
      anyDuplicated(runs$link_id)) {
    stop("Could not assign canonical selector-to-PME link IDs.",
         call. = FALSE)
  }
  runs <- runs[order(runs$link_id), , drop = FALSE]
  phi <- phi[order(phi$link_id, phi$variable), , drop = FALSE]
  temporary_runs <- paste0(run_path, ".tmp-", Sys.getpid())
  temporary_phi <- paste0(phi_path, ".tmp-", Sys.getpid())
  verified_write_csv(runs, temporary_runs, compressed = endsWith(run_path, ".gz"), row.names = FALSE)
  verified_write_csv(phi, temporary_phi, compressed = endsWith(phi_path, ".gz"), row.names = FALSE)
  if (!file.rename(temporary_runs, run_path) ||
      !file.rename(temporary_phi, phi_path)) {
    unlink(c(temporary_runs, temporary_phi))
    stop("Could not write canonical selector-to-PME links.", call. = FALSE)
  }
}

main_run_path <- file.path(output_directory, "selection_to_pme.csv")
main_phi_path <- file.path(output_directory, "selection_to_pme_phi.csv")
main_mapping <- canonical_link_map(read.csv(
  main_run_path, stringsAsFactors = FALSE
))
apply_link_map(main_run_path, main_phi_path, main_mapping)
apply_link_map(
  "results/provenance/paired_baselines/selection_to_pme_competitors.csv.gz",
  "results/provenance/paired_baselines/selection_to_pme_competitors_phi.csv.gz",
  main_mapping
)

write.csv(
  data.frame(
    benchmark = "ackley3",
    variable = 1:3,
    pme = selection_pme_reference("ackley3", 3L),
    source = "exact by permutation symmetry",
    stringsAsFactors = FALSE
  ),
  "results/summary/ackley3_pme_population_truth.csv",
  row.names = FALSE
)

cat("Refreshed", nrow(ackley_rows),
    "Ackley-3 selector-to-PME runs.\n")
