source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/03_pme/00_common.R")
source("experiments/09_observed_data/00_common.R")

profile <- Sys.getenv("PRIME_OBSERVED_PROFILE", "smoke")
benchmark_exact_dimension_budget <- 12L
configuration <- switch(
  profile,
  smoke = list(datasets = "airfoil", replications = 1L, controls = 5L),
  main = list(
    datasets = c("concrete", "airfoil", "auto_mpg"),
    replications = 20L,
    controls = 20L
  ),
  stop("Unknown observed-data profile: ", profile, call. = FALSE)
)

tasks <- expand.grid(
  dataset = configuration$datasets,
  replication = seq_len(configuration$replications),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$task_id <- seq_len(nrow(tasks))
tasks$seed <- 202608190L +
  100000L * match(tasks$dataset, configuration$datasets) +
  tasks$replication

output_directory <- file.path(
  "results/raw", paste0("observed_v052_", profile)
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(task_id) {
  task <- tasks[task_id, , drop = FALSE]
  checkpoint <- file.path(
    output_directory, sprintf("task_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    if (identical(as.integer(existing$row$task_id), as.integer(task_id))) {
      return(existing)
    }
  }

  data <- load_observed_data(task$dataset)
  augmented <- inject_permuted_controls(
    data$X, configuration$controls, task$seed
  )
  error <- ""
  selection_start <- proc.time()[["elapsed"]]
  screen <- tryCatch(
    prime_screen(augmented$X, data$y),
    error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    }
  )
  selection_elapsed <- proc.time()[["elapsed"]] - selection_start
  selected <- if (is.null(screen)) integer() else screen$selected

  allocation_start <- proc.time()[["elapsed"]]
  allocation <- if (
    is.null(screen) || length(selected) > benchmark_exact_dimension_budget
  ) {
    if (!is.null(screen) &&
        length(selected) > benchmark_exact_dimension_budget) {
      error <- paste0("selected dimension ", length(selected),
                      " exceeds the benchmark exact-enumeration budget of ",
                      benchmark_exact_dimension_budget)
    }
    list(
      effects = rep(NA_real_, ncol(augmented$X)),
      value_calls = NA_integer_,
      feasible = FALSE
    )
  } else if (!length(selected)) {
    error <- "empty_screen: allocation unavailable"
    list(
      effects = rep(NA_real_, ncol(augmented$X)),
      value_calls = NA_integer_,
      feasible = FALSE
    )
  } else {
    tryCatch({
      result <- estimate_pme_with_value(
        make_llr_value(augmented$X, data$y),
        selected, ncol(augmented$X)
      )
      result$feasible <- TRUE
      result
    }, error = function(condition) {
      error <<- conditionMessage(condition)
      list(
        effects = rep(NA_real_, ncol(augmented$X)),
        value_calls = NA_integer_,
        feasible = FALSE
      )
    })
  }
  allocation_elapsed <- proc.time()[["elapsed"]] - allocation_start
  valid <- allocation$feasible && all(is.finite(allocation$effects))

  row <- data.frame(
    task_id = task_id,
    dataset = task$dataset,
    dataset_label = data$label,
    replication = task$replication,
    seed = task$seed,
    n = nrow(data$X),
    original_p = ncol(data$X),
    control_p = configuration$controls,
    selected = paste(selected, collapse = ","),
    selected_size = length(selected),
    original_selected = sum(selected %in% augmented$original),
    controls_selected = sum(selected %in% augmented$controls),
    control_mass = if (valid) {
      sum(allocation$effects[augmented$controls])
    } else {
      NA_real_
    },
    feasible = allocation$feasible,
    value_calls = allocation$value_calls,
    selection_elapsed = selection_elapsed,
    allocation_elapsed = allocation_elapsed,
    elapsed_sec = selection_elapsed + allocation_elapsed,
    error = error,
    stringsAsFactors = FALSE
  )
  role <- rep("permuted control", ncol(augmented$X))
  role[augmented$original] <- "original"
  phi <- data.frame(
    task_id = task_id,
    dataset = task$dataset,
    dataset_label = data$label,
    replication = task$replication,
    variable = seq_len(ncol(augmented$X)),
    variable_name = colnames(augmented$X),
    role = role,
    selected = seq_len(ncol(augmented$X)) %in% selected,
    estimate = allocation$effects,
    stringsAsFactors = FALSE
  )
  object <- list(row = row, phi = phi)
  atomic_save_observed(object, checkpoint)
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
    tasks$task_id, run_task,
    mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(tasks$task_id, run_task)
}
rows <- do.call(rbind, lapply(objects, `[[`, "row"))
phi <- do.call(rbind, lapply(objects, `[[`, "phi"))
rows <- rows[order(rows$dataset, rows$replication), , drop = FALSE]
phi <- phi[order(
  phi$dataset, phi$replication, phi$variable
), , drop = FALSE]
write.csv(
  rows, file.path(output_directory, "observed_runs.csv"), row.names = FALSE
)
write.csv(
  phi, file.path(output_directory, "observed_phi.csv"), row.names = FALSE
)
cat("Completed", nrow(rows), "observed-data runs.\n")
