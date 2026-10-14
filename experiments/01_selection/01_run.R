source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/01_selection/00_common.R")

profile <- Sys.getenv("PRIME_SELECTION_PROFILE", "smoke")
configuration <- selection_profile(profile)
methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
methods_to_run <- if (profile == "main") "PRIME screen" else methods
p <- 100L
rho <- 0.5

method_packages <- c(
  "PRIME screen" = "MASS",
  "FIRST" = "first",
  "LASSO" = "glmnet",
  "FOCI" = "FOCI",
  "FORD" = "FORD"
)
required <- unique(unname(method_packages[methods_to_run]))
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

tasks <- expand.grid(
  benchmark = configuration$benchmarks,
  n = configuration$n,
  design = configuration$designs,
  replication = seq_len(configuration$replications),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
tasks$task_id <- seq_len(nrow(tasks))
tasks$seed <- 202608080L +
  1000000L * match(tasks$benchmark, configuration$benchmarks) +
  100000L * match(tasks$design, configuration$designs) +
  1000L * match(tasks$n, configuration$n) + tasks$replication
if (anyDuplicated(tasks$seed)) {
  stop("Selection task seeds must be unique.", call. = FALSE)
}

checkpoint_directory <- file.path(
  "results/raw",
  paste0("selection_v052_", profile, "_checkpoints")
)
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(task_id) {
  checkpoint <- file.path(
    checkpoint_directory,
    sprintf("task_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    keep <- existing$method %in% methods_to_run
    valid <- sum(keep) == length(methods_to_run) &&
      setequal(existing$method[keep], methods_to_run) &&
      all(existing$task_id[keep] == task_id)
    if (valid) return(existing[keep, , drop = FALSE])
  }
  task <- tasks[task_id, , drop = FALSE]
  data <- simulate_selection_data(
    task$benchmark,
    task$n,
    p,
    task$design,
    rho,
    task$seed
  )
  rows <- lapply(methods_to_run, function(method) {
    fit <- run_selection_method(method, data$X, data$y, task$seed)
    metrics <- selection_metrics(fit$selected, data$active, p)
    data.frame(
      task,
      p = p,
      rho = if (task$design == "independent") 0 else rho,
      active = paste(data$active, collapse = ","),
      proxies = paste(data$proxies, collapse = ","),
      method = method,
      metrics,
      elapsed_sec = fit$elapsed_sec,
      error = fit$error,
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      prime_version = as.character(utils::packageVersion("prime")),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  selection_atomic_save(result, checkpoint)
  result
}

cores <- min(
  8L,
  suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "8")))
)
task_order <- order(-tasks$n, tasks$task_id)
results <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    task_order,
    run_task,
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  lapply(task_order, run_task)
}

current <- do.call(rbind, results)
current <- current[current$method == "PRIME screen" | profile != "main", ]
if (profile == "main") {
  baseline <- read.csv(
    "results/provenance/paired_baselines/selection_competitors.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline <- baseline[
    baseline$method %in% setdiff(methods, "PRIME screen"),
    ,
    drop = FALSE
  ]
  expected_baseline <- nrow(tasks) * (length(methods) - 1L)
  if (nrow(baseline) != expected_baseline) {
    stop(
      "Expected ", expected_baseline,
      " reusable baseline rows; found ", nrow(baseline),
      call. = FALSE
    )
  }
  baseline_keys <- paste(
    baseline$task_id, baseline$benchmark, baseline$n, baseline$design,
    baseline$replication, baseline$seed, sep = "|"
  )
  expected_keys <- rep(
    paste(
      tasks$task_id, tasks$benchmark, tasks$n, tasks$design,
      tasks$replication, tasks$seed, sep = "|"
    ),
    each = length(methods) - 1L
  )
  if (!setequal(baseline_keys, expected_keys)) {
    stop("Reusable baseline task keys do not match current tasks.",
         call. = FALSE)
  }
  current$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  baseline$result_origin[
    is.na(baseline$result_origin) | !nzchar(baseline$result_origin)
  ] <- "reused paired baseline"
  columns <- unique(c(names(current), names(baseline)))
  current[setdiff(columns, names(current))] <- NA
  baseline[setdiff(columns, names(baseline))] <- NA
  raw <- rbind(
    current[columns],
    baseline[columns]
  )
} else {
  raw <- current
  raw$result_origin <- "current profile run"
}
raw <- raw[order(raw$task_id, match(raw$method, methods)), , drop = FALSE]
output <- selection_output_path(profile)
write.csv(raw, output, row.names = FALSE)
cat("Completed", nrow(tasks), "selection tasks and", nrow(raw),
    "method runs.\n")
