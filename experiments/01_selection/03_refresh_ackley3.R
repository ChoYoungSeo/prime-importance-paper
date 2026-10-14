source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/01_selection/00_common.R")

required_packages <- c("first", "glmnet", "FOCI", "FORD")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(
    "Missing required selector packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

configuration <- selection_profile("main")
methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
p <- 100L
rho <- 0.5

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
tasks <- tasks[tasks$benchmark == "ackley3", , drop = FALSE]

checkpoint_directory <-
  "results/raw/selection_v052_main_checkpoints"
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)

package_version_or_na <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

run_task <- function(index) {
  task <- tasks[index, , drop = FALSE]
  checkpoint <- file.path(
    checkpoint_directory,
    sprintf("task_%04d.rds", task$task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    if (identical(attr(existing, "benchmark_definition"),
                  "ackley3_v1") &&
        setequal(existing$method, methods) &&
        !any(nzchar(ifelse(is.na(existing$error), "", existing$error)))) {
      return(existing)
    }
  }

  data <- simulate_selection_data(
    task$benchmark, task$n, p, task$design, rho, task$seed
  )
  rows <- lapply(methods, function(method) {
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
      firstllr_version = package_version_or_na("firstllr"),
      ford_version = package_version_or_na("FORD"),
      result_origin = if (method == "PRIME screen") {
        paste("prime", as.character(utils::packageVersion("prime")), "rerun")
      } else {
        "paired Ackley-3 rerun"
      },
      prime_version = if (method == "PRIME screen") {
        as.character(utils::packageVersion("prime"))
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  attr(result, "benchmark_definition") <- "ackley3_v1"
  selection_atomic_save(result, checkpoint)
  result
}

cores <- min(
  6L,
  suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
)
if (!is.finite(cores) || cores < 1L) cores <- 1L
Sys.setenv(PRIME_NUM_THREADS = "1")
order_index <- order(-tasks$n, tasks$task_id)
objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    order_index,
    run_task,
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  lapply(order_index, run_task)
}
ackley3 <- do.call(rbind, objects)
ackley3 <- ackley3[
  order(ackley3$task_id, match(ackley3$method, methods)),
]

replace_benchmark_rows <- function(path, replacement, keep_methods = NULL) {
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
      !current$benchmark %in% c("moon", "sobol_g3", "hartmann3", "ackley3"),
      names(current),
      drop = FALSE
    ],
    replacement[names(current)]
  )
  updated <- updated[order(
    updated$task_id, match(updated$method, methods)
  ), , drop = FALSE]
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  verified_write_csv(updated, temporary, compressed = endsWith(path, ".gz"), row.names = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not replace ", path, call. = FALSE)
  }
}

replace_benchmark_rows(
  "results/raw/selection_v052_runs.csv",
  ackley3
)
replace_benchmark_rows(
  "results/provenance/paired_baselines/selection_competitors.csv.gz",
  ackley3,
  setdiff(methods, "PRIME screen")
)

cat("Refreshed", nrow(tasks), "Ackley-3 datasets and", nrow(ackley3),
    "method runs.\n")
