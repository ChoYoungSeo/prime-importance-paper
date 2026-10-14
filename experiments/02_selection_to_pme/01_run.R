source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/02_selection_to_pme/00_common.R")

profile <- Sys.getenv("PRIME_SELECTION_LINK_PROFILE", "smoke")
replications <- switch(
  profile,
  smoke = 1L,
  pilot = 3L,
  main = 30L,
  stop("Unknown profile: ", profile, call. = FALSE)
)
methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
methods_to_run <- if (profile == "main") "PRIME screen" else methods
selection_file <- switch(
  profile,
  smoke = "results/raw/selection_v052_smoke.csv",
  pilot = "results/raw/selection_v052_pilot.csv",
  main = "results/raw/selection_v052_runs.csv"
)
selection <- read.csv(
  selection_file,
  stringsAsFactors = FALSE
)
tasks <- selection[
  selection$n == 100L &
    selection$design %in% c("independent", "balanced_proxy") &
    selection$replication <= replications &
    selection$method %in% methods,
  ,
  drop = FALSE
]
tasks <- tasks[
  order(
    tasks$benchmark,
    tasks$design,
    tasks$replication,
    match(tasks$method, methods)
  ),
  ,
  drop = FALSE
]
tasks$link_id <- seq_len(nrow(tasks))
expected_benchmarks <- 3L
expected <- expected_benchmarks * 2L * replications * length(methods)
if (nrow(tasks) != expected) {
  stop("Expected ", expected, " tasks; found ", nrow(tasks), call. = FALSE)
}
all_tasks <- tasks
run_tasks <- tasks[tasks$method %in% methods_to_run, , drop = FALSE]

output_directory <- file.path(
  "results/raw", paste0("selection_to_pme_v052_", profile)
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(run_index) {
  task <- run_tasks[run_index, , drop = FALSE]
  link_id <- task$link_id
  checkpoint <- file.path(
    output_directory,
    sprintf("task_%04d.rds", link_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    valid <- identical(as.character(existing$row$method), task$method) &&
      identical(as.integer(existing$row$link_id), as.integer(link_id)) &&
      identical(as.integer(existing$row$task_id), as.integer(task$task_id)) &&
      identical(as.character(existing$row$benchmark), task$benchmark) &&
      identical(as.character(existing$row$design), task$design) &&
      identical(as.integer(existing$row$replication),
                as.integer(task$replication)) &&
      identical(as.integer(existing$row$seed), as.integer(task$seed))
    if (valid) return(existing)
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
  row <- data.frame(
    link_id = link_id,
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
    stringsAsFactors = FALSE
  )
  phi <- data.frame(
    link_id = link_id,
    benchmark = task$benchmark,
    design = task$design,
    replication = task$replication,
    method = task$method,
    variable = seq_len(task$p),
    active = seq_len(task$p) %in% data$active,
    proxy = seq_len(task$p) %in% data$proxies,
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
results <- if (.Platform$OS.type == "unix" && cores > 1L) {
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
bind_fill <- function(objects) {
  columns <- unique(unlist(lapply(objects, names), use.names = FALSE))
  objects <- lapply(objects, function(object) {
    object[setdiff(columns, names(object))] <- NA
    object[columns]
  })
  do.call(rbind, objects)
}
rows <- bind_fill(lapply(results, `[[`, "row"))
phi <- bind_fill(lapply(results, `[[`, "phi"))
if (profile == "main") {
  baseline_rows <- read.csv(
    "results/provenance/paired_baselines/selection_to_pme_competitors.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_phi <- read.csv(
    "results/provenance/paired_baselines/selection_to_pme_competitors_phi.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_rows <- baseline_rows[
    baseline_rows$method != "PRIME screen" &
      baseline_rows$benchmark %in% unique(all_tasks$benchmark),
    , drop = FALSE
  ]
  baseline_phi <- baseline_phi[
    baseline_phi$method != "PRIME screen" &
      baseline_phi$benchmark %in% unique(all_tasks$benchmark),
    , drop = FALSE
  ]
  expected_baseline <- all_tasks[
    all_tasks$method != "PRIME screen", , drop = FALSE
  ]
  key <- c("benchmark", "design", "replication", "method")
  observed_keys <- do.call(paste, c(baseline_rows[key], sep = "|"))
  expected_keys <- do.call(paste, c(expected_baseline[key], sep = "|"))
  if (!setequal(observed_keys, expected_keys)) {
    stop("Reusable selector-to-PME baseline keys do not match.",
         call. = FALSE)
  }
  baseline_check <- merge(
    baseline_rows[c(key, "seed", "selected")],
    expected_baseline[c(key, "seed", "selected")],
    by = key,
    suffixes = c("_old", "_new"),
    sort = FALSE
  )
  selected_old <- ifelse(
    is.na(baseline_check$selected_old), "", baseline_check$selected_old
  )
  selected_new <- ifelse(
    is.na(baseline_check$selected_new), "", baseline_check$selected_new
  )
  if (nrow(baseline_check) != nrow(expected_baseline) ||
      !isTRUE(all.equal(
        baseline_check$seed_old, baseline_check$seed_new,
        check.attributes = FALSE
      )) || !identical(as.character(selected_old), as.character(selected_new))) {
    stop("Reusable selector-to-PME baselines are not paired exactly.",
         call. = FALSE)
  }
  rows$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  phi$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  baseline_rows$result_origin[
    is.na(baseline_rows$result_origin) |
      !nzchar(baseline_rows$result_origin)
  ] <- "reused paired baseline"
  baseline_phi$result_origin[
    is.na(baseline_phi$result_origin) |
      !nzchar(baseline_phi$result_origin)
  ] <- "reused paired baseline"
  rows <- rbind(rows[names(baseline_rows)], baseline_rows)
  phi <- rbind(phi[names(baseline_phi)], baseline_phi)
}
rows <- rows[order(
  rows$benchmark, rows$design, rows$replication,
  match(rows$method, methods)
), , drop = FALSE]
phi <- phi[order(
  phi$benchmark, phi$design, phi$replication,
  match(phi$method, methods), phi$variable
), , drop = FALSE]
write.csv(
  rows,
  file.path(output_directory, "selection_to_pme.csv"),
  row.names = FALSE
)
write.csv(
  phi,
  file.path(output_directory, "selection_to_pme_phi.csv"),
  row.names = FALSE
)
cat("Completed", nrow(rows), "selector-to-PME tasks.\n")
