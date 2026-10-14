source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/02_selection_to_pme/00_common.R")
selection_helpers <- new.env(parent = globalenv())
sys.source(
  "experiments/01_selection/00_common.R",
  envir = selection_helpers
)

profile <- Sys.getenv("PRIME_FRIEDMAN_DETAIL_PROFILE", "smoke")
replications <- switch(
  profile,
  smoke = 2L,
  main = 30L,
  stop("Unknown Friedman-detail profile: ", profile, call. = FALSE)
)
methods <- c("PRIME screen", "FIRST", "FOCI", "FORD")
selection <- if (profile == "main") {
  read.csv(
    "results/raw/selection_v052_runs.csv",
    stringsAsFactors = FALSE
  )
} else {
  rows <- lapply(seq_len(replications), function(replication) {
    seed <- 202609020L + replication
    data <- simulate_selection_task(
      "friedman", 200L, 100L, "balanced_proxy", 0.5, seed
    )
    do.call(rbind, lapply(methods, function(method) {
      fit <- selection_helpers$run_selection_method(
        method, data$X, data$y, seed
      )
      data.frame(
        task_id = replication,
        benchmark = "friedman",
        n = 200L,
        p = 100L,
        design = "balanced_proxy",
        rho = 0.5,
        replication = replication,
        seed = seed,
        method = method,
        selection_helpers$selection_metrics(
          fit$selected, data$active, 100L
        ),
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}
tasks <- selection[
  selection$benchmark == "friedman" &
    selection$design == "balanced_proxy" &
    selection$n == 200 &
    selection$replication <= replications &
    selection$method %in% methods,
  ,
  drop = FALSE
]
tasks <- tasks[order(
  tasks$replication, match(tasks$method, methods)
), ]
tasks$detail_id <- seq_len(nrow(tasks))
expected <- replications * length(methods)
if (nrow(tasks) != expected) {
  stop("Expected ", expected, " Friedman-detail tasks; found ", nrow(tasks),
       call. = FALSE)
}

output_directory <- file.path(
  "results/raw", paste0("friedman_detail_v052_", profile)
)
checkpoint_directory <- file.path(output_directory, "checkpoints")
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(detail_id) {
  task <- tasks[tasks$detail_id == detail_id, , drop = FALSE]
  checkpoint <- file.path(
    checkpoint_directory, sprintf("task_%04d.rds", detail_id)
  )
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  data <- simulate_selection_task(
    "friedman", task$n, task$p, "balanced_proxy", task$rho, task$seed
  )
  selected <- parse_selected(task$selected)
  truth <- selection_pme_reference("friedman", task$p)
  error <- ""
  result <- tryCatch(
    selected_set_pme(data$X, data$y, selected, task$p),
    error = function(condition) {
      error <<- conditionMessage(condition)
      list(
        effects = rep(NA_real_, task$p), value_calls = NA_integer_,
        feasible = FALSE
      )
    }
  )
  valid <- isTRUE(result$feasible) && all(is.finite(result$effects))
  row <- data.frame(
    detail_id = detail_id,
    task_id = task$task_id,
    replication = task$replication,
    seed = task$seed,
    method = task$method,
    selected = task$selected,
    selected_size = task$selected_size,
    active_inclusion = task$inclusion,
    exact = task$exact,
    feasible = valid,
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
    error = error,
    stringsAsFactors = FALSE
  )
  phi <- data.frame(
    detail_id = detail_id,
    replication = task$replication,
    method = task$method,
    variable = seq_len(task$p),
    role = ifelse(
      seq_len(task$p) %in% data$active,
      "active",
      ifelse(
        seq_len(task$p) %in% data$proxies,
        "correlated proxy",
        "independent noise"
      )
    ),
    estimate = result$effects,
    truth = truth,
    stringsAsFactors = FALSE
  )
  object <- list(row = row, phi = phi)
  atomic_save_rds(object, checkpoint)
  object
}

cores <- suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
if (!length(cores) || !is.finite(cores) || cores < 1L) cores <- 1L
Sys.setenv(PRIME_NUM_THREADS = "1")
objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    tasks$detail_id, run_task,
    mc.cores = min(cores, 6L), mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(tasks$detail_id, run_task)
}
runs <- do.call(rbind, lapply(objects, function(object) object$row))
phi <- do.call(rbind, lapply(objects, function(object) object$phi))
runs <- runs[order(runs$replication, match(runs$method, methods)), ]
phi <- phi[order(
  phi$replication, match(phi$method, methods), phi$variable
), ]
write.csv(
  runs, file.path(output_directory, "friedman_detail_runs.csv"),
  row.names = FALSE
)
write.csv(
  phi, file.path(output_directory, "friedman_detail_phi.csv"),
  row.names = FALSE
)
cat("Completed", nrow(runs), "Friedman-detail runs.\n")
