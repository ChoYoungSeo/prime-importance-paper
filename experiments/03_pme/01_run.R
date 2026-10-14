source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/03_pme/00_common.R")

profile <- Sys.getenv("PRIME_PROFILE", "smoke")
specification <- switch(
  profile,
  smoke = list(
    benchmarks = "branin",
    n = 100L,
    rho = c(0, 0.9),
    replications = 1L,
    methods = c(
      "PRIME", "Full local-linear PME", "PME-kNN",
      "NCMD-PME", "Oracle local-linear PME"
    )
  ),
  pilot = list(
    benchmarks = c("branin", "bratley"),
    n = c(100L, 250L),
    rho = c(0, 0.9),
    replications = 3L,
    methods = c(
      "PRIME", "Full local-linear PME", "PME-kNN",
      "NCMD-PME", "Oracle local-linear PME"
    )
  ),
  main = list(
    benchmarks = c("branin", "bratley"),
    n = c(100L, 250L, 500L),
    rho = c(0, 0.5, 0.9),
    replications = 30L,
    methods = c(
      "PRIME", "Full local-linear PME", "PME-kNN",
      "NCMD-PME", "Oracle local-linear PME"
    )
  ),
  stop("Unknown profile: ", profile, call. = FALSE)
)
methods <- specification$methods
methods_to_run <- if (profile == "main") {
  "PRIME"
} else {
  methods
}

tasks <- expand.grid(
  benchmark = specification$benchmarks,
  n = specification$n,
  rho = specification$rho,
  replication = seq_len(specification$replications),
  stringsAsFactors = FALSE
)
output_directory <- file.path("results/raw", paste0("pme_v052_", profile))
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

run_task <- function(task_id) {
  checkpoint <- file.path(
    output_directory,
    sprintf("task_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) {
    existing <- readRDS(checkpoint)
    keep <- existing$rows$method %in% methods_to_run
    valid <- sum(keep) == length(methods_to_run) &&
      setequal(existing$rows$method[keep], methods_to_run) &&
      all(existing$rows$task_id[keep] == task_id)
    if (valid) {
      return(list(
        rows = existing$rows[keep, , drop = FALSE],
        phi = existing$phi[existing$phi$method %in% methods_to_run,
                           , drop = FALSE]
      ))
    }
  }
  task <- tasks[task_id, ]
  benchmark_id <- match(task$benchmark, names(proxy_registry()))
  seed <- 20260816L + 100000L * benchmark_id +
    1000L * task$n + 100L * round(10 * task$rho) +
    task$replication
  data <- simulate_proxy_data(
    task$benchmark, task$n, task$rho, seed
  )
  truth <- reference_pme(task$benchmark, data$p)
  rows <- list()
  phi <- list()
  for (run_method_id in seq_along(methods_to_run)) {
    method <- methods_to_run[run_method_id]
    method_id <- match(method, methods)
    result <- run_pme_method(
      method, data, seed + 10000L * method_id
    )
    metrics <- if (all(is.finite(result$effects))) {
      pme_metrics(result$effects, truth, data$active, data$proxies)
    } else {
      data.frame(
        half_l1 = NA_real_,
        active_mae = NA_real_,
        inactive_mass = NA_real_,
        proxy_mass = NA_real_
      )
    }
    rows[[run_method_id]] <- data.frame(
      task_id = task_id,
      benchmark = task$benchmark,
      n = task$n,
      p = data$p,
      rho = task$rho,
      replication = task$replication,
      seed = seed,
      method = method,
      selected = paste(result$selected, collapse = ","),
      selected_size = length(result$selected),
      value_calls = result$value_calls,
      elapsed_sec = result$elapsed_sec,
      metrics,
      error = result$error,
      stringsAsFactors = FALSE
    )
    phi[[run_method_id]] <- data.frame(
      task_id = task_id,
      benchmark = task$benchmark,
      n = task$n,
      rho = task$rho,
      replication = task$replication,
      method = method,
      variable = seq_len(data$p),
      active = seq_len(data$p) %in% data$active,
      proxy = seq_len(data$p) %in% data$proxies,
      estimate = result$effects,
      truth = truth,
      stringsAsFactors = FALSE
    )
  }
  object <- list(rows = do.call(rbind, rows), phi = do.call(rbind, phi))
  atomic_save_rds(object, checkpoint)
  object
}

cores <- min(
  6L,
  suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
)
if (!is.finite(cores) || cores < 1L) cores <- 1L
if (.Platform$OS.type == "unix" && cores > 1L) {
  Sys.setenv(PRIME_NUM_THREADS = "1")
  objects <- parallel::mclapply(
    seq_len(nrow(tasks)),
    run_task,
    mc.cores = cores,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  objects <- lapply(seq_len(nrow(tasks)), run_task)
}
rows <- do.call(rbind, lapply(objects, `[[`, "rows"))
phi <- do.call(rbind, lapply(objects, `[[`, "phi"))
if (profile == "main") {
  baseline_rows <- read.csv(
    "results/provenance/paired_baselines/pme_controls.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_phi <- read.csv(
    "results/provenance/paired_baselines/pme_controls_phi.csv.gz",
    stringsAsFactors = FALSE
  )
  baseline_rows <- baseline_rows[
    baseline_rows$method != "PRIME", , drop = FALSE
  ]
  baseline_phi <- baseline_phi[
    baseline_phi$method != "PRIME", , drop = FALSE
  ]
  key <- c("benchmark", "n", "rho", "replication", "method")
  expected <- expand.grid(
    benchmark = specification$benchmarks,
    n = specification$n,
    rho = specification$rho,
    replication = seq_len(specification$replications),
    method = setdiff(methods, "PRIME"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  expected$seed <- 20260816L +
    100000L * match(expected$benchmark, names(proxy_registry())) +
    1000L * expected$n + 100L * round(10 * expected$rho) +
    expected$replication
  observed_keys <- do.call(paste, c(baseline_rows[key], sep = "|"))
  expected_keys <- do.call(paste, c(expected[key], sep = "|"))
  if (!setequal(observed_keys, expected_keys)) {
    stop("Reusable low-dimensional PME baseline keys do not match.",
         call. = FALSE)
  }
  seed_check <- merge(
    baseline_rows[c(key, "seed")], expected[c(key, "seed")],
    by = key, suffixes = c("_old", "_expected"), sort = FALSE
  )
  if (nrow(seed_check) != nrow(expected) || !isTRUE(all.equal(
    seed_check$seed_old, seed_check$seed_expected,
    check.attributes = FALSE
  ))) {
    stop("Reusable low-dimensional PME seeds do not match.", call. = FALSE)
  }
  rows$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  phi$result_origin <- paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  baseline_rows$result_origin <- "reused paired baseline"
  baseline_phi$result_origin <- "reused paired baseline"
  rows <- rbind(rows[names(baseline_rows)], baseline_rows)
  phi <- rbind(phi[names(baseline_phi)], baseline_phi)
}
rows <- rows[order(
  rows$benchmark, rows$n, rows$rho, rows$replication,
  match(rows$method, methods)
), , drop = FALSE]
phi <- phi[order(
  phi$benchmark, phi$n, phi$rho, phi$replication,
  match(phi$method, methods), phi$variable
), , drop = FALSE]
write.csv(
  rows,
  file.path(output_directory, "pme_runs.csv"),
  row.names = FALSE
)
write.csv(
  phi,
  file.path(output_directory, "pme_phi.csv"),
  row.names = FALSE
)
