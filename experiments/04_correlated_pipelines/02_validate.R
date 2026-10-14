profile <- Sys.getenv("PRIME_CORRELATED_PROFILE", "smoke")
directory <- file.path(
  "results/raw", paste0("correlated_pipelines_v052_", profile)
)
runs <- read.csv(
  file.path(directory, "correlated_pipeline_runs.csv"),
  stringsAsFactors = FALSE
)
phi <- read.csv(
  file.path(directory, "correlated_pipeline_phi.csv"),
  stringsAsFactors = FALSE
)
replications <- switch(profile, smoke = 1L, pilot = 3L, main = 20L)
n_values <- switch(profile, smoke = 100L, pilot = c(100L, 500L),
                   main = c(100L, 200L, 500L))
source("experiments/lib/pme_pipelines.R")
methods <- complete_pipeline_methods()
expected_runs <- 2L * length(n_values) * replications * length(methods)
stopifnot(
  nrow(runs) == expected_runs,
  nrow(phi) == expected_runs * 20L,
  !anyDuplicated(runs[c("benchmark", "n", "replication", "method")]),
  !anyDuplicated(phi[c(
    "benchmark", "n", "replication", "method", "variable"
  )]),
  setequal(unique(runs$method), methods),
  all(runs$selected_size >= 0L),
  all(runs$value_calls[runs$feasible] >= 0L),
  all(is.finite(runs$half_l1[runs$feasible])),
  all(is.na(runs$half_l1[!runs$feasible])),
  all(is.na(runs$error[runs$feasible]) |
      runs$error[runs$feasible] == "")
)
if (profile == "main") {
  proposed_name <- "PRIME"
  current_origins <- c(
    paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  )
  stopifnot(
    all(runs$result_origin[runs$method == proposed_name] %in%
        current_origins),
    all(runs$result_origin[runs$method != proposed_name] %in%
        c("reused paired baseline", "paired baseline rerun")),
    all(phi$result_origin[phi$method == proposed_name] %in%
        current_origins),
    all(phi$result_origin[phi$method != proposed_name] %in%
        c("reused paired baseline", "paired baseline rerun"))
  )
}

effect_sums <- aggregate(
  estimate ~ benchmark + n + replication + method,
  phi,
  function(x) sum(x, na.rm = TRUE)
)
truth_sums <- aggregate(
  truth ~ benchmark + n + replication + method,
  phi,
  sum
)
coordinate_counts <- aggregate(
  variable ~ benchmark + n + replication + method,
  phi,
  length
)
effect_sums <- merge(
  effect_sums,
  runs[c(
    "benchmark", "n", "replication", "method",
    "selected_size", "feasible"
  )],
  by = c("benchmark", "n", "replication", "method"),
  all.x = TRUE
)
nonempty <- effect_sums$feasible & effect_sums$selected_size > 0L
empty <- effect_sums$feasible & effect_sums$selected_size == 0L
stopifnot(
  all(coordinate_counts$variable == 20L),
  max(abs(truth_sums$truth - 1)) < 1e-7,
  max(abs(effect_sums$estimate[nonempty] - 1)) < 1e-7,
  !any(empty) || max(abs(effect_sums$estimate[empty])) < 1e-12
)
proposed <- runs$method == "PRIME"
if (any(!runs$feasible[proposed])) {
  stop("The proposed pipeline has an infeasible run.", call. = FALSE)
}
cat("Correlated-pipeline validation: PASS\n")
print(aggregate(
  cbind(half_l1, active_inclusion, exact_active_set,
        proxy_count, feasible, elapsed_sec) ~ benchmark + n + method,
  runs,
  function(x) mean(x, na.rm = TRUE)
))
