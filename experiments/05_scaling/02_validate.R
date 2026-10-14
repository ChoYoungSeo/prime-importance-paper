profile <- Sys.getenv("PRIME_SCALING_PROFILE", "smoke")
directory <- file.path("results/raw", paste0("scaling_v052_", profile))
runs <- read.csv(
  file.path(directory, "scaling_runs.csv"),
  stringsAsFactors = FALSE
)
phi <- read.csv(
  file.path(directory, "scaling_phi.csv"),
  stringsAsFactors = FALSE
)
configuration <- switch(
  profile,
  smoke = list(benchmarks = 1L, p = 1L, replications = 1L),
  pilot = list(benchmarks = 5L, p = 2L, replications = 1L),
  main = list(benchmarks = 5L, p = 5L, replications = 30L)
)
source("experiments/lib/pme_pipelines.R")
methods <- complete_pipeline_methods()
expected <- configuration$benchmarks * configuration$p *
  configuration$replications * length(methods)
stopifnot(
  nrow(runs) == expected,
  !anyDuplicated(runs[c("benchmark", "p", "replication", "method")]),
  !anyDuplicated(phi[c(
    "benchmark", "p", "replication", "method", "variable"
  )]),
  setequal(unique(runs$method), methods),
  all(runs$selected_size >= 0L),
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
  estimate ~ benchmark + p + replication + method,
  phi,
  function(x) sum(x, na.rm = TRUE)
)
truth_sums <- aggregate(
  truth ~ benchmark + p + replication + method,
  phi,
  sum
)
coordinate_counts <- aggregate(
  variable ~ benchmark + p + replication + method,
  phi,
  length
)
effect_sums <- merge(
  effect_sums,
  runs[c(
    "benchmark", "p", "replication", "method",
    "selected_size", "feasible"
  )],
  by = c("benchmark", "p", "replication", "method"),
  all.x = TRUE
)
nonempty <- effect_sums$feasible & effect_sums$selected_size > 0L
empty <- effect_sums$feasible & effect_sums$selected_size == 0L
stopifnot(
  all(coordinate_counts$variable == coordinate_counts$p),
  max(abs(truth_sums$truth - 1)) < 1e-7,
  max(abs(effect_sums$estimate[nonempty] - 1)) < 1e-7,
  !any(empty) || max(abs(effect_sums$estimate[empty])) < 1e-12
)
proposed <- runs$method == "PRIME"
if (any(!runs$feasible[proposed])) {
  stop("The proposed scaling pipeline has an infeasible run.", call. = FALSE)
}
cat("Scaling validation: PASS\n")
print(aggregate(
  cbind(half_l1, active_inclusion, exact_active_set,
        false_positives, feasible, elapsed_sec) ~ benchmark + p + method,
  runs,
  function(x) mean(x, na.rm = TRUE)
))
