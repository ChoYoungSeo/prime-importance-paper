profile <- Sys.getenv("PRIME_OBSERVED_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(datasets = 1L, replications = 1L, controls = 5L),
  main = list(datasets = 3L, replications = 20L, controls = 20L),
  stop("Unknown observed-data profile: ", profile, call. = FALSE)
)
directory <- file.path(
  "results/raw", paste0("observed_v052_", profile)
)
runs <- read.csv(
  file.path(directory, "observed_runs.csv"), stringsAsFactors = FALSE
)
phi <- read.csv(
  file.path(directory, "observed_phi.csv"), stringsAsFactors = FALSE
)
expected_runs <- configuration$datasets * configuration$replications
stopifnot(
  nrow(runs) == expected_runs,
  !anyDuplicated(runs[c("dataset", "replication")]),
  all(runs$control_p == configuration$controls),
  all(runs$selected_size ==
      runs$original_selected + runs$controls_selected),
  all(runs$controls_selected >= 0L),
  all(runs$control_mass[runs$feasible] >= 0),
  all(runs$control_mass[runs$feasible] <= 1),
  all(is.na(runs$error[runs$feasible]) |
      runs$error[runs$feasible] == "")
)
coordinate_counts <- aggregate(
  variable ~ dataset + replication, phi, length
)
expected_counts <- merge(
  coordinate_counts,
  runs[c("dataset", "replication", "original_p", "control_p")],
  by = c("dataset", "replication")
)
stopifnot(all(
  expected_counts$variable ==
    expected_counts$original_p + expected_counts$control_p
))
effect_sums <- aggregate(
  estimate ~ dataset + replication, phi, sum, na.rm = TRUE
)
effect_sums <- merge(
  effect_sums,
  runs[c("dataset", "replication", "feasible", "selected_size")],
  by = c("dataset", "replication")
)
nonempty <- effect_sums$feasible & effect_sums$selected_size > 0L
empty <- effect_sums$feasible & effect_sums$selected_size == 0L
stopifnot(
  !any(nonempty) || max(abs(effect_sums$estimate[nonempty] - 1)) < 1e-7,
  !any(empty) || max(abs(effect_sums$estimate[empty])) < 1e-12
)
cat("Observed-data validation: PASS\n")
print(aggregate(
  cbind(original_selected, controls_selected, control_mass,
        elapsed_sec) ~ dataset,
  runs,
  function(x) c(mean = mean(x, na.rm = TRUE),
                median = median(x, na.rm = TRUE))
))
