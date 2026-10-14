profile <- Sys.getenv("PRIME_FRIEDMAN_DETAIL_PROFILE", "smoke")
replications <- switch(
  profile,
  smoke = 2L,
  main = 30L,
  stop("Unknown Friedman-detail profile: ", profile, call. = FALSE)
)
methods <- c("PRIME screen", "FIRST", "FOCI", "FORD")
directory <- file.path(
  "results/raw", paste0("friedman_detail_v052_", profile)
)
runs <- read.csv(
  file.path(directory, "friedman_detail_runs.csv"),
  stringsAsFactors = FALSE
)
phi <- read.csv(
  file.path(directory, "friedman_detail_phi.csv"),
  stringsAsFactors = FALSE
)
expected_runs <- replications * length(methods)
stopifnot(
  nrow(runs) == expected_runs,
  nrow(phi) == expected_runs * 100L,
  !anyDuplicated(runs[c("replication", "method")]),
  !anyDuplicated(phi[c("replication", "method", "variable")]),
  setequal(unique(runs$method), methods),
  all(is.finite(phi$truth)),
  all(is.finite(runs$half_l1[runs$feasible])),
  all(runs$error[runs$feasible] == "" | is.na(runs$error[runs$feasible]))
)
effect_sums <- aggregate(
  estimate ~ replication + method, phi, sum, na.rm = TRUE
)
effect_sums <- merge(
  effect_sums,
  runs[c("replication", "method", "feasible")],
  by = c("replication", "method")
)
stopifnot(max(abs(effect_sums$estimate[effect_sums$feasible] - 1)) < 1e-7)
if (profile == "main") {
  exact_counts <- aggregate(exact ~ method, runs, sum)
  exact_counts <- setNames(exact_counts$exact, exact_counts$method)
  feasible_counts <- aggregate(feasible ~ method, runs, sum)
  feasible_counts <- setNames(
    feasible_counts$feasible, feasible_counts$method
  )
  stopifnot(
    all(exact_counts[methods] == c(29L, 11L, 0L, 8L)),
    all(feasible_counts[methods] == c(30L, 30L, 29L, 30L))
  )
}
cat("Friedman coordinate-allocation validation: PASS\n")
print(aggregate(cbind(exact, active_inclusion, half_l1) ~ method, runs, mean,
                na.rm = TRUE))
