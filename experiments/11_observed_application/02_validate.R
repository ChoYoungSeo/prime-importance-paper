source("experiments/09_observed_data/00_common.R")

profile <- Sys.getenv("PRIME_OBSERVED_APPLICATION_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(bootstraps = 4L, cv_repeats = 1L, cv_folds = 3L),
  main = list(bootstraps = 200L, cv_repeats = 5L, cv_folds = 5L),
  stop("Unknown observed-application profile: ", profile, call. = FALSE)
)
datasets <- c("airfoil", "auto_mpg")
methods <- c("PRIME", "Full LL-PME", "PME-kNN")
directory <- file.path(
  "results/raw", paste0("observed_application_v052_", profile)
)
runs <- read.csv(file.path(directory, "bootstrap_runs.csv"), stringsAsFactors = FALSE)
phi <- read.csv(file.path(directory, "bootstrap_phi.csv"), stringsAsFactors = FALSE)
cv <- read.csv(file.path(directory, "nested_cv.csv"), stringsAsFactors = FALSE)

expected_runs <- length(datasets) * configuration$bootstraps * length(methods)
stopifnot(
  nrow(runs) == expected_runs,
  !anyDuplicated(runs[c("dataset", "replication", "method")]),
  setequal(unique(runs$dataset), datasets),
  setequal(unique(runs$method), methods),
  all(runs$feasible),
  all(runs$error == "" | is.na(runs$error)),
  !anyDuplicated(phi[c("dataset", "replication", "method", "variable")]),
  all(is.finite(phi$estimate))
)
coordinate_count <- aggregate(
  variable ~ dataset + replication + method, phi, length
)
coordinate_count <- merge(
  coordinate_count,
  unique(runs[c("dataset", "replication", "method", "p")]),
  by = c("dataset", "replication", "method")
)
stopifnot(all(coordinate_count$variable == coordinate_count$p))
effect_sum <- aggregate(
  estimate ~ dataset + replication + method, phi, sum
)
stopifnot(max(abs(effect_sum$estimate - 1)) < 1e-7)

expected_cv <- length(datasets) * configuration$cv_repeats *
  configuration$cv_folds * 2L
stopifnot(
  nrow(cv) == expected_cv,
  !anyDuplicated(cv[c("dataset", "repeat_id", "fold", "feature_set")]),
  setequal(
    unique(cv$feature_set), c("All coordinates", "Selected coordinates")
  ),
  all(is.finite(cv$rmse) & cv$rmse > 0),
  all(cv$selected_size > 0)
)

cat("Observed application validation: PASS\n")
print(aggregate(rmse ~ dataset + feature_set, cv, function(x) {
  c(mean = mean(x), median = median(x), q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75))
}))
