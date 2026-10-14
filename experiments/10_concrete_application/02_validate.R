profile <- Sys.getenv("PRIME_CONCRETE_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(bootstraps = 4L, cv_repeats = 1L, cv_folds = 3L),
  main = list(bootstraps = 200L, cv_repeats = 5L, cv_folds = 5L),
  stop("Unknown concrete profile: ", profile, call. = FALSE)
)
directory <- file.path("results/raw", paste0("concrete_v052_", profile))
runs <- read.csv(file.path(directory, "bootstrap_runs.csv"), stringsAsFactors = FALSE)
phi <- read.csv(file.path(directory, "bootstrap_phi.csv"), stringsAsFactors = FALSE)
cv <- read.csv(file.path(directory, "nested_cv.csv"), stringsAsFactors = FALSE)

expected_methods <- 4L
stopifnot(
  nrow(runs) == configuration$bootstraps * expected_methods,
  !anyDuplicated(runs[c("replication", "representation", "method")]),
  all(runs$feasible),
  all(runs$error == "" | is.na(runs$error)),
  nrow(phi) == configuration$bootstraps * expected_methods * 8L,
  !anyDuplicated(phi[c(
    "replication", "representation", "method", "variable"
  )]),
  all(is.finite(phi$estimate)),
  max(abs(
    aggregate(estimate ~ replication + representation + method, phi, sum)$estimate - 1
  )) < 1e-7,
  nrow(cv) == configuration$cv_repeats * configuration$cv_folds * 3L,
  !anyDuplicated(cv[c("repeat_id", "fold", "feature_set")]),
  all(is.finite(cv$rmse) & cv$rmse > 0)
)
cat("Concrete application validation: PASS\n")
print(aggregate(rmse ~ feature_set, cv, function(x) {
  c(mean = mean(x), median = median(x), q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75))
}))
