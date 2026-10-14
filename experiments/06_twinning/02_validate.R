profile <- Sys.getenv("PRIME_TWINNING_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(benchmarks = 1L, p = 1L, rho = 1L,
               replications = 1L, ratios = 2L),
  main = list(benchmarks = 3L, p = 2L, rho = 2L,
              replications = 5L, ratios = 4L),
  stop("Unknown Data-Twinning profile: ", profile, call. = FALSE)
)
path <- file.path(
  "results/raw", paste0("twinning_v052_", profile),
  "twinning_runs.csv"
)
runs <- read.csv(path, stringsAsFactors = FALSE)
expected <- with(
  configuration,
  benchmarks * p * rho * replications * ratios
)
stopifnot(
  nrow(runs) == expected,
  !anyDuplicated(runs[c(
    "benchmark", "p", "rho", "replication", "twin_r"
  )]),
  all(is.na(runs$error) | runs$error == ""),
  all(runs$query_count > 0L),
  all(runs$selected_size >= 0L),
  all(runs$active_recall >= 0 & runs$active_recall <= 1),
  all(runs$candidate_jaccard >= 0 & runs$candidate_jaccard <= 1),
  all(runs$selected_jaccard >= 0 & runs$selected_jaccard <= 1),
  all(is.finite(runs$speedup) & runs$speedup > 0)
)
all_query <- runs$twin_r == "all"
stopifnot(
  all(runs$candidate_jaccard[all_query] == 1),
  all(runs$selected_jaccard[all_query] == 1),
  all(abs(runs$speedup[all_query] - 1) < 1e-12)
)
cat("Data-Twinning validation: PASS\n")
print(aggregate(
  cbind(active_recall, selected_jaccard, speedup) ~ twin_r,
  runs,
  median
))
