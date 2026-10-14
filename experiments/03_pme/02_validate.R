profile <- Sys.getenv("PRIME_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(benchmarks = c(branin = 8L), n = 1L, rho = 2L,
               replications = 1L),
  pilot = list(benchmarks = c(branin = 8L, bratley = 10L), n = 2L,
               rho = 2L, replications = 3L),
  main = list(benchmarks = c(branin = 8L, bratley = 10L), n = 3L,
              rho = 3L, replications = 30L),
  stop("Unknown profile: ", profile, call. = FALSE)
)
methods <- c(
  "PRIME", "Full local-linear PME", "PME-kNN",
  "NCMD-PME", "Oracle local-linear PME"
)
directory <- file.path("results/raw", paste0("pme_v052_", profile))
runs <- read.csv(file.path(directory, "pme_runs.csv"), stringsAsFactors = FALSE)
phi <- read.csv(file.path(directory, "pme_phi.csv"), stringsAsFactors = FALSE)

dataset_count <- length(configuration$benchmarks) * configuration$n *
  configuration$rho * configuration$replications
expected_runs <- dataset_count * length(methods)
expected_phi <- configuration$n * configuration$rho *
  configuration$replications * length(methods) *
  sum(configuration$benchmarks)

stopifnot(
  nrow(runs) == expected_runs,
  nrow(phi) == expected_phi,
  !anyDuplicated(runs[c(
    "benchmark", "n", "rho", "replication", "method"
  )]),
  !anyDuplicated(phi[c(
    "benchmark", "n", "rho", "replication", "method", "variable"
  )]),
  setequal(unique(runs$method), methods),
  all(is.na(runs$error) | runs$error == ""),
  all(is.finite(runs$half_l1)),
  all(runs$half_l1 >= 0 & runs$half_l1 <= 1),
  all(is.finite(phi$estimate)),
  all(is.finite(phi$truth))
)
if (profile == "main") {
  proposed <- "PRIME"
  current_origins <- c(
    paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  )
  stopifnot(
    all(runs$result_origin[runs$method == proposed] %in% current_origins),
    all(runs$result_origin[runs$method != proposed] ==
        "reused paired baseline"),
    all(phi$result_origin[phi$method == proposed] %in% current_origins),
    all(phi$result_origin[phi$method != proposed] ==
        "reused paired baseline")
  )
}

effect_sums <- aggregate(
  estimate ~ benchmark + n + rho + replication + method,
  phi,
  sum
)
truth_sums <- aggregate(
  truth ~ benchmark + n + rho + replication + method,
  phi,
  sum
)
stopifnot(
  max(abs(effect_sums$estimate - 1)) < 1e-7,
  max(abs(truth_sums$truth - 1)) < 1e-7
)

cat("Low-dimensional PME validation: PASS\n")
print(aggregate(
  cbind(half_l1, inactive_mass, proxy_mass, elapsed_sec) ~ method,
  runs,
  function(x) c(mean = mean(x), median = median(x))
))
