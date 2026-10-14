profile <- Sys.getenv("PRIME_SELECTION_LINK_PROFILE", "smoke")
replications <- switch(
  profile,
  smoke = 1L,
  pilot = 3L,
  main = 30L,
  stop("Unknown profile: ", profile, call. = FALSE)
)
directory <- file.path(
  "results/raw", paste0("selection_to_pme_v052_", profile)
)
runs <- read.csv(
  file.path(directory, "selection_to_pme.csv"),
  stringsAsFactors = FALSE
)
phi <- read.csv(
  file.path(directory, "selection_to_pme_phi.csv"),
  stringsAsFactors = FALSE
)
methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
expected_runs <- 3L * 2L * replications * length(methods)
expected_phi <- expected_runs * 100L

stopifnot(
  nrow(runs) == expected_runs,
  nrow(phi) == expected_phi,
  !anyDuplicated(runs$link_id),
  setequal(runs$link_id, seq_len(expected_runs)),
  all(table(phi$link_id) == 100L),
  !anyDuplicated(runs[c("benchmark", "design", "replication", "method")]),
  !anyDuplicated(phi[c(
    "benchmark", "design", "replication", "method", "variable"
  )]),
  setequal(unique(runs$method), methods),
  all(is.finite(runs$half_l1[runs$feasible])),
  all(runs$half_l1[runs$feasible] >= 0 &
      runs$half_l1[runs$feasible] <= 1),
  all(is.finite(phi$truth)),
  all(is.na(runs$half_l1[!runs$feasible])),
  all(!runs$feasible[runs$selected_size == 0L]),
  all(runs$feasible[
    runs$method == "PRIME screen" &
      runs$selected_size > 0L & runs$selected_size <= 12L
  ]),
  all(!runs$feasible[
    runs$method == "PRIME screen" & runs$selected_size > 12L
  ])
)
ackley_truth <- unique(phi[
  phi$benchmark == "ackley3", c("variable", "truth")
])
ackley_truth <- ackley_truth[order(ackley_truth$variable), ]
references <- readRDS("results/provenance/pme_reference_objects.rds")
expected_ackley <- numeric(100L)
expected_ackley[1:3] <- references$ackley3_rho0.0$pme
ackley_reference <- readRDS(file.path(
  "results/provenance/reference_cache",
  references$ackley3_rho0.0$source
))
stopifnot(
  nrow(ackley_truth) == 100L,
  max(abs(ackley_truth$truth - expected_ackley)) < 1e-10,
  max(abs(expected_ackley[1:3] - 1 / 3)) < 1e-12,
  ackley_reference$n_ref == 500000L,
  is.finite(ackley_reference$max_half_pme_difference),
  ackley_reference$max_half_pme_difference < 0.01
)
if (profile == "main") {
  current_origins <- c(
    paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  )
  stopifnot(
    all(runs$result_origin[runs$method == "PRIME screen"] %in%
        current_origins),
    all(runs$result_origin[runs$method != "PRIME screen"] %in%
        c("reused paired baseline", "paired Ackley-3 rerun")),
    all(phi$result_origin[phi$method == "PRIME screen"] %in%
        current_origins),
    all(phi$result_origin[phi$method != "PRIME screen"] %in%
        c("reused paired baseline", "paired Ackley-3 rerun")),
    all(runs$result_origin[
      runs$benchmark == "ackley3" & runs$method != "PRIME screen"
    ] == "paired Ackley-3 rerun"),
    all(phi$result_origin[
      phi$benchmark == "ackley3" & phi$method != "PRIME screen"
    ] == "paired Ackley-3 rerun")
  )
}

effect_sums <- aggregate(
  estimate ~ benchmark + design + replication + method,
  phi,
  function(x) sum(x, na.rm = TRUE)
)
effect_sums <- merge(
  effect_sums,
  runs[c(
    "benchmark", "design", "replication", "method",
    "selected_size", "feasible"
  )],
  by = c("benchmark", "design", "replication", "method"),
  all.x = TRUE
)
nonempty <- effect_sums$feasible & effect_sums$selected_size > 0L
empty <- effect_sums$feasible & effect_sums$selected_size == 0L
stopifnot(
  max(abs(effect_sums$estimate[nonempty] - 1)) < 1e-7
)
if (any(empty)) {
  stopifnot(max(abs(effect_sums$estimate[empty])) < 1e-12)
}
if (any(runs$feasible)) {
  stopifnot(all(is.na(runs$error[runs$feasible]) |
                runs$error[runs$feasible] == ""))
}

cat("Selector-to-PME validation: PASS\n")
print(table(runs$method, runs$feasible, useNA = "ifany"))
print(aggregate(half_l1 ~ method, runs, mean, na.rm = TRUE))
print(aggregate(cbind(active_inclusion, feasible) ~ method, runs, mean))
