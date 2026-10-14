source("experiments/01_selection/00_common.R")

profile <- Sys.getenv("PRIME_SELECTION_PROFILE", "smoke")
configuration <- selection_profile(profile)
methods <- c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
raw <- read.csv(selection_output_path(profile), stringsAsFactors = FALSE)

expected_tasks <- length(configuration$benchmarks) * length(configuration$n) *
  length(configuration$designs) * configuration$replications
expected_rows <- expected_tasks * length(methods)
stopifnot(
  nrow(raw) == expected_rows,
  !anyDuplicated(paste(
    raw$benchmark,
    raw$n,
    raw$design,
    raw$replication,
    raw$method,
    sep = "|"
  )),
  all(table(raw$task_id) == length(methods)),
  !any(nzchar(ifelse(is.na(raw$error), "", raw$error))),
  all(is.finite(raw$elapsed_sec) & raw$elapsed_sec >= 0)
)
if (profile == "main") {
  current_origins <- c(
    paste("prime", as.character(utils::packageVersion("prime")), "rerun")
  )
  stopifnot(
    all(raw$result_origin[raw$method == "PRIME screen"] %in% current_origins),
    all(raw$prime_version[raw$method == "PRIME screen"] ==
        as.character(utils::packageVersion("prime"))),
    all(raw$result_origin[raw$method != "PRIME screen"] %in%
        c("reused paired baseline", "paired Ackley-3 rerun")),
    all(raw$result_origin[
      raw$benchmark == "ackley3" & raw$method != "PRIME screen"
    ] == "paired Ackley-3 rerun")
  )
}

parse_indices <- function(value) {
  if (is.na(value) || !nzchar(trimws(value))) return(integer())
  sort(unique(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])))
}
selected <- lapply(raw$selected, parse_indices)
active <- lapply(raw$active, parse_indices)
false_positive <- mapply(
  function(selected_set, active_set) {
    length(setdiff(selected_set, active_set))
  },
  selected,
  active
)
false_negative <- mapply(
  function(selected_set, active_set) {
    length(setdiff(active_set, selected_set))
  },
  selected,
  active
)
stopifnot(
  all(raw$selected_size == lengths(selected)),
  all(raw$false_positives == false_positive),
  all(raw$false_negatives == false_negative),
  all(raw$inclusion == (false_negative == 0L)),
  all(raw$exact == mapply(identical, selected, active)),
  all(vapply(selected, function(index) {
    all(index >= 1L & index <= 100L)
  }, logical(1)))
)

task_map <- unique(raw[
  c("task_id", "benchmark", "n", "design", "replication", "seed")
])
task_map <- task_map[order(task_map$task_id), , drop = FALSE]
row.names(task_map) <- NULL
expected_map <- expand.grid(
  benchmark = configuration$benchmarks,
  n = configuration$n,
  design = configuration$designs,
  replication = seq_len(configuration$replications),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_map$task_id <- seq_len(nrow(expected_map))
expected_map$seed <- 202608080L +
  1000000L * match(expected_map$benchmark, configuration$benchmarks) +
  100000L * match(expected_map$design, configuration$designs) +
  1000L * match(expected_map$n, configuration$n) +
  expected_map$replication
expected_map <- expected_map[
  c("task_id", "benchmark", "n", "design", "replication", "seed")
]
stopifnot(identical(task_map, expected_map))

for (benchmark in configuration$benchmarks) {
  dimension <- selection_registry()[[benchmark]]$active_dimension
  for (design in configuration$designs) {
    expected_active <- selection_active_positions(100L, dimension, design)
    cell <- raw[raw$benchmark == benchmark & raw$design == design, ]
    stopifnot(all(cell$active == paste(expected_active, collapse = ",")))
    covariance <- selection_covariance(100L, expected_active, design, 0.5)
    if (design == "spaced_ar1") {
      block <- covariance[expected_active, expected_active, drop = FALSE]
      diag(block) <- 0
      stopifnot(max(abs(block)) < 1e-5)
    }
    if (design == "balanced_proxy") {
      proxies <- max(expected_active) + seq_along(expected_active)
      stopifnot(
        all(abs(diag(covariance[expected_active, proxies]) - 0.5) < 1e-12),
        all(cell$proxies == paste(proxies, collapse = ","))
      )
    }
  }
}

cat("Selection validation: PASS\n")
cat("  tasks:", expected_tasks, "\n")
cat("  rows:", expected_rows, "\n")
