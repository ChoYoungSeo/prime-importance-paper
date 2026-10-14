summary_dir <- "results/summary"
required_files <- c(
  "selection_summary.csv",
  "selection_calibration_summary.csv",
  "selection_to_pme_summary.csv",
  "pme_estimator_calibration_summary.csv",
  "pme_correlated_pipeline_summary.csv",
  "pme_scaling_summary.csv",
  "pme_reference_allocations.csv",
  "observed_application_summary.csv",
  "observed_application_effects.csv",
  "twinning_v052_summary.csv",
  "shapley_pme_proxy_population.csv"
)

required_paths <- file.path(summary_dir, required_files)
missing_files <- required_paths[!file.exists(required_paths)]
if (length(missing_files)) {
  stop(
    "Missing tracked summary files: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

summaries <- lapply(required_paths, function(path) {
  data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(data) || anyDuplicated(names(data))) {
    stop("Invalid summary structure: ", path, call. = FALSE)
  }
  numeric_columns <- vapply(data, is.numeric, logical(1))
  if (any(numeric_columns) && any(is.infinite(as.matrix(data[numeric_columns])), na.rm = TRUE)) {
    stop("Infinite numeric value in: ", path, call. = FALSE)
  }
  data
})
names(summaries) <- required_files

required_methods <- list(
  selection_summary.csv = c("PRIME screen", "FIRST", "FOCI", "FORD"),
  selection_to_pme_summary.csv = c("PRIME screen", "FIRST", "FOCI", "FORD"),
  pme_estimator_calibration_summary.csv = c("PRIME", "PME-kNN", "NCMD-PME"),
  pme_correlated_pipeline_summary.csv = "PRIME",
  pme_scaling_summary.csv = "PRIME"
)

for (file in names(required_methods)) {
  methods <- unique(summaries[[file]]$method)
  missing_methods <- setdiff(required_methods[[file]], methods)
  if (length(missing_methods)) {
    stop(
      "Missing methods in ", file, ": ",
      paste(missing_methods, collapse = ", "),
      call. = FALSE
    )
  }
}

bounded_columns <- grep(
  "(^|_)(pme|effect|median|q05|q95|selection_frequency)$",
  names(summaries[["observed_application_effects.csv"]]),
  value = TRUE
)
observed_values <- unlist(
  summaries[["observed_application_effects.csv"]][bounded_columns],
  use.names = FALSE
)
if (any(observed_values < 0 | observed_values > 1, na.rm = TRUE)) {
  stop("Observed-data allocation summaries fall outside [0, 1].", call. = FALSE)
}

reference <- summaries[["pme_reference_allocations.csv"]]
reference_sums <- aggregate(pme ~ benchmark + rho, reference, sum)
if (max(abs(reference_sums$pme - 1)) > 1e-8) {
  stop("Reference PME allocations do not sum to one.", call. = FALSE)
}

required_provenance <- c(
  "results/provenance/validation_report.md",
  "results/provenance/theory_audit.md",
  "results/provenance/reference_audit.md",
  "results/provenance/r_session_info.txt"
)
if (any(!file.exists(required_provenance))) {
  stop("The tracked provenance package is incomplete.", call. = FALSE)
}

cat(
  "Repository snapshot validation passed.\n",
  "Summaries checked:", length(required_files), "\n",
  "Reference PME cells:", nrow(reference_sums), "\n"
)
