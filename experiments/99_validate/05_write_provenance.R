classify_file <- function(path) {
  if (grepl("^software/prime/", path)) return("software source")
  if (grepl("^experiments/", path)) return("experiment source")
  if (grepl(
    "^manuscript/sections/|^manuscript/.*\\.tex$|^manuscript/references\\.bib$",
    path
  )) {
    return("manuscript source")
  }
  if (grepl("^data/", path)) return("observed-data snapshot")
  if (grepl("^results/revisions/", path)) return("verified release input")
  if (grepl("^results/raw/", path)) return("raw result")
  if (grepl("^results/summary/", path)) return("result summary")
  if (grepl("^results/provenance/inherited/", path)) {
    return("inherited data snapshot")
  }
  if (grepl("^results/provenance/paired_baselines/", path)) {
    return("paired competitor baseline")
  }
  if (grepl("^results/provenance/", path)) return("validation record")
  if (grepl(
    "^manuscript/figures/|^manuscript/tables/",
    path
  )) return("paper artifact")
  if (grepl("^manuscript/(main|supplement)\\.pdf$", path)) {
    return("compiled manuscript")
  }
  "project specification"
}

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  if (length(output) != 1L) {
    stop("Could not hash ", path, call. = FALSE)
  }
  sub("[[:space:]].*$", "", output)
}

source_files <- c(
  "README.md", "VALIDATION.md", "CITATION.cff", "METHOD_SPEC.md", "Makefile",
  ".github/workflows/R-CMD-check.yaml", ".gitignore",
  "config/experiment_design.md",
  list.files("software/prime", recursive = TRUE, full.names = TRUE),
  list.files("experiments", pattern = "\\.([Rr]|py)$", recursive = TRUE,
             full.names = TRUE),
  list.files("manuscript/sections", pattern = "\\.tex$",
             recursive = TRUE,
             full.names = TRUE),
  "manuscript/main.tex",
  "manuscript/supplement.tex",
  "manuscript/references.bib",
  list.files("data", recursive = TRUE, full.names = TRUE)
)
release_manifest <- jsonlite::fromJSON("results/current_release.json")
stopifnot(release_manifest$schema_version == 2L,
          release_manifest$status == "validated_public_release")
result_files <- c(
  "results/current_release.json",
  list.files(release_manifest$data_root, recursive = TRUE, full.names = TRUE),
  "results/raw/selection_v052_runs.csv",
  list.files("results/raw", pattern = "\\.csv$", recursive = TRUE,
             full.names = TRUE),
  list.files("results/summary", pattern = "\\.csv$", recursive = TRUE,
             full.names = TRUE),
  list.files("results/provenance/historical", recursive = TRUE,
             full.names = TRUE),
  list.files("results/provenance/release_lineage", recursive = TRUE,
             full.names = TRUE),
  list.files("results/provenance/inherited", recursive = TRUE,
             full.names = TRUE),
  list.files("results/provenance/paired_baselines", recursive = TRUE,
             full.names = TRUE),
  list.files("results/provenance/reference_cache", recursive = TRUE,
             full.names = TRUE),
  setdiff(
    list.files("results/provenance", recursive = FALSE, full.names = TRUE),
    c(
      "results/provenance/file_manifest.csv",
      "results/provenance/r_session_info.txt"
    )
  ),
  "results/provenance/pme_knn_subset_value_parity.csv",
  "results/provenance/selection_implementation_parity.csv",
  "results/provenance/pme_reference_objects.rds",
  "results/provenance/R_CMD_check.log",
  "results/provenance/validation_report.md"
)
paper_files <- c(
  list.files("manuscript/figures", recursive = TRUE,
             full.names = TRUE),
  list.files("manuscript/tables", recursive = TRUE,
             full.names = TRUE),
  c("manuscript/main.pdf",
    "manuscript/supplement.pdf")
)
files <- sort(unique(c(source_files, result_files, paper_files)))
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[!grepl("\\.(o|so)$", files)]

manifest <- data.frame(
  path = files,
  role = vapply(files, classify_file, character(1)),
  bytes = unname(file.info(files)$size),
  sha256 = vapply(files, sha256, character(1)),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  "results/provenance/file_manifest.csv",
  row.names = FALSE
)

session <- capture.output(sessionInfo())
writeLines(session, "results/provenance/r_session_info.txt")

cat("Wrote provenance for", nrow(manifest), "files.\n")
