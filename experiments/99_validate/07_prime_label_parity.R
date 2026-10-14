source_files <- c(
  Sys.glob("results/raw/*v030*/*.csv"),
  Sys.glob("results/raw/*v030*.csv")
)
source_files <- source_files[!grepl("checkpoints", source_files, fixed = TRUE)]
target_files <- sub("v030", "v040", source_files, fixed = TRUE)

if (!length(source_files) || any(!file.exists(target_files))) {
  stop("The PRIME label-migration file set is incomplete.", call. = FALSE)
}

prime_label <- function(path) {
  if (grepl(
    "selection_to_pme|selection_v030|friedman_detail",
    path
  )) {
    "PRIME screen"
  } else {
    "PRIME"
  }
}

normalize <- function(data, path) {
  if ("method" %in% names(data)) {
    data$method[data$method == "EVSA"] <- prime_label(path)
  }
  data$result_origin <- NULL
  data$evsa_version <- NULL
  data$prime_version <- NULL
  data
}

for (index in seq_along(source_files)) {
  before <- normalize(
    read.csv(source_files[index], stringsAsFactors = FALSE),
    source_files[index]
  )
  after <- normalize(
    read.csv(target_files[index], stringsAsFactors = FALSE),
    target_files[index]
  )
  comparison <- all.equal(
    before,
    after,
    tolerance = 0,
    check.attributes = TRUE
  )
  if (!isTRUE(comparison)) {
    stop(
      "PRIME numerical parity failed for ", source_files[index], ": ",
      paste(comparison, collapse = "; "),
      call. = FALSE
    )
  }
}

report <- c(
  "PRIME label migration validation",
  paste0("Files compared: ", length(source_files)),
  "Selection-only labels: EVSA -> PRIME screen",
  "Allocation and complete-pipeline labels: EVSA -> PRIME",
  "Numerical and structural parity after label normalization: exact"
)
writeLines(report, "results/provenance/prime_label_migration_validation.txt")
cat(paste(report, collapse = "\n"), "\n")
