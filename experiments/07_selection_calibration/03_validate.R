profile <- Sys.getenv("PRIME_SELECTION_CALIBRATION_PROFILE", "smoke")
output_directory <- file.path(
  "results", "raw", paste0("selection_calibration_v052_", profile)
)
buffer_path <- file.path(output_directory, "buffer_ablation.csv")
level_path <- file.path(output_directory, "clipped_level.csv")
if (!file.exists(buffer_path) || !file.exists(level_path)) {
  stop("Run both selection-calibration experiments before validation.",
       call. = FALSE)
}

buffer <- read.csv(buffer_path, stringsAsFactors = FALSE)
level <- read.csv(level_path, stringsAsFactors = FALSE)
task_count <- switch(
  profile,
  smoke = 54L,
  pilot = 720L,
  main = 2250L,
  stop("Unknown selection-calibration profile: ", profile, call. = FALSE)
)
stopifnot(
  nrow(buffer) == 4L * task_count,
  nrow(level) == 2L * task_count,
  length(unique(buffer$task_id)) == task_count,
  length(unique(level$task_id)) == task_count
)

buffer_methods <- c(
  "buffered forward + buffered backward",
  "buffered forward + positive-part backward",
  "positive-part forward + buffered backward",
  "positive-part forward + positive-part backward"
)
level_methods <- c(
  "buffered paired contrast",
  "FIRST-style clipped level"
)
stopifnot(
  setequal(unique(buffer$method), buffer_methods),
  setequal(unique(level$method), level_methods)
)

paired_keys <- c("benchmark", "n", "design", "replication", "seed")
buffer_key <- buffer[
  buffer$method == "buffered forward + buffered backward",
  c(paired_keys, "selected")
]
level_key <- level[
  level$method == "buffered paired contrast",
  c(paired_keys, "selected")
]
comparison <- merge(
  buffer_key,
  level_key,
  by = paired_keys,
  suffixes = c("_buffer", "_level"),
  sort = FALSE
)
stopifnot(
  nrow(comparison) == task_count,
  identical(comparison$selected_buffer, comparison$selected_level)
)

summarize_method <- function(data, method, rule, forward_rule, backward_rule) {
  cell <- data[data$method == method, , drop = FALSE]
  data.frame(
    rule = rule,
    forward_rule = forward_rule,
    backward_rule = backward_rule,
    candidate_size = mean(cell$candidate_size),
    selected_size = mean(cell$selected_size),
    tpr = mean(cell$tpr),
    fpr = mean(cell$fpr),
    false_positives = mean(cell$false_positives),
    inclusion = mean(cell$inclusion),
    exact = mean(cell$exact),
    over12 = mean(cell$selected_size > 12L),
    stringsAsFactors = FALSE
  )
}

summary <- rbind(
  summarize_method(
    buffer, "buffered forward + buffered backward",
    "buffered / buffered", "buffered forward", "buffered backward"
  ),
  summarize_method(
    buffer, "buffered forward + positive-part backward",
    "buffered / positive-part",
    "buffered forward", "positive-part backward"
  ),
  summarize_method(
    buffer, "positive-part forward + buffered backward",
    "positive-part / buffered",
    "positive-part forward", "buffered backward"
  ),
  summarize_method(
    buffer, "positive-part forward + positive-part backward",
    "positive-part / positive-part",
    "positive-part forward", "positive-part backward"
  ),
  summarize_method(
    level, "FIRST-style clipped level",
    "FIRST-style clipped level",
    "FIRST-style clipped level", "FIRST-style clipped level"
  )
)

if (profile == "main") {
  dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)
  write.csv(
    summary,
    "results/summary/selection_calibration_summary.csv",
    row.names = FALSE
  )
}
print(summary, row.names = FALSE)
cat("Selection-calibration validation passed for", task_count, "paired tasks.\n")
