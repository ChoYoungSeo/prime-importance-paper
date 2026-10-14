matched <- read.csv(
  "results/raw/ridge_matched/matched_runs.csv", stringsAsFactors = FALSE
)
validation <- jsonlite::fromJSON("results/raw/ridge_matched/validation.json")
stopifnot(
  isTRUE(validation$complete), validation$validated_rows == 870L,
  validation$unique_datasets == 290L, nrow(matched) == 870L,
  length(unique(matched$dataset_id)) == 290L,
  !anyDuplicated(matched$task_id),
  all(matched$status == "complete"),
  all(matched$forward_buffer), all(matched$backward_buffer),
  all(matched$seed_rule %in% c("fixed", "bic", "empty")),
  all(table(matched$dataset_id, matched$seed_rule) == 1L)
)
general <- matched[matched$family == "benchmark", , drop = FALSE]
suppressor <- matched[matched$family == "suppressor", , drop = FALSE]
rules <- c("fixed", "bic", "empty")
labels <- c(fixed = "Fixed ridge", bic = "BIC-gated", empty = "Empty")

general <- general[general$seed_rule %in% rules, , drop = FALSE]
suppressor <- suppressor[suppressor$seed_rule %in% rules, , drop = FALSE]
general_summary <- aggregate(
  inclusion ~ benchmark + seed_rule, general, mean
)
suppressor_summary <- aggregate(
  inclusion ~ n + noise_sd + seed_rule, suppressor, mean
)

extract_row <- function(data, condition, label) {
  cell <- data[condition, , drop = FALSE]
  value <- setNames(cell$inclusion, cell$seed_rule)
  data.frame(
    scenario = label,
    fixed = value[["fixed"]],
    bic = value[["bic"]],
    empty = value[["empty"]],
    stringsAsFactors = FALSE
  )
}

table <- rbind(
  extract_row(
    general_summary,
    general_summary$benchmark == "ackley3",
    "Ackley--3, $n=p=100$"
  ),
  extract_row(
    general_summary,
    general_summary$benchmark == "friedman",
    "Friedman, $n=p=100$"
  ),
  extract_row(
    suppressor_summary,
    suppressor_summary$n == 100 & suppressor_summary$noise_sd == 0.1,
    "Suppressor, $n=100$, $\\sigma=0.1$"
  ),
  extract_row(
    suppressor_summary,
    suppressor_summary$n == 500 & suppressor_summary$noise_sd == 0.5,
    "Suppressor, $n=500$, $\\sigma=0.5$"
  )
)

format_value <- function(value) sprintf("%.3f", value)
lines <- c(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Active-set coverage under alternative initializations, with the same calibrated forward and backward rules in every column. The benchmark rows pool 10 replications over three predictor laws; the suppressor rows use 50 replications.}",
  "\\label{tab:ridge-seed-diagnostic}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  paste(c("Scenario", unname(labels)), collapse = " & "),
  "\\\\",
  "\\midrule",
  vapply(seq_len(nrow(table)), function(index) {
    paste0(
      table$scenario[index], " & ",
      format_value(table$fixed[index]), " & ",
      format_value(table$bic[index]), " & ",
      format_value(table$empty[index]), " \\\\"
    )
  }, character(1)),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
dir.create("manuscript/tables", recursive = TRUE, showWarnings = FALSE)
writeLines(lines, "manuscript/tables/ridge_seed_diagnostic.tex")
dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)
write.csv(table, "results/summary/ridge_seed_diagnostic.csv", row.names = FALSE)
