old <- read.csv(
  "results/raw/selection_v020_runs.csv",
  stringsAsFactors = FALSE
)
current <- read.csv(
  "results/raw/selection_v040_runs.csv",
  stringsAsFactors = FALSE
)
old$method[old$method == "FIRST-LLR"] <- "PRIME screen"

keys <- c("benchmark", "n", "design", "replication", "method", "seed")
stopifnot(
  !anyDuplicated(old[keys]),
  !anyDuplicated(current[keys]),
  setequal(
    do.call(paste, c(old[keys], sep = "|")),
    do.call(paste, c(current[keys], sep = "|"))
  )
)

old_prime <- old[old$method == "PRIME screen", ]
current_prime <- current[current$method == "PRIME screen", ]
paired <- merge(
  old_prime,
  current_prime,
  by = keys,
  suffixes = c("_cutoff", "_positive"),
  sort = FALSE
)
stopifnot(nrow(paired) == nrow(old_prime))

old_selected <- ifelse(
  is.na(paired$selected_cutoff), "", paired$selected_cutoff
)
current_selected <- ifelse(
  is.na(paired$selected_positive), "", paired$selected_positive
)
paired$selected_set_changed <- old_selected != current_selected

summary <- aggregate(
  cbind(
    inclusion_cutoff,
    inclusion_positive,
    exact_cutoff,
    exact_positive,
    tpr_cutoff,
    tpr_positive,
    false_positives_cutoff,
    false_positives_positive,
    selected_set_changed
  ) ~ benchmark + n + design,
  paired,
  mean
)
write.csv(
  summary,
  "results/provenance/strict_positive_gain_main_comparison.csv",
  row.names = FALSE
)

old_competitors <- old[old$method != "PRIME screen", ]
current_competitors <- current[current$method != "PRIME screen", ]
competitors <- merge(
  old_competitors,
  current_competitors,
  by = keys,
  suffixes = c("_old", "_current"),
  sort = FALSE
)
stopifnot(
  nrow(competitors) == nrow(old_competitors),
  identical(
    ifelse(is.na(competitors$selected_old), "", competitors$selected_old),
    ifelse(
      is.na(competitors$selected_current), "", competitors$selected_current
    )
  ),
  isTRUE(all.equal(
    competitors$tpr_old,
    competitors$tpr_current,
    check.attributes = FALSE
  )),
  identical(competitors$inclusion_old, competitors$inclusion_current),
  identical(competitors$exact_old, competitors$exact_current)
)

cat("Forward-rule comparison: PASS\n")
cat("  paired PRIME-screen tasks:", nrow(paired), "\n")
cat("  changed final sets:", sum(paired$selected_set_changed), "\n")
cat(
  "  exact recovery, cutoff/positive:",
  mean(paired$exact_cutoff), "/", mean(paired$exact_positive), "\n"
)
cat(
  "  active inclusion, cutoff/positive:",
  mean(paired$inclusion_cutoff), "/", mean(paired$inclusion_positive), "\n"
)
