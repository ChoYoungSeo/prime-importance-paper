table_directory <- "manuscript/tables"
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

escape_tex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&#_$])", "\\\\\\1", x)
  x
}

write_lines <- function(lines, path) {
  writeLines(lines, path, useBytes = TRUE)
}

calibration <- read.csv(
  "results/summary/selection_calibration_summary.csv",
  stringsAsFactors = FALSE
)
calibration_order <- c(
  "buffered / buffered",
  "buffered / positive-part",
  "positive-part / buffered",
  "positive-part / positive-part",
  "FIRST-style clipped level"
)
calibration <- calibration[
  match(calibration_order, calibration$rule),
  ,
  drop = FALSE
]
calibration_forward <- c(
  "Buffered", "Buffered", "Raw", "Raw", "Clipped level"
)
calibration_backward <- c(
  "Buffered", "Raw", "Buffered", "Raw", "---"
)
calibration_rows <- vapply(seq_len(nrow(calibration)), function(index) {
  row <- calibration[index, ]
  sprintf(
    "%s & %s & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
    calibration_forward[index],
    calibration_backward[index],
    row$selected_size,
    row$inclusion,
    row$exact,
    row$false_positives,
    row$over12
  )
}, character(1))
write_lines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    "\\caption{Matched selection-rule ablation over 2,250 simulated datasets.}",
    "\\label{tab:selection-calibration-ablation}",
    "\\begin{tabular}{llccccc}",
    "\\toprule",
    "Forward & Backward & Mean $|\\widehat S|$ & Inclusion & Exact & Mean FP & $\\Pr(|\\widehat S|>12)$ \\\\",
    "\\midrule",
    calibration_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "selection_calibration_ablation.tex")
)

low <- read.csv(
  "results/raw/pme_v052_main/pme_runs.csv",
  stringsAsFactors = FALSE
)
low$method <- factor(
  low$method,
  levels = c(
    "PRIME",
    "Full local-linear PME",
    "PME-kNN",
    "NCMD-PME",
    "Oracle local-linear PME"
  ),
  labels = c(
    "PRIME", "Full-coordinate local-linear PME", "PME-kNN", "NCMD-PME",
    "Oracle-coordinate LL-PME"
  )
)
low_summary <- do.call(
  rbind,
  lapply(split(low, low$method, drop = TRUE), function(cell) {
    data.frame(
      method = as.character(cell$method[1]),
      mean_error = mean(cell$half_l1),
      median_error = median(cell$half_l1),
      high_rho_error = mean(cell$half_l1[cell$rho == 0.9]),
      mean_time = mean(cell$elapsed_sec),
      stringsAsFactors = FALSE
    )
  })
)
low_summary <- low_summary[order(low_summary$method), ]
low_summary <- low_summary[
  match(levels(low$method), low_summary$method),
]
low_rows <- vapply(seq_len(nrow(low_summary)), function(index) {
  row <- low_summary[index, ]
  sprintf(
    "%s & %.3f & %.3f & %.3f & %.2f \\\\",
    escape_tex(row$method),
    row$mean_error,
    row$median_error,
    row$high_rho_error,
    row$mean_time
  )
}, character(1))
write_lines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{PME estimator performance pooled over two response functions, three sample sizes, three dependence levels, and 30 paired replications per cell.  Error is half-$L_1$ allocation loss.}",
    "\\label{tab:pme-estimator-diagnostics}",
    "\\small",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Method & Mean error & Median error & Mean error, $\\rho=0.9$ & Mean time (s) \\\\",
    "\\midrule",
    low_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "pme_estimator_diagnostics.tex")
)

scaling <- read.csv(
  "results/raw/scaling_v052_main/scaling_runs.csv",
  stringsAsFactors = FALSE
)
scaling <- scaling[scaling$p == max(scaling$p), ]
method_order <- c(
  "PRIME",
  "FIRST screen + PME-kNN",
  "NNVS screen + PME-kNN",
  "LASSO screen + PME-kNN",
  "FOCI screen + PME-kNN",
  "FORD screen + PME-kNN"
)
method_label <- c(
  "PRIME",
  "FIRST + PME-kNN",
  "NNVS + PME-kNN",
  "LASSO + PME-kNN",
  "FOCI + PME-kNN",
  "FORD + PME-kNN"
)
benchmark_order <- c("otl", "piston", "borehole")
benchmark_label <- c("OTL", "Piston", "Borehole")

scaling_rows <- character()
for (benchmark_index in seq_along(benchmark_order)) {
  benchmark <- benchmark_order[benchmark_index]
  for (method_index in seq_along(method_order)) {
    cell <- scaling[
      scaling$benchmark == benchmark &
        scaling$method == method_order[method_index],
    ]
    scaling_rows <- c(
      scaling_rows,
      sprintf(
        "%s & %s & %.3f & %.2f & %.1f & %d/%d \\\\",
        if (method_index == 1L) benchmark_label[benchmark_index] else "",
        method_label[method_index],
        median(cell$half_l1, na.rm = TRUE),
        median(cell$elapsed_sec),
        median(cell$selected_size),
        sum(cell$feasible),
        nrow(cell)
      )
    )
  }
  if (benchmark_index < length(benchmark_order)) {
    scaling_rows <- c(scaling_rows, "\\addlinespace")
  }
}
write_lines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Performance at $n=p=500$.  Loss is the median half-$L_1$ error among completed allocations; total runtime and selected dimension are medians over all 30 paired runs.}",
    "\\label{tab:pme-scaling-pmax}",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{7pt}",
    "\\begin{tabular}{llrrrr}",
    "\\toprule",
    "Response & Method & Median loss & Total time (s) & Selected $d$ & Completed \\\\",
    "\\midrule",
    scaling_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "pme_scaling_pmax.tex")
)

correlated <- read.csv(
  "results/raw/correlated_pipelines_v052_main/correlated_pipeline_runs.csv",
  stringsAsFactors = FALSE
)
correlated_method_order <- c(
  "PRIME",
  "FIRST screen + PME-kNN",
  "NNVS screen + PME-kNN",
  "LASSO screen + PME-kNN",
  "FOCI screen + PME-kNN",
  "FORD screen + PME-kNN"
)
correlated_method_label <- c(
  "PRIME", "FIRST + PME-kNN", "NNVS + PME-kNN",
  "LASSO + PME-kNN", "FOCI + PME-kNN", "FORD + PME-kNN"
)
correlated_rows <- vapply(seq_along(correlated_method_order), function(index) {
  cell <- correlated[correlated$method == correlated_method_order[index], ]
  sprintf(
    "%s & %.3f & %.3f & %.4f & %d/%d \\\\",
    correlated_method_label[index],
    mean(cell$half_l1, na.rm = TRUE),
    mean(cell$active_inclusion),
    mean(cell$proxy_mass, na.rm = TRUE),
    sum(cell$feasible),
    nrow(cell)
  )
}, character(1))
write_lines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Correlated-proxy performance pooled over Branin and Bratley, three sample sizes, and 20 paired replications per cell.}",
    "\\label{tab:correlated-pipeline}",
    "\\small",
    "\\begin{tabular}{lrrrr}",
    "\\toprule",
    "Method & Mean loss & Active-set coverage & Proxy allocation & Completed \\\\",
    "\\midrule",
    correlated_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "pme_correlated_pipeline.tex")
)

twinning <- read.csv(
  "results/raw/twinning_v052_main/twinning_runs.csv",
  stringsAsFactors = FALSE
)
twinning <- twinning[twinning$twin_r != "all", , drop = FALSE]
twinning_ratios <- c("2", "3", "5")
twinning_rows <- vapply(twinning_ratios, function(ratio) {
  cell <- twinning[twinning$twin_r == ratio, ]
  speedup_interval <- stats::quantile(
    cell$speedup, c(0.25, 0.75), na.rm = TRUE
  )
  sprintf(
    "%s & %.0f & %.3f & %.3f & %.2f [%.2f, %.2f] \\\\",
    ratio,
    stats::median(cell$query_count),
    stats::median(cell$selected_jaccard),
    stats::median(cell$active_recall),
    stats::median(cell$speedup),
    speedup_interval[1],
    speedup_interval[2]
  )
}, character(1))
write_lines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Representative-query approximation at $n=500$.  For split ratio $n/q$, the $q$ twinned rows serve as query locations while all observations remain available as neighbors.  Speedup is the median [interquartile range] over 60 paired datasets per ratio.}",
    "\\label{tab:twinning-diagnostic}",
    "\\small",
    "\\begin{tabular}{rrrrr}",
    "\\toprule",
    "Split ratio $n/q$ & Queries $q$ & Selected-set Jaccard & Active recall & Speedup \\\\",
    "\\midrule",
    twinning_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "twinning_diagnostic.tex")
)
