library(ggplot2)

figure_directory <- "manuscript/figures"
table_directory <- "manuscript/tables"
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

paper_theme <- function() {
  theme_classic(base_size = 10) +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "black", size = 8.5),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      plot.margin = margin(6, 8, 5, 6)
    )
}

observed_phi <- read.csv(
  "results/raw/observed_application_v052_main/bootstrap_phi.csv",
  stringsAsFactors = FALSE
)
observed_runs <- read.csv(
  "results/raw/observed_application_v052_main/bootstrap_runs.csv",
  stringsAsFactors = FALSE
)
observed_cv <- read.csv(
  "results/raw/observed_application_v052_main/nested_cv.csv",
  stringsAsFactors = FALSE
)
concrete_phi <- read.csv(
  "results/raw/concrete_v052_main/bootstrap_phi.csv",
  stringsAsFactors = FALSE
)
concrete_runs <- read.csv(
  "results/raw/concrete_v052_main/bootstrap_runs.csv",
  stringsAsFactors = FALSE
)
concrete_cv <- read.csv(
  "results/raw/concrete_v052_main/nested_cv.csv",
  stringsAsFactors = FALSE
)

concrete_phi <- concrete_phi[
  concrete_phi$representation == "engineering", , drop = FALSE
]
concrete_runs <- concrete_runs[
  concrete_runs$representation == "engineering", , drop = FALSE
]

allocation_distance <- function(data, ids) {
  data <- data[data$method %in% c("PRIME", "Full LL-PME"), ]
  wide <- reshape(
    data[c(ids, "variable", "method", "estimate")],
    idvar = c(ids, "variable"), timevar = "method", direction = "wide"
  )
  wide$absolute_difference <- abs(
    wide[["estimate.PRIME"]] - wide[["estimate.Full LL-PME"]]
  )
  formula <- stats::as.formula(paste("absolute_difference ~", paste(ids, collapse = "+")))
  out <- aggregate(formula, wide, sum)
  names(out)[ncol(out)] <- "total_variation"
  out$total_variation <- out$total_variation / 2
  out
}

observed_tv <- allocation_distance(
  observed_phi, c("dataset", "replication")
)
concrete_tv <- allocation_distance(concrete_phi, "replication")
concrete_tv$dataset <- "concrete"
tv <- rbind(
  observed_tv[c("dataset", "replication", "total_variation")],
  concrete_tv[c("dataset", "replication", "total_variation")]
)

cv_ratio <- function(data, ids) {
  data <- data[data$feature_set %in% c(
    "All coordinates", "Selected coordinates"
  ), ]
  wide <- reshape(
    data[c(ids, "feature_set", "rmse", "selected_size")],
    idvar = ids, timevar = "feature_set", direction = "wide"
  )
  wide$rmse_ratio <-
    wide[["rmse.Selected coordinates"]] / wide[["rmse.All coordinates"]]
  wide$selected_dimension <- wide[["selected_size.Selected coordinates"]]
  wide$full_rmse <- wide[["rmse.All coordinates"]]
  wide$selected_rmse <- wide[["rmse.Selected coordinates"]]
  wide
}

observed_ratio <- cv_ratio(
  observed_cv, c("dataset", "repeat_id", "fold")
)
concrete_ratio <- cv_ratio(concrete_cv, c("repeat_id", "fold"))
concrete_ratio$dataset <- "concrete"
ratio <- rbind(
  observed_ratio[c(
    "dataset", "repeat_id", "fold", "rmse_ratio", "selected_dimension",
    "full_rmse", "selected_rmse"
  )],
  concrete_ratio[c(
    "dataset", "repeat_id", "fold", "rmse_ratio", "selected_dimension",
    "full_rmse", "selected_rmse"
  )]
)

dataset_order <- c("concrete", "airfoil", "auto_mpg")
dataset_labels <- c(
  concrete = "Concrete", airfoil = "Airfoil", auto_mpg = "Auto MPG"
)
dataset_palette <- c(
  concrete = "#0072B2", airfoil = "#009E73", auto_mpg = "#D55E00"
)
tv$dataset <- factor(tv$dataset, levels = dataset_order)
ratio$dataset <- factor(ratio$dataset, levels = dataset_order)

observed_allocation <- observed_phi[
  observed_phi$method %in% c("PRIME", "Full LL-PME"),
  c("dataset", "variable_name", "method", "estimate")
]
concrete_allocation <- concrete_phi[
  concrete_phi$method %in% c("PRIME", "Full LL-PME"),
  c("variable_name", "method", "estimate")
]
concrete_allocation$dataset <- "concrete"
allocation_source <- rbind(
  observed_allocation,
  concrete_allocation[c("dataset", "variable_name", "method", "estimate")]
)
allocation_median <- aggregate(
  estimate ~ dataset + variable_name + method,
  allocation_source,
  median
)
allocation_points <- reshape(
  allocation_median,
  idvar = c("dataset", "variable_name"),
  timevar = "method",
  direction = "wide"
)
allocation_points$full_pme <- allocation_points[["estimate.Full LL-PME"]]
allocation_points$screened_pme <-
  allocation_points[["estimate.PRIME"]]
allocation_points$dataset <- factor(
  allocation_points$dataset, levels = dataset_order
)

p_allocation <- ggplot(
  allocation_points,
  aes(full_pme, screened_pme, color = dataset)
) +
  geom_abline(
    slope = 1, intercept = 0, color = "grey45", linewidth = 0.45
  ) +
  geom_point(size = 2.1, alpha = 0.9) +
  facet_wrap(
    ~dataset, nrow = 1, labeller = as_labeller(dataset_labels)
  ) +
  scale_color_manual(values = dataset_palette) +
  scale_x_sqrt(
    limits = c(0, 0.7), breaks = c(0, 0.1, 0.25, 0.5),
    expand = expansion(mult = c(0.035, 0.02))
  ) +
  scale_y_sqrt(
    limits = c(0, 0.7), breaks = c(0, 0.05, 0.1, 0.25, 0.5, 0.7),
    expand = expansion(mult = c(0.035, 0.02))
  ) +
  coord_equal() +
  labs(
    title = "A  Variable-level allocation",
    x = "Full-coordinate LL-PME", y = "PRIME"
  ) +
  paper_theme() +
  theme(
    strip.background = element_rect(
      fill = "grey94", color = "grey35", linewidth = 0.35
    ),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = grid::unit(0.7, "lines")
  )

prediction_long <- rbind(
  data.frame(
    ratio[c("dataset", "repeat_id", "fold")],
    feature_set = "All variables",
    rmse = ratio$full_rmse
  ),
  data.frame(
    ratio[c("dataset", "repeat_id", "fold")],
    feature_set = "PRIME",
    rmse = ratio$selected_rmse
  )
)
prediction_long$feature_set <- factor(
  prediction_long$feature_set,
  levels = c("All variables", "PRIME")
)
prediction_long$pair_id <- interaction(
  prediction_long$dataset,
  prediction_long$repeat_id,
  prediction_long$fold,
  drop = TRUE
)
prediction_median <- aggregate(
  rmse ~ dataset + feature_set,
  prediction_long,
  median
)
p_prediction <- ggplot(
  prediction_long,
  aes(feature_set, rmse, group = pair_id)
) +
  geom_line(
    color = "grey70", linewidth = 0.35, alpha = 0.65
  ) +
  geom_point(
    aes(color = dataset), size = 1.35, alpha = 0.62
  ) +
  geom_line(
    data = prediction_median,
    aes(feature_set, rmse, group = dataset),
    inherit.aes = FALSE,
    linewidth = 0.75, color = "black"
  ) +
  geom_point(
    data = prediction_median,
    aes(feature_set, rmse),
    inherit.aes = FALSE,
    shape = 23, size = 2.6, fill = "white", color = "black", stroke = 0.55
  ) +
  facet_wrap(
    ~dataset, nrow = 1, scales = "free_y",
    labeller = as_labeller(dataset_labels)
  ) +
  scale_color_manual(values = dataset_palette) +
  labs(
    title = "B  Held-out prediction",
    x = NULL,
    y = "Test RMSE"
  ) +
  paper_theme() +
  theme(
    strip.background = element_rect(
      fill = "grey94", color = "grey35", linewidth = 0.35
    ),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = grid::unit(0.7, "lines")
  )

grDevices::cairo_pdf(
  file.path(figure_directory, "observed_application_main.pdf"),
  width = 8.0, height = 5.8
)
grid::grid.newpage()
layout <- grid::grid.layout(
  nrow = 2L, ncol = 1L, heights = grid::unit(c(1.35, 1), "null")
)
grid::pushViewport(grid::viewport(layout = layout))
print(p_allocation, vp = grid::viewport(layout.pos.row = 1L))
print(p_prediction, vp = grid::viewport(layout.pos.row = 2L))
grid::popViewport()
grDevices::dev.off()

all_runs <- rbind(
  observed_runs[c("dataset", "method", "selected_size", "value_calls", "elapsed_sec")],
  transform(
    concrete_runs[c("method", "selected_size", "value_calls", "elapsed_sec")],
    dataset = "concrete"
  )[c("dataset", "method", "selected_size", "value_calls", "elapsed_sec")]
)
dataset_dimensions <- data.frame(
  dataset = dataset_order,
  n = c(1030L, 1503L, 392L),
  p = c(8L, 5L, 8L),
  stringsAsFactors = FALSE
)
summary_rows <- lapply(dataset_order, function(dataset) {
  run_cell <- all_runs[all_runs$dataset == dataset, ]
  screen <- run_cell[run_cell$method == "PRIME", ]
  full <- run_cell[run_cell$method == "Full LL-PME", ]
  tv_cell <- tv[as.character(tv$dataset) == dataset, ]
  ratio_cell <- ratio[as.character(ratio$dataset) == dataset, ]
  dimensions <- dataset_dimensions[dataset_dimensions$dataset == dataset, ]
  data.frame(
    dataset = dataset,
    n = dimensions$n,
    p = dimensions$p,
    selected_dimension = stats::median(screen$selected_size),
    tv = stats::median(tv_cell$total_variation),
    screened_calls = stats::median(screen$value_calls),
    full_calls = stats::median(full$value_calls),
    rmse_ratio = stats::median(ratio_cell$rmse_ratio),
    stringsAsFactors = FALSE
  )
})
application_summary <- do.call(rbind, summary_rows)
write.csv(
  application_summary,
  "results/summary/observed_application_summary.csv",
  row.names = FALSE
)

english_rows <- vapply(seq_len(nrow(application_summary)), function(index) {
  row <- application_summary[index, ]
  sprintf(
    "%s & %d & %d & %.0f & %.0f & %.0f & %.3f \\\\",
    dataset_labels[row$dataset], row$n, row$p, row$selected_dimension,
    row$screened_calls, row$full_calls, row$rmse_ratio
  )
}, character(1))
writeLines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Observed-data comparison.  Selected dimension and subset-value counts are bootstrap medians.  The test-RMSE ratio is PRIME divided by the full-variable RMSE, summarized by the median over 25 held-out folds.}",
    "\\label{tab:observed-application}",
    "\\small",
    "\\begin{tabular}{lrrrrrr}",
    "\\toprule",
    "Data & $n$ & $p$ & Selected $d$ & \\multicolumn{2}{c}{Subset-value evaluations} & Test-RMSE ratio \\\\",
    "\\cmidrule(lr){5-6}",
    " & & & & PRIME & Full-coordinate & \\\\",
    "\\midrule",
    english_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "observed_application_main.tex")
)

effect_summary <- function(data) {
  pieces <- lapply(
    split(data, interaction(
      data$dataset, data$method, data$variable_name,
      drop = TRUE, lex.order = TRUE
    )),
    function(cell) data.frame(
      dataset = cell$dataset[1L],
      method = cell$method[1L],
      variable_name = cell$variable_name[1L],
      selection_frequency = mean(cell$selected),
      median = stats::median(cell$estimate),
      q05 = unname(stats::quantile(cell$estimate, 0.05)),
      q95 = unname(stats::quantile(cell$estimate, 0.95)),
      stringsAsFactors = FALSE
    )
  )
  do.call(rbind, pieces)
}
observed_effects <- effect_summary(observed_phi[
  observed_phi$method %in% c("PRIME", "Full LL-PME"),
])
write.csv(
  observed_effects,
  "results/summary/observed_application_effects.csv",
  row.names = FALSE
)

variable_labels <- c(
  Frequency = "Frequency", Angle = "Angle of attack",
  Chord_Length = "Chord length", Velocity = "Velocity",
  Displacement = "Displacement thickness", cylinders = "Cylinders",
  displacement = "Displacement", horsepower = "Horsepower",
  weight = "Weight", acceleration = "Acceleration",
  model_year = "Model year", origin_2 = "Origin: Europe",
  origin_3 = "Origin: Japan"
)
format_effect_table <- function(
    dataset, caption, label,
    header = "Variable & Selection frequency & PRIME & Full-coordinate LL-PME \\\\"
) {
  data <- observed_effects[observed_effects$dataset == dataset, ]
  screened <- data[data$method == "PRIME", ]
  full <- data[data$method == "Full LL-PME", ]
  screened <- screened[order(screened$variable_name), ]
  full <- full[match(screened$variable_name, full$variable_name), ]
  rows <- vapply(seq_len(nrow(screened)), function(index) {
    name <- variable_labels[screened$variable_name[index]]
    if (is.na(name)) name <- screened$variable_name[index]
    sprintf(
      "%s & %.2f & %.3f [%.3f, %.3f] & %.3f [%.3f, %.3f] \\\\",
      name, screened$selection_frequency[index],
      screened$median[index], screened$q05[index], screened$q95[index],
      full$median[index], full$q05[index], full$q95[index]
    )
  }, character(1))
  c(
    "\\begin{table}[H]", "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), "\\small",
    "\\begin{tabular}{lrrr}", "\\toprule",
    header,
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"
  )
}
observed_tables <- c(
  format_effect_table(
    "airfoil", "Airfoil PME over 200 bootstrap samples.  Allocation entries are medians [5th, 95th percentiles].",
    "tab:airfoil-allocation"
  ),
  format_effect_table(
    "auto_mpg", "Auto MPG PME over 200 bootstrap samples.  Allocation entries are medians [5th, 95th percentiles].",
    "tab:auto-mpg-allocation"
  )
)
writeLines(
  observed_tables,
  file.path(table_directory, "observed_allocation_supplement.tex")
)
