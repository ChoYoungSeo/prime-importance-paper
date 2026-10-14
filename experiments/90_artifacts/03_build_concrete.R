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
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.15, "lines"),
      plot.margin = margin(5, 7, 5, 5)
    )
}

bootstrap <- read.csv(
  "results/raw/concrete_v052_main/bootstrap_phi.csv",
  stringsAsFactors = FALSE
)
runs <- read.csv(
  "results/raw/concrete_v052_main/bootstrap_runs.csv",
  stringsAsFactors = FALSE
)

summarize_effects <- function(data) {
  pieces <- lapply(
    split(data, interaction(
      data$representation, data$method, data$variable_name,
      drop = TRUE, lex.order = TRUE
    )),
    function(cell) {
      data.frame(
        representation = cell$representation[1L],
        method = cell$method[1L],
        variable_name = cell$variable_name[1L],
        median = stats::median(cell$estimate),
        q05 = unname(stats::quantile(cell$estimate, 0.05)),
        q95 = unname(stats::quantile(cell$estimate, 0.95)),
        selection_frequency = mean(cell$selected),
        stringsAsFactors = FALSE
      )
    }
  )
  do.call(rbind, pieces)
}
effect_summary <- summarize_effects(bootstrap)
write.csv(
  effect_summary,
  "results/summary/concrete_pme_bootstrap_summary.csv",
  row.names = FALSE
)

engineering <- effect_summary[
  effect_summary$representation == "engineering", , drop = FALSE
]
method_order <- c("PRIME", "Full LL-PME", "PME-kNN")
engineering$method <- factor(
  engineering$method,
  levels = method_order,
  labels = c("PRIME", "Full-coordinate LL-PME", "PME-kNN")
)
variable_order <- c(
  "fineagg", "coarseagg", "slag", "ash", "superplastic", "cement",
  "water_binder_ratio", "age"
)
variable_label <- c(
  fineagg = "Fine aggregate",
  coarseagg = "Coarse aggregate",
  slag = "Blast-furnace slag",
  ash = "Fly ash",
  superplastic = "Superplasticizer",
  cement = "Cement",
  water_binder_ratio = "Water-to-binder ratio",
  age = "Age"
)
engineering$variable_index <- match(engineering$variable_name, variable_order)
method_offset <- c(
  "PRIME" = -0.20,
  "Full-coordinate LL-PME" = 0,
  "PME-kNN" = 0.20
)
engineering$plot_y <- engineering$variable_index +
  unname(method_offset[as.character(engineering$method)])
method_palette <- c(
  "PRIME" = "#0072B2",
  "Full-coordinate LL-PME" = "#7A5195",
  "PME-kNN" = "#D55E00"
)
p_concrete <- ggplot(engineering, aes(color = method)) +
  geom_segment(
    aes(x = q05, xend = q95, y = plot_y, yend = plot_y),
    linewidth = 0.65
  ) +
  geom_point(aes(x = median, y = plot_y), size = 2.1) +
  scale_color_manual(values = method_palette) +
  scale_x_continuous(
    limits = c(0, 0.82), breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = seq_along(variable_order),
    labels = unname(variable_label[variable_order]),
    limits = c(0.55, 8.45)
  ) +
  labs(x = "PME allocation", y = NULL) +
  paper_theme()
ggsave(
  file.path(figure_directory, "concrete_pme_main.pdf"),
  p_concrete, width = 7.3, height = 4.5, units = "in",
  device = grDevices::cairo_pdf
)

raw <- effect_summary[
  effect_summary$representation == "raw" &
    effect_summary$method == "PRIME",
  ,
  drop = FALSE
]
raw <- raw[match(c(
  "cement", "slag", "ash", "water", "superplastic", "coarseagg",
  "fineagg", "age"
), raw$variable_name), ]
raw_labels <- c(
  cement = "Cement", slag = "Blast-furnace slag", ash = "Fly ash",
  water = "Water", superplastic = "Superplasticizer",
  coarseagg = "Coarse aggregate", fineagg = "Fine aggregate", age = "Age"
)
raw_rows <- vapply(seq_len(nrow(raw)), function(index) {
  sprintf(
    "%s & %.3f & %.3f [%.3f, %.3f] \\\\",
    raw_labels[raw$variable_name[index]], raw$selection_frequency[index],
    raw$median[index], raw$q05[index], raw$q95[index]
  )
}, character(1))
writeLines(
  c(
    "\\begin{table}[H]",
    "\\centering",
    "\\caption{Concrete analysis in the original coordinate system over 200 bootstrap samples.}",
    "\\label{tab:concrete-raw}",
    "\\small",
    "\\begin{tabular}{lrr}",
    "\\toprule",
    "Variable & Selection frequency & PME median [90\\% interval] \\\\",
    "\\midrule",
    raw_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ),
  file.path(table_directory, "concrete_raw_supplement.tex")
)

engineering_runs <- runs[runs$representation == "engineering", ]
call_summary <- aggregate(value_calls ~ method, engineering_runs, stats::median)
time_summary <- aggregate(elapsed_sec ~ method, engineering_runs, stats::median)
write.csv(
  merge(call_summary, time_summary, by = "method"),
  "results/summary/concrete_computation_summary.csv",
  row.names = FALSE
)
