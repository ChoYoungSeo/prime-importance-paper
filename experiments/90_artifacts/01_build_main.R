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
      strip.background = element_rect(fill = "grey94", color = "grey35",
                                      linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.25, "lines"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_figure <- function(plot, filename, width, height) {
  ggsave(
    file.path(figure_directory, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}

summary_table <- function(data, groups, metrics) {
  split_key <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(data, split_key), function(cell) {
    identifiers <- cell[1L, groups, drop = FALSE]
    values <- unlist(lapply(metrics, function(metric) {
      x <- cell[[metric]]
      c(
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q25 = unname(quantile(x, 0.25, na.rm = TRUE)),
        q75 = unname(quantile(x, 0.75, na.rm = TRUE))
      )
    }))
    names(values) <- unlist(lapply(metrics, function(metric) {
      paste0(metric, "_", c("mean", "median", "q25", "q75"))
    }))
    cbind(
      identifiers,
      as.data.frame(as.list(values)),
      replications = nrow(cell)
    )
  })
  do.call(rbind, pieces)
}

conditional_variance <- read.csv(
  "results/provenance/inherited/conditional_variance_summary.csv",
  stringsAsFactors = FALSE
)
conditional_wide <- reshape(
  conditional_variance[c(
    "benchmark", "n", "rho", "method", "rmse_normalized_tau"
  )],
  idvar = c("benchmark", "n", "rho"),
  timevar = "method",
  direction = "wide"
)
conditional_wide$rmse_ratio <-
  conditional_wide[["rmse_normalized_tau.Local linear"]] /
  conditional_wide[["rmse_normalized_tau.Local constant"]]
conditional_wide$n_label <- factor(
  conditional_wide$n,
  levels = c(100, 300),
  labels = c("n = 100", "n = 300")
)
conditional_wide$rho_label <- factor(
  conditional_wide$rho,
  levels = c(0, 0.5),
  labels = c("Independent", "AR(1), rho = 0.5")
)
p_conditional <- ggplot(
  conditional_wide,
  aes(
    x = `rmse_normalized_tau.Local constant`,
    y = `rmse_normalized_tau.Local linear`
  )
) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.45) +
  geom_point(shape = 16, size = 2.4, color = "black") +
  facet_wrap(~n_label, nrow = 1) +
  scale_x_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_equal() +
  labs(
    x = "Local-constant normalized RMSE",
    y = "Local-linear normalized RMSE"
  ) +
  paper_theme() +
  theme(legend.position = "none")
save_figure(
  p_conditional,
  "conditional_variance_rmse.pdf",
  6.8,
  3.4
)

selection <- read.csv(
  "results/raw/selection_v052_runs.csv",
  stringsAsFactors = FALSE
)
selection$method <- factor(
  selection$method,
  levels = c("PRIME screen", "FIRST", "FORD", "FOCI", "LASSO"),
  labels = c("PRIME screen", "FIRST", "FORD", "FOCI", "LASSO")
)
selection$benchmark <- factor(
  selection$benchmark,
  levels = c("branin", "ackley3", "friedman"),
  labels = c("Branin", "Ackley-3", "Friedman")
)
selection$design <- factor(
  selection$design,
  levels = c("independent", "spaced_ar1", "balanced_proxy"),
  labels = c("Independent", "Spaced AR(1)", "Correlated proxies")
)
selection$n_label <- factor(
  selection$n, levels = c(100, 200, 500, 1000, 2000)
)
selection_cells <- aggregate(
  cbind(inclusion, exact, selected_size) ~
    benchmark + design + n_label + method,
  selection,
  mean
)
selection_palette <- c(
  "PRIME screen" = "#0072B2",
  "FIRST" = "#D55E00",
  "FORD" = "#009E73",
  "FOCI" = "#CC79A7",
  "LASSO" = "#4D4D4D"
)
selection_plot <- function(metric, y_label) {
  ggplot(
    selection_cells,
    aes(
      n_label, y = .data[[metric]], color = method, shape = method,
      linetype = method, group = method
    )
  ) +
    geom_line(linewidth = 0.52) +
    geom_point(size = 1.75, stroke = 0.45) +
    facet_grid(design ~ benchmark) +
    scale_color_manual(values = selection_palette) +
    scale_shape_manual(values = c(16, 1, 17, 18, 15)) +
    scale_linetype_manual(values = c(
      "solid", "dashed", "dotdash", "dotted", "longdash"
    )) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(x = "Sample size", y = y_label) +
    paper_theme() +
    theme(legend.position = "top")
}
p_selection <- selection_plot("inclusion", "Active-set coverage")
save_figure(p_selection, "selection_coverage_main.pdf", 8.2, 6.7)
p_selection_exact <- selection_plot("exact", "Exact active-set recovery")
save_figure(
  p_selection_exact, "selection_exact_supplement.pdf", 8.2, 6.7
)
write.csv(
  summary_table(
    selection,
    c("benchmark", "design", "n_label", "method"),
    c("tpr", "inclusion", "exact", "false_positives", "selected_size",
      "elapsed_sec")
  ),
  "results/summary/selection_summary.csv",
  row.names = FALSE
)

selector <- read.csv(
  "results/raw/selection_to_pme_v052_main/selection_to_pme.csv",
  stringsAsFactors = FALSE
)
selector$method <- factor(
  selector$method,
  levels = c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD"),
  labels = c("PRIME screen", "FIRST", "LASSO", "FOCI", "FORD")
)
selector$benchmark <- factor(
  selector$benchmark,
  levels = c("branin", "ackley3", "friedman"),
  labels = c("Branin", "Ackley-3", "Friedman")
)
selector$design <- factor(
  selector$design,
  levels = c("independent", "balanced_proxy"),
  labels = c("Independent", "Correlated proxies")
)
selector_plot <- droplevels(selector[selector$method != "LASSO", ])
selector_palette <- c(
  "PRIME screen" = "#0072B2",
  "FIRST" = "#E69F00",
  "FOCI" = "#CC79A7",
  "FORD" = "#D55E00"
)
p_selector <- ggplot(
  selector_plot[selector_plot$feasible, ],
  aes(method, half_l1, fill = method)
) +
  geom_boxplot(width = 0.62, outlier.shape = 21, outlier.size = 0.9,
               linewidth = 0.35) +
  facet_grid(design ~ benchmark) +
  scale_fill_manual(values = selector_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = expression(frac(1, 2) * L[1]~PME~error)) +
  paper_theme() +
  theme(axis.text.x = element_text(angle = 24, hjust = 1))
save_figure(p_selector, "selection_to_pme_main.pdf", 8.0, 5.0)
write.csv(
  summary_table(selector, c("benchmark", "design", "method"),
                c("half_l1", "inactive_mass", "selected_size", "feasible")),
  "results/summary/selection_to_pme_summary.csv",
  row.names = FALSE
)

low <- read.csv("results/raw/pme_v052_main/pme_runs.csv", stringsAsFactors = FALSE)
low <- low[low$method %in% c(
  "PRIME", "Full local-linear PME", "PME-kNN", "NCMD-PME"
), ]
low$method <- factor(
  low$method,
  levels = c(
    "PRIME", "Full local-linear PME", "PME-kNN", "NCMD-PME"
  ),
  labels = c("PRIME", "Full-coordinate LL-PME", "PME-kNN", "NCMD-PME")
)
low$benchmark <- factor(low$benchmark, levels = c("branin", "bratley"),
                        labels = c("Branin", "Bratley"))
low$rho_label <- factor(
  low$rho,
  levels = c(0, 0.5, 0.9),
  labels = c(expression(rho == 0), expression(rho == 0.5), expression(rho == 0.9))
)
low$n_label <- factor(low$n, levels = c(100, 250, 500))
low_palette <- c(
  "PRIME" = "#0072B2",
  "Full-coordinate LL-PME" = "#7A5195",
  "PME-kNN" = "#D55E00",
  "NCMD-PME" = "#009E73"
)
p_low <- ggplot(low, aes(n_label, half_l1, fill = method)) +
  geom_boxplot(position = position_dodge2(width = 0.78, preserve = "single"),
               width = 0.68, outlier.shape = NA, linewidth = 0.32) +
  facet_grid(rho_label ~ benchmark, labeller = label_parsed) +
  scale_fill_manual(values = low_palette) +
  scale_y_sqrt(breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1),
               expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "A  Allocation accuracy",
    x = "Sample size",
    y = expression(frac(1, 2) * L[1]~PME~error)
  ) +
  paper_theme() +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))

low_runtime <- summary_table(
  low,
  c("benchmark", "n_label", "method"),
  "elapsed_sec"
)
p_low_runtime <- ggplot(
  low_runtime,
  aes(n_label, elapsed_sec_median, color = method,
      shape = method, group = method)
) +
  geom_linerange(
    aes(ymin = elapsed_sec_q25, ymax = elapsed_sec_q75),
    position = position_dodge(width = 0.08), linewidth = 0.42
  ) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 1.8, stroke = 0.45) +
  facet_wrap(~benchmark, nrow = 1) +
  scale_color_manual(values = low_palette) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10)) +
  labs(
    title = "B  Elapsed time",
    x = "Sample size",
    y = "Seconds (log scale)"
  ) +
  paper_theme() +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    legend.position = "none"
  )

grDevices::cairo_pdf(
  file.path(figure_directory, "pme_estimator_calibration.pdf"),
  width = 8.0,
  height = 8.5
)
grid::grid.newpage()
layout <- grid::grid.layout(
  nrow = 2L,
  ncol = 1L,
  heights = grid::unit(c(2.25, 1), "null")
)
grid::pushViewport(grid::viewport(layout = layout))
print(p_low, vp = grid::viewport(layout.pos.row = 1L))
print(
  p_low_runtime,
  vp = grid::viewport(layout.pos.row = 2L)
)
grid::popViewport()
grDevices::dev.off()
write.csv(
  summary_table(low, c("benchmark", "n", "rho", "method"),
                c("half_l1", "inactive_mass", "proxy_mass", "elapsed_sec")),
  "results/summary/pme_estimator_calibration_summary.csv",
  row.names = FALSE
)

corr <- read.csv(
  "results/raw/correlated_pipelines_v052_main/correlated_pipeline_runs.csv",
  stringsAsFactors = FALSE
)
method_levels <- c(
  "PRIME",
  "FIRST screen + PME-kNN",
  "NNVS screen + PME-kNN",
  "LASSO screen + PME-kNN",
  "FOCI screen + PME-kNN",
  "FORD screen + PME-kNN"
)
method_labels <- c(
  "PRIME", "FIRST + PME-kNN", "NNVS + PME-kNN",
  "LASSO + PME-kNN", "FOCI + PME-kNN", "FORD + PME-kNN"
)
corr$method <- factor(corr$method, levels = method_levels, labels = method_labels)
corr$benchmark <- factor(corr$benchmark, levels = c("branin", "bratley"),
                         labels = c("Branin", "Bratley"))
corr$n_label <- factor(corr$n, levels = c(100, 200, 500))
pipeline_palette <- c(
  "PRIME" = "#0072B2",
  "FIRST + PME-kNN" = "#D55E00",
  "NNVS + PME-kNN" = "#009E73",
  "LASSO + PME-kNN" = "#CC79A7",
  "FOCI + PME-kNN" = "#E69F00",
  "FORD + PME-kNN" = "#4D4D4D"
)
write.csv(
  summary_table(corr, c("benchmark", "n", "method"),
                c("half_l1", "proxy_mass", "selected_size", "elapsed_sec")),
  "results/summary/pme_correlated_pipeline_summary.csv",
  row.names = FALSE
)

friedman_phi <- read.csv(
  "results/raw/friedman_detail_v052_main/friedman_detail_phi.csv",
  stringsAsFactors = FALSE
)
friedman_method_levels <- c("PRIME screen", "FIRST", "FOCI", "FORD")
friedman_method_labels <- c("PRIME screen", "FIRST", "FOCI", "FORD")
friedman_phi <- friedman_phi[is.finite(friedman_phi$estimate), , drop = FALSE]
friedman_phi$method <- factor(
  friedman_phi$method,
  levels = friedman_method_levels,
  labels = friedman_method_labels
)
friedman_active <- friedman_phi[friedman_phi$role == "active", ]
friedman_active$coordinate <- paste0("X", friedman_active$variable)
aggregate_friedman_role <- function(role, label) {
  rows <- friedman_phi[friedman_phi$role == role, ]
  out <- aggregate(cbind(estimate, truth) ~ detail_id + method, rows, sum)
  out$coordinate <- label
  out
}
friedman_allocation <- rbind(
  friedman_active[c(
    "detail_id", "method", "estimate", "truth", "coordinate"
  )],
  aggregate_friedman_role("correlated proxy", "Proxies"),
  aggregate_friedman_role("independent noise", "Other")
)
friedman_allocation$coordinate <- factor(
  friedman_allocation$coordinate,
  levels = c(paste0("X", 1:5), "Proxies", "Other")
)
friedman_truth <- unique(friedman_allocation[c("coordinate", "truth")])
p_friedman_allocation <- ggplot(
  friedman_allocation,
  aes(coordinate, estimate)
) +
  geom_boxplot(
    width = 0.58, fill = "grey80", outlier.shape = NA, linewidth = 0.34
  ) +
  geom_point(
    data = friedman_truth,
    aes(coordinate, truth),
    inherit.aes = FALSE,
    shape = 23, size = 2.6, fill = "#E69F00", color = "black",
    stroke = 0.55
  ) +
  facet_wrap(~method, nrow = 1) +
  scale_y_sqrt(
    limits = c(0, 1),
    breaks = c(0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = NULL, y = "PME allocation") +
  paper_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 24, hjust = 1),
    panel.spacing = grid::unit(1.0, "lines")
  )
save_figure(
  p_friedman_allocation,
  "pme_friedman_allocation_main.pdf",
  9.0,
  3.7
)

scaling_file <- "results/raw/scaling_v052_main/scaling_runs.csv"
if (file.exists(scaling_file)) {
  scaling <- read.csv(scaling_file, stringsAsFactors = FALSE)
  scaling$method <- factor(scaling$method, levels = method_levels,
                           labels = method_labels)
  scaling$benchmark <- factor(
    scaling$benchmark,
    levels = c("branin", "friedman", "otl", "piston", "borehole"),
    labels = c("Branin", "Friedman", "OTL", "Piston", "Borehole")
  )
  cells <- summary_table(
    scaling,
    c("benchmark", "p", "method"),
    c("half_l1", "elapsed_sec", "selected_size")
  )
  cells$p <- as.numeric(as.character(cells$p))
  plot_cells <- cells[
    as.character(cells$benchmark) %in% c("OTL", "Piston", "Borehole"),
    , drop = FALSE
  ]
  p_scaling_error <- ggplot(
    plot_cells,
    aes(p, half_l1_median, color = method, shape = method, group = method)
  ) +
    geom_linerange(aes(ymin = half_l1_q25, ymax = half_l1_q75),
                   position = position_dodge(width = 0.06), linewidth = 0.45) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 2.0, stroke = 0.45) +
    facet_wrap(~benchmark, nrow = 1, scales = "free_y") +
    scale_x_log10(breaks = c(20, 50, 100, 250, 500)) +
    scale_color_manual(values = pipeline_palette) +
    scale_shape_manual(values = c(16, 17, 15, 18, 3, 4)) +
    labs(
      x = "Ambient dimension",
      y = expression(frac(1, 2) * L[1]~PME~error)
    ) +
    paper_theme() +
    theme(
      legend.key.width = grid::unit(0.85, "lines")
    ) +
    guides(
      color = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE)
    )
  save_figure(p_scaling_error, "pme_scaling_main.pdf", 8.2, 3.2)
  write.csv(cells, "results/summary/pme_scaling_summary.csv", row.names = FALSE)
}

observed_directory <- "results/raw/observed_v052_main"
if (file.exists(file.path(observed_directory, "observed_runs.csv"))) {
  observed_runs <- read.csv(
    file.path(observed_directory, "observed_runs.csv"),
    stringsAsFactors = FALSE
  )
  write.csv(
    summary_table(
      observed_runs, "dataset",
      c("original_selected", "controls_selected", "control_mass", "elapsed_sec")
    ),
    "results/summary/observed_pme_summary.csv",
    row.names = FALSE
  )
}

twinning_file <- "results/raw/twinning_v052_main/twinning_runs.csv"
if (file.exists(twinning_file)) {
  twinning <- read.csv(twinning_file, stringsAsFactors = FALSE)
  twinning <- twinning[twinning$twin_r != "all", , drop = FALSE]
  twinning$twin_r <- factor(twinning$twin_r, levels = c("2", "3", "5"))
  write.csv(
    summary_table(
      twinning, "twin_r",
      c("active_recall", "selected_jaccard", "candidate_jaccard", "speedup")
    ),
    "results/summary/twinning_v052_summary.csv",
    row.names = FALSE
  )
}
