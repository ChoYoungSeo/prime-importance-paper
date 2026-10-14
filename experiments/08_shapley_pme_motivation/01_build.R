library(ggplot2)

source("software/prime/R/constants.R")
source("software/prime/R/pme.R")

exact_shapley <- function(values) {
  values <- as.numeric(values)
  dimension <- as.integer(round(log2(length(values))))
  effects <- numeric(dimension)
  for (coordinate in seq_len(dimension)) {
    bit <- bitwShiftL(1L, coordinate - 1L)
    for (mask in 0:(2^dimension - 1L)) {
      if (bitwAnd(mask, bit) != 0L) next
      size <- sum(as.logical(intToBits(mask)[seq_len(dimension)]))
      weight <- factorial(size) * factorial(dimension - size - 1L) /
        factorial(dimension)
      effects[coordinate] <- effects[coordinate] + weight *
        (values[bitwOr(mask, bit) + 1L] - values[mask + 1L])
    }
  }
  effects
}

closed_effects <- function(beta, rho) {
  total_variance <- 2 + beta^2 + 2 * rho * beta
  c(
    0,
    (1 + beta * rho)^2,
    1,
    (1 + beta * rho)^2 + 1,
    (rho + beta)^2,
    1 + beta^2 + 2 * rho * beta,
    1 + (rho + beta)^2,
    total_variance
  ) / total_variance
}

total_effects <- function(closed) {
  dimension <- as.integer(round(log2(length(closed))))
  full_mask <- 2^dimension - 1L
  vapply(0:full_mask, function(mask) {
    1 - closed[bitwXor(full_mask, mask) + 1L]
  }, numeric(1))
}

analytic_allocations <- function(beta, rho) {
  total_variance <- 2 + beta^2 + 2 * rho * beta
  shared_variance <- 1 + beta^2 + 2 * rho * beta
  shapley <- c(
    1 + rho * beta - 0.5 * rho^2 * (1 - beta^2),
    1,
    beta^2 + rho * beta + 0.5 * rho^2 * (1 - beta^2)
  ) / total_variance
  pme <- c(
    shared_variance / (1 + beta^2),
    1,
    beta^2 * shared_variance / (1 + beta^2)
  ) / total_variance
  list(shapley = shapley, pme = pme)
}

rho_values <- c(0, 0.3, 0.6, 0.9)
beta_grid <- sort(unique(c(seq(0, 1, length.out = 501L), 0.001, 0.01, 0.1)))

result <- do.call(rbind, lapply(rho_values, function(rho) {
  do.call(rbind, lapply(beta_grid, function(beta) {
    allocation <- analytic_allocations(beta, rho)
    data.frame(
      beta = beta,
      rho = rho,
      method = c("Shapley effects", "PME"),
      effect = c(allocation$shapley[3], allocation$pme[3])
    )
  }))
}))

write.csv(
  result,
  "results/summary/shapley_pme_proxy_population.csv",
  row.names = FALSE
)

validation_grid <- expand.grid(
  beta = seq(0, 1, by = 0.05),
  rho = rho_values
)
validation_error <- apply(validation_grid, 1L, function(row) {
  beta <- row[["beta"]]
  rho <- row[["rho"]]
  closed <- closed_effects(beta, rho)
  formula <- analytic_allocations(beta, rho)
  max(abs(c(
    exact_shapley(closed) - formula$shapley,
    as.numeric(prime_allocate(total_effects(closed))) - formula$pme
  )))
})

for (rho in rho_values) {
  endpoint_zero <- analytic_allocations(0, rho)
  endpoint_one <- analytic_allocations(1, rho)
  pme_x3 <- result[
    result$method == "PME" & result$rho == rho, "effect"
  ]
  stopifnot(
    abs(endpoint_zero$pme[3]) < 1e-12,
    abs(endpoint_zero$shapley[3] - rho^2 / 4) < 1e-12,
    max(abs(endpoint_one$shapley - endpoint_one$pme)) < 1e-12,
    min(diff(pme_x3)) >= -1e-12
  )
}
stopifnot(max(validation_error) < 1e-10)

weak_effect <- analytic_allocations(beta = 1e-3, rho = 0.9)
stopifnot(
  abs(weak_effect$shapley[3] - 0.202767705181) < 1e-12,
  abs(weak_effect$pme[3] - 5.00449344465e-7) < 1e-18
)

rho_labels <- sprintf("%.1f", rho_values)
result$rho_label <- factor(
  sprintf("%.1f", result$rho),
  levels = rho_labels
)

rho_palette <- setNames(
  c("#333333", "#0072B2", "#D55E00", "#009E73"),
  rho_labels
)
rho_linetype <- setNames(c("solid", "dashed", "dotdash", "longdash"), rho_labels)

make_panel <- function(method, title, show_legend) {
  panel_data <- result[result$method == method, , drop = FALSE]
  ggplot(
    panel_data,
    aes(
      x = beta,
      y = effect,
      color = rho_label,
      linetype = rho_label,
      group = rho_label
    )
  ) +
    geom_line(linewidth = 0.82) +
    scale_color_manual(values = rho_palette, name = expression(rho)) +
    scale_linetype_manual(values = rho_linetype, name = expression(rho)) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      breaks = seq(0, 0.45, 0.1),
      limits = c(0, 0.45),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = expression(beta),
      y = expression("Allocation to " * X[3])
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "left"
    )
}

shapley_panel <- make_panel("Shapley effects", "A  Shapley effects", TRUE)
pme_panel <- make_panel("PME", "B  Proportional marginal effects", FALSE) +
  labs(y = NULL)

grDevices::cairo_pdf(
  "manuscript/figures/shapley_pme_proxy_population.pdf",
  width = 6.8,
  height = 3.0
)
grid::grid.newpage()
layout <- grid::grid.layout(nrow = 1L, ncol = 2L)
grid::pushViewport(grid::viewport(layout = layout))
print(shapley_panel, vp = grid::viewport(layout.pos.col = 1L))
print(pme_panel, vp = grid::viewport(layout.pos.col = 2L))
grid::popViewport()
grDevices::dev.off()

cat("Built and validated the Gaussian Shapley-PME transition figure.\n")
