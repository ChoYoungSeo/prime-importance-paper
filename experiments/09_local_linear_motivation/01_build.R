library(ggplot2)

mean_function <- function(x) x + x^2

h0 <- 0.65
x_grid <- seq(-0.85, 0.85, length.out = 401L)
geometry <- data.frame(
  x = x_grid,
  response = mean_function(x_grid),
  local_constant = h0^2 / 3,
  local_linear = x_grid + h0^2 / 3
)

residual_x <- seq(-0.6, 0.6, by = 0.2)
residuals <- data.frame(
  x = residual_x,
  response = mean_function(residual_x),
  local_constant = h0^2 / 3,
  local_linear = residual_x + h0^2 / 3
)

base_geometry <- function() {
  list(
    annotate(
      "rect",
      xmin = -h0,
      xmax = h0,
      ymin = -Inf,
      ymax = Inf,
      fill = "#E8E8E8",
      alpha = 0.55
    ),
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey45"),
    scale_x_continuous(
      breaks = c(-0.5, 0, 0.5),
      expand = expansion(mult = c(0, 0.01))
    ),
    scale_y_continuous(
      breaks = c(-0.25, 0.25, 0.75, 1.25),
      expand = expansion(mult = c(0, 0.01))
    ),
    coord_cartesian(xlim = c(-0.85, 0.85), ylim = c(-0.4, 1.55)),
    theme_classic(base_size = 9.5),
    theme(
      plot.title = element_text(face = "bold", size = 9.5),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 8.5),
      legend.key.width = grid::unit(0.9, "lines"),
      plot.margin = margin(4, 6, 2, 4)
    )
  )
}

panel_a <- ggplot(geometry, aes(x = x)) +
  base_geometry() +
  geom_segment(
    data = residuals,
    aes(x = x, xend = x, y = response, yend = local_constant),
    inherit.aes = FALSE,
    color = "#D55E00",
    alpha = 0.48,
    linewidth = 0.55
  ) +
  geom_line(aes(y = response, color = "mean"), linewidth = 0.85) +
  geom_line(
    aes(y = local_constant, color = "constant fit"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(
    data = residuals,
    aes(x = x, y = response),
    inherit.aes = FALSE,
    size = 1.1,
    color = "#222222"
  ) +
  scale_color_manual(
    values = c("mean" = "#222222", "constant fit" = "#D55E00"),
    breaks = c("mean", "constant fit"),
    labels = c(expression(m(u)), "constant fit"),
    name = NULL
  ) +
  labs(
    title = "A  Local-constant residuals",
    x = "Centered displacement",
    y = "Conditional mean"
  )

panel_b <- ggplot(geometry, aes(x = x)) +
  base_geometry() +
  geom_segment(
    data = residuals,
    aes(x = x, xend = x, y = response, yend = local_linear),
    inherit.aes = FALSE,
    color = "#0072B2",
    alpha = 0.48,
    linewidth = 0.55
  ) +
  geom_line(aes(y = response, color = "mean"), linewidth = 0.85) +
  geom_line(
    aes(y = local_linear, color = "linear fit"),
    linewidth = 0.8,
    linetype = "dotdash"
  ) +
  geom_point(
    data = residuals,
    aes(x = x, y = response),
    inherit.aes = FALSE,
    size = 1.1,
    color = "#222222"
  ) +
  scale_color_manual(
    values = c("mean" = "#222222", "linear fit" = "#0072B2"),
    breaks = c("mean", "linear fit"),
    labels = c(expression(m(u)), "linear fit"),
    name = NULL
  ) +
  labs(
    title = "B  Local-linear residuals",
    x = "Centered displacement",
    y = NULL
  )

h_grid <- seq(0.02, 1, length.out = 400L)
drift <- rbind(
  data.frame(
    h = h_grid,
    method = "Local constant",
    contribution = h_grid^2 / 3 + 4 * h_grid^4 / 45
  ),
  data.frame(
    h = h_grid,
    method = "Local linear",
    contribution = 4 * h_grid^4 / 45
  )
)

panel_c <- ggplot(
  drift,
  aes(x = h, y = contribution, color = method, linetype = method)
) +
  geom_line(linewidth = 0.85) +
  scale_color_manual(
    values = c("Local constant" = "#D55E00", "Local linear" = "#0072B2"),
    breaks = c("Local constant", "Local linear"),
    labels = c("constant", "linear"),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Local constant" = "dashed", "Local linear" = "dotdash"),
    breaks = c("Local constant", "Local linear"),
    labels = c("constant", "linear"),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.4, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.43)) +
  labs(
    title = "C  Integrated mean drift",
    x = expression("Neighborhood radius " * h),
    y = "Mean squared residual"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 9.5),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 8.5),
    legend.key.width = grid::unit(0.9, "lines"),
    plot.margin = margin(4, 4, 2, 6)
  )

grDevices::cairo_pdf(
  "manuscript/figures/local_linear_drift.pdf",
  width = 7.2,
  height = 2.75
)
grid::grid.newpage()
layout <- grid::grid.layout(nrow = 1L, ncol = 3L)
grid::pushViewport(grid::viewport(layout = layout))
print(panel_a, vp = grid::viewport(layout.pos.col = 1L))
print(panel_b, vp = grid::viewport(layout.pos.col = 2L))
print(panel_c, vp = grid::viewport(layout.pos.col = 3L))
grid::popViewport()
grDevices::dev.off()

stopifnot(
  all(drift$contribution >= 0),
  all(
    drift$contribution[drift$method == "Local linear"] <=
      drift$contribution[drift$method == "Local constant"]
  ),
  max(abs(
    drift$contribution[drift$method == "Local constant"] -
      drift$contribution[drift$method == "Local linear"] -
      h_grid^2 / 3
  )) < 1e-12
)

cat("Built the local-linear conditional-mean drift illustration.\n")
