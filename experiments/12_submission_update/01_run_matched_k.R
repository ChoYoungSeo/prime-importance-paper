args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1] else
  "results/generated/matched_k_submission_20260920"
workers <- as.integer(Sys.getenv("PRIME_CORES", "2"))
if (!is.finite(workers) || workers < 1L) workers <- 1L

suppressPackageStartupMessages(library(prime))
source("experiments/00_reference/benchmarks.R")
source("experiments/03_pme/00_common.R")
source("experiments/12_submission_update/matched_k_helpers.R")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

population <- read.csv(
  "results/submission_20260920/friedman/population_allocations.csv",
  stringsAsFactors = FALSE
)
grid <- expand.grid(
  n = c(100L, 200L, 500L),
  rho = c(0, 0.5, 0.9),
  replication = seq_len(30L),
  KEEP.OUT.ATTRS = FALSE
)
grid <- grid[order(grid$n, grid$rho, grid$replication), ]
n_index <- match(grid$n, c(100L, 200L, 500L))
rho_index <- match(grid$rho, c(0, 0.5, 0.9))
grid$seed <- 202609160L + 10000L * n_index + 1000L * rho_index +
  grid$replication

run_one <- function(task) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(task$seed)
  z <- matrix(rnorm(task$n * 5L), ncol = 5L)
  e <- matrix(rnorm(task$n * 5L), ncol = 5L)
  X <- cbind(
    pnorm(z),
    pnorm(task$rho * z + sqrt(1 - task$rho^2) * e)
  )
  colnames(X) <- paste0("X", seq_len(10L))
  y <- friedman_function(X) + rnorm(task$n, sd = 1)
  fit <- matched_k_local_constant_game(X, y)
  truth <- population[abs(population$rho - task$rho) < 1e-12, ]
  stopifnot(
    nrow(truth) == 10L,
    identical(as.character(truth$variable), colnames(X)),
    length(fit$pme) == ncol(X),
    length(fit$shapley) == ncol(X)
  )
  data.frame(
    task_id = sprintf(
      "friedman_n%d_rho%g_r%03d",
      task$n, task$rho, task$replication
    ),
    n = task$n,
    rho = task$rho,
    replication = task$replication,
    seed = task$seed,
    pme_half_l1 = 0.5 * sum(abs(fit$pme - truth$PME)),
    shapley_half_l1 =
      0.5 * sum(abs(fit$shapley - truth$Shapley_full10)),
    pme_proxy_mass = sum(fit$pme[6:10]),
    shapley_proxy_mass = sum(fit$shapley[6:10]),
    stringsAsFactors = FALSE
  )
}

tasks <- split(grid, seq_len(nrow(grid)))
rows <- parallel::mclapply(
  tasks,
  function(task) run_one(task[1, ]),
  mc.cores = workers,
  mc.preschedule = FALSE
)
summary <- do.call(rbind, rows)
write.csv(summary, file.path(output, "full_summary.csv"), row.names = FALSE)
metrics <- c(
  "pme_half_l1", "shapley_half_l1",
  "pme_proxy_mass", "shapley_proxy_mass"
)
cells <- aggregate(
  summary[, metrics],
  by = list(n = summary$n, rho = summary$rho),
  FUN = mean
)
write.csv(cells, file.path(output, "cell_means.csv"), row.names = FALSE)

reference <- read.csv(
  "results/submission_20260920/matched_k/cell_means.csv",
  stringsAsFactors = FALSE
)
columns <- c("n", "rho", metrics)
stopifnot(
  identical(cells[, c("n", "rho")], reference[, c("n", "rho")]),
  max(abs(as.matrix(cells[, metrics]) -
          as.matrix(reference[, metrics]))) < 1e-12
)
cat("Matched-K Friedman ablation: PASS (270/270 datasets)\n")
