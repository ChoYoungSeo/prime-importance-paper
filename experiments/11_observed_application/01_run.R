source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/11_observed_application/00_common.R")

profile <- Sys.getenv("PRIME_OBSERVED_APPLICATION_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(bootstraps = 4L, cv_repeats = 1L, cv_folds = 3L),
  main = list(bootstraps = 200L, cv_repeats = 5L, cv_folds = 5L),
  stop("Unknown observed-application profile: ", profile, call. = FALSE)
)
datasets <- c("airfoil", "auto_mpg")
methods <- c("PRIME", "Full LL-PME", "PME-kNN")
output_directory <- file.path(
  "results/raw", paste0("observed_application_v052_", profile)
)
bootstrap_directory <- file.path(output_directory, "bootstrap_checkpoints")
cv_directory <- file.path(output_directory, "cv_checkpoints")
dir.create(bootstrap_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(cv_directory, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PRIME_NUM_THREADS = "1")

bootstrap_tasks <- expand.grid(
  dataset = datasets,
  replication = seq_len(configuration$bootstraps),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
bootstrap_tasks$task_id <- seq_len(nrow(bootstrap_tasks))

run_bootstrap <- function(task_id) {
  task <- bootstrap_tasks[task_id, , drop = FALSE]
  checkpoint <- file.path(
    bootstrap_directory, sprintf("bootstrap_%04d.rds", task_id)
  )
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  dataset_index <- match(task$dataset, datasets)
  seed <- 202611000L + 10000L * dataset_index + task$replication
  data <- load_observed_data(task$dataset)
  set.seed(seed)
  rows <- sample.int(length(data$y), replace = TRUE)
  X <- data$X[rows, , drop = FALSE]
  y <- data$y[rows]
  run_rows <- vector("list", length(methods))
  phi_rows <- vector("list", length(methods))
  for (index in seq_along(methods)) {
    method <- methods[index]
    error <- ""
    start <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      observed_application_method(method, X, y),
      error = function(condition) {
        error <<- conditionMessage(condition)
        list(
          effects = rep(NA_real_, ncol(X)), selected = integer(),
          value_calls = NA_integer_
        )
      }
    )
    elapsed <- proc.time()[["elapsed"]] - start
    selected <- if (is.null(fit$selected)) integer() else fit$selected
    feasible <- all(is.finite(fit$effects))
    run_rows[[index]] <- data.frame(
      task_id = task_id,
      dataset = task$dataset,
      dataset_label = data$label,
      replication = task$replication,
      seed = seed,
      n = nrow(X),
      p = ncol(X),
      method = method,
      selected = paste(selected, collapse = ","),
      selected_size = length(selected),
      value_calls = fit$value_calls,
      feasible = feasible,
      elapsed_sec = elapsed,
      error = error,
      stringsAsFactors = FALSE
    )
    phi_rows[[index]] <- data.frame(
      task_id = task_id,
      dataset = task$dataset,
      dataset_label = data$label,
      replication = task$replication,
      seed = seed,
      method = method,
      variable = seq_len(ncol(X)),
      variable_name = colnames(X),
      selected = seq_len(ncol(X)) %in% selected,
      estimate = fit$effects,
      stringsAsFactors = FALSE
    )
  }
  object <- list(
    runs = do.call(rbind, run_rows),
    phi = do.call(rbind, phi_rows)
  )
  observed_application_atomic_save(object, checkpoint)
  object
}

cores <- suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
if (!length(cores) || !is.finite(cores) || cores < 1L) cores <- 1L
bootstrap_objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    bootstrap_tasks$task_id, run_bootstrap,
    mc.cores = min(cores, 6L), mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(bootstrap_tasks$task_id, run_bootstrap)
}
bootstrap_runs <- do.call(rbind, lapply(bootstrap_objects, `[[`, "runs"))
bootstrap_phi <- do.call(rbind, lapply(bootstrap_objects, `[[`, "phi"))
bootstrap_runs <- bootstrap_runs[order(
  bootstrap_runs$dataset, bootstrap_runs$replication, bootstrap_runs$method
), , drop = FALSE]
bootstrap_phi <- bootstrap_phi[order(
  bootstrap_phi$dataset, bootstrap_phi$replication,
  bootstrap_phi$method, bootstrap_phi$variable
), , drop = FALSE]
write.csv(
  bootstrap_runs, file.path(output_directory, "bootstrap_runs.csv"),
  row.names = FALSE
)
write.csv(
  bootstrap_phi, file.path(output_directory, "bootstrap_phi.csv"),
  row.names = FALSE
)

cv_tasks <- expand.grid(
  dataset = datasets,
  repeat_id = seq_len(configuration$cv_repeats),
  fold = seq_len(configuration$cv_folds),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
cv_tasks$task_id <- seq_len(nrow(cv_tasks))
k_grid <- c(5L, 10L, 20L, 30L, 40L)

run_cv <- function(task_id) {
  task <- cv_tasks[task_id, , drop = FALSE]
  checkpoint <- file.path(cv_directory, sprintf("cv_%04d.rds", task_id))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  data <- load_observed_data(task$dataset)
  X <- data$X
  y <- data$y
  dataset_index <- match(task$dataset, datasets)
  set.seed(202612000L + 1000L * dataset_index + task$repeat_id)
  outer_fold <- sample(rep(
    seq_len(configuration$cv_folds), length.out = nrow(X)
  ))
  test <- which(outer_fold == task$fold)
  train <- which(outer_fold != task$fold)
  screen <- prime_screen(X[train, , drop = FALSE], y[train])
  selected <- screen$selected
  if (!length(selected)) selected <- seq_len(ncol(X))
  feature_sets <- list(
    `All coordinates` = seq_len(ncol(X)),
    `Selected coordinates` = selected
  )
  set.seed(
    202613000L + 10000L * dataset_index +
      100L * task$repeat_id + task$fold
  )
  inner_fold <- sample(rep(seq_len(5L), length.out = length(train)))
  rows <- lapply(names(feature_sets), function(feature_set) {
    variables <- feature_sets[[feature_set]]
    tuned_k <- concrete_tune_k(
      X[train, , drop = FALSE], y[train], variables, inner_fold, k_grid
    )
    scaled <- concrete_standardize(
      X[train, variables, drop = FALSE],
      X[test, variables, drop = FALSE]
    )
    prediction <- concrete_knn_predict(
      scaled$train, y[train], scaled$test, tuned_k
    )
    data.frame(
      task_id = task_id,
      dataset = task$dataset,
      dataset_label = data$label,
      repeat_id = task$repeat_id,
      fold = task$fold,
      feature_set = feature_set,
      selected = paste(selected, collapse = ","),
      selected_size = length(selected),
      k = tuned_k,
      rmse = sqrt(mean((y[test] - prediction)^2)),
      stringsAsFactors = FALSE
    )
  })
  object <- do.call(rbind, rows)
  observed_application_atomic_save(object, checkpoint)
  object
}

cv_objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    cv_tasks$task_id, run_cv,
    mc.cores = min(cores, 6L), mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(cv_tasks$task_id, run_cv)
}
cv <- do.call(rbind, cv_objects)
cv <- cv[order(cv$dataset, cv$repeat_id, cv$fold, cv$feature_set), ]
write.csv(cv, file.path(output_directory, "nested_cv.csv"), row.names = FALSE)

cat(
  "Completed", nrow(bootstrap_tasks), "bootstrap tasks and",
  nrow(cv_tasks), "nested-CV tasks.\n"
)
