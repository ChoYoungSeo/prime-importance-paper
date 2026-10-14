source("experiments/lib/verified_io.R")
verified_require_context()
library(prime)
source("experiments/10_concrete_application/00_common.R")

profile <- Sys.getenv("PRIME_CONCRETE_PROFILE", "smoke")
configuration <- switch(
  profile,
  smoke = list(bootstraps = 4L, cv_repeats = 1L, cv_folds = 3L),
  main = list(bootstraps = 200L, cv_repeats = 5L, cv_folds = 5L),
  stop("Unknown concrete profile: ", profile, call. = FALSE)
)
output_directory <- file.path(
  "results/raw", paste0("concrete_v052_", profile)
)
checkpoint_directory <- file.path(output_directory, "bootstrap_checkpoints")
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)

data <- concrete_representations()
methods <- list(
  raw = "PRIME",
  engineering = c("PRIME", "Full LL-PME", "PME-kNN")
)
Sys.setenv(PRIME_NUM_THREADS = "1")

run_bootstrap <- function(replication) {
  checkpoint <- file.path(
    checkpoint_directory, sprintf("bootstrap_%04d.rds", replication)
  )
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  seed <- 202608300L + replication
  set.seed(seed)
  rows <- sample.int(length(data$y), replace = TRUE)
  run_rows <- list()
  phi_rows <- list()
  index <- 0L
  for (representation in names(methods)) {
    X <- data[[representation]][rows, , drop = FALSE]
    y <- data$y[rows]
    for (method in methods[[representation]]) {
      index <- index + 1L
      error <- ""
      start <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        concrete_method(method, X, y),
        error = function(condition) {
          error <<- conditionMessage(condition)
          list(
            effects = rep(NA_real_, ncol(X)), selected = integer(),
            value_calls = NA_integer_
          )
        }
      )
      elapsed <- proc.time()[["elapsed"]] - start
      feasible <- all(is.finite(fit$effects))
      run_rows[[index]] <- data.frame(
        replication = replication,
        seed = seed,
        representation = representation,
        method = method,
        selected = paste(fit$selected, collapse = ","),
        selected_size = length(fit$selected),
        value_calls = fit$value_calls,
        feasible = feasible,
        elapsed_sec = elapsed,
        error = error,
        stringsAsFactors = FALSE
      )
      phi_rows[[index]] <- data.frame(
        replication = replication,
        seed = seed,
        representation = representation,
        method = method,
        variable = seq_len(ncol(X)),
        variable_name = colnames(X),
        selected = seq_len(ncol(X)) %in% fit$selected,
        estimate = fit$effects,
        stringsAsFactors = FALSE
      )
    }
  }
  object <- list(
    runs = do.call(rbind, run_rows),
    phi = do.call(rbind, phi_rows)
  )
  concrete_atomic_save(object, checkpoint)
  object
}

cores <- suppressWarnings(as.integer(Sys.getenv("PRIME_CORES", "6")))
if (!length(cores) || !is.finite(cores) || cores < 1L) cores <- 1L
bootstraps <- seq_len(configuration$bootstraps)
objects <- if (.Platform$OS.type == "unix" && cores > 1L) {
  parallel::mclapply(
    bootstraps, run_bootstrap, mc.cores = min(cores, 6L),
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
} else {
  lapply(bootstraps, run_bootstrap)
}
bootstrap_runs <- do.call(rbind, lapply(objects, `[[`, "runs"))
bootstrap_phi <- do.call(rbind, lapply(objects, `[[`, "phi"))
write.csv(
  bootstrap_runs, file.path(output_directory, "bootstrap_runs.csv"),
  row.names = FALSE
)
write.csv(
  bootstrap_phi, file.path(output_directory, "bootstrap_phi.csv"),
  row.names = FALSE
)

X <- data$engineering
y <- data$y
n <- nrow(X)
k_grid <- c(5L, 10L, 20L, 30L, 40L)
cv_rows <- list()
cv_index <- 0L
for (repeat_index in seq_len(configuration$cv_repeats)) {
  set.seed(202609000L + repeat_index)
  outer_fold <- sample(rep(seq_len(configuration$cv_folds), length.out = n))
  for (outer in seq_len(configuration$cv_folds)) {
    test <- which(outer_fold == outer)
    train <- which(outer_fold != outer)
    screen <- prime_screen(X[train, , drop = FALSE], y[train])
    selected <- screen$selected
    if (!length(selected)) selected <- seq_len(ncol(X))
    sets <- list(
      `All coordinates` = seq_len(ncol(X)),
      `Selected coordinates` = selected,
      `Ratio and age` = match(c("water_binder_ratio", "age"), colnames(X))
    )
    set.seed(202610000L + 100L * repeat_index + outer)
    inner_fold <- sample(rep(seq_len(5L), length.out = length(train)))
    for (set_name in names(sets)) {
      variables <- sets[[set_name]]
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
      cv_index <- cv_index + 1L
      cv_rows[[cv_index]] <- data.frame(
        repeat_id = repeat_index,
        fold = outer,
        feature_set = set_name,
        selected = paste(selected, collapse = ","),
        selected_size = length(selected),
        k = tuned_k,
        rmse = sqrt(mean((y[test] - prediction)^2)),
        stringsAsFactors = FALSE
      )
    }
  }
}
cv <- do.call(rbind, cv_rows)
write.csv(cv, file.path(output_directory, "nested_cv.csv"), row.names = FALSE)
cat(
  "Completed", configuration$bootstraps,
  "bootstrap samples and", nrow(cv), "nested-CV rows.\n"
)
