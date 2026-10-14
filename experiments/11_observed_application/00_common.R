source("experiments/10_concrete_application/00_common.R")

observed_application_method <- function(method, X, y) {
  switch(
    method,
    "PRIME" = {
      fit <- prime::prime(X, y)
      if (!length(fit$selected)) {
        stop("empty_screen: allocation unavailable", call. = FALSE)
      }
      list(
        effects = fit$effects,
        selected = fit$selected,
        value_calls = if (is.null(fit$subset_values)) 0L else nrow(fit$subset_values)
      )
    },
    "Full LL-PME" = estimate_pme_with_value(
      make_llr_value(X, y), seq_len(ncol(X)), ncol(X)
    ),
    "PME-kNN" = estimate_pme_with_value(
      make_local_constant_value(X, y, neighbors = 2L),
      seq_len(ncol(X)), ncol(X)
    ),
    stop("Unknown observed-data method: ", method, call. = FALSE)
  )
}

observed_application_atomic_save <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save observed-data checkpoint: ", path, call. = FALSE)
  }
}
