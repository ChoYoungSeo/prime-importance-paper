# Isolated reporting adapter for the unmodified official FIRST supplementary source.
# The caller owns the RNG seed. No installed first package is loaded.
load_first_reference <- function(path) {
  expected <- "759c110ef16e4cce1fcebf17e1d1ed03035a90c92878e02715ff038f1afe1ee4"
  if (!file.exists(path)) stop("FIRST reference source is unavailable: ", path)
  for (pkg in c("FNN", "digest")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("FIRST reference dependency is unavailable: ", pkg)
  }
  sha <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(sha, expected)) stop("FIRST reference SHA-256 mismatch")
  env <- new.env(parent = baseenv())
  # The reference uses unqualified runif only for its documented tie jitter.
  env$runif <- stats::runif
  sys.source(path, envir = env)
  capture_return <- function(fun, original, replacement) {
    b <- body(fun)
    if (!identical(b[[length(b)]], original))
      stop("FIRST reference final return differs from the pinned contract")
    b[[length(b)]] <- replacement
    body(fun) <- b
    fun
  }
  env$.first_capture <- capture_return(
    env$first, quote(return(imp)),
    quote(return(list(selected = subset, total_importance = imp))))
  env$.shap_capture <- capture_return(
    env$shapleysobol.knn, quote(return(ssi)),
    quote(return(list(
      effects = ssi, denominator = y.var - noise.var, risk_full = noise.var,
      game = if (exists("x.var", inherits = FALSE)) x.var else NULL))))
  env$source_sha256 <- sha
  env$source_path <- normalizePath(path, mustWork = TRUE)
  env
}

run_first_shap <- function(env, X, y) {
  unavailable <- function(reason) stop("FIRST-SHAP unavailable: ", reason,
                                     call. = FALSE)
  if (!is.environment(env) || !is.function(env$.first_capture) ||
      !is.function(env$.shap_capture)) unavailable("invalid_reference_environment")
  if (!(is.matrix(X) || is.data.frame(X)) || ncol(X) < 1L || nrow(X) < 2L)
    unavailable("invalid_input_dimensions")
  if (is.matrix(y)) y <- c(y)
  if (!is.numeric(y) || length(y) != nrow(X) || any(!is.finite(y)))
    unavailable("invalid_response")
  p <- ncol(X)
  screen <- env$.first_capture(
    X = X, y = y, n.knn = 2, n.mc = nrow(X), twin.mc = FALSE,
    rescale = TRUE, n.forward = 2, parl = NULL, verbose = FALSE)
  selected <- screen$selected
  if (!length(selected)) unavailable("empty_selection")
  if (anyNA(selected) || any(selected != as.integer(selected)) ||
      any(selected < 1L | selected > p) || anyDuplicated(selected))
    unavailable("invalid_selected_set")
  selected <- as.integer(selected)
  allocation <- env$.shap_capture(
    X = X[, selected, drop = FALSE], y = y, noise = TRUE, n.knn = 2,
    n.mc = nrow(X), twin.mc = FALSE, rescale = TRUE, parl = NULL)
  denominator <- allocation$denominator
  if (length(denominator) != 1L || !is.finite(denominator) || denominator <= 0)
    unavailable("nonpositive_or_nonfinite_denominator")
  if (length(allocation$risk_full) != 1L || !is.finite(allocation$risk_full) ||
      length(allocation$effects) != length(selected) ||
      any(!is.finite(allocation$effects))) unavailable("nonfinite_allocation")
  if (abs(sum(allocation$effects) - 1) > 1e-8)
    unavailable("shapley_efficiency_failure")
  if (!is.null(allocation$game) && any(!is.finite(allocation$game)))
    unavailable("nonfinite_game")
  effects <- numeric(p)
  effects[selected] <- allocation$effects  # Signed estimates are retained.
  list(effects = effects, selected = selected, denominator = denominator,
       risk_full = allocation$risk_full, game = allocation$game,
       total_importance = screen$total_importance)
}
