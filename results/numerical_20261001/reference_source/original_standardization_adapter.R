# Shared preprocessing and public-argument orchestration only.
# All estimator source bodies remain frozen. Call standardize_input once per input.
.standardized_state <- new.env(parent = emptyenv())
sha256_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sha256_object <- function(x) digest::digest(x, algo = "sha256")
atomic_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid()); saveRDS(x, tmp, version = 3)
  stopifnot(file.rename(tmp, path)); invisible(path)
}
atomic_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid()); write.csv(x, tmp, row.names = FALSE, na = "NA")
  stopifnot(file.rename(tmp, path)); invisible(path)
}
atomic_json <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  jsonlite::write_json(x, tmp, pretty = TRUE, auto_unbox = TRUE, digits = 17, null = "null", na = "null")
  stopifnot(file.rename(tmp, path)); invisible(path)
}
standardize_input <- function(X) {
  raw_hash <- sha256_object(X)
  if (!(is.matrix(X) || is.data.frame(X)) || !all(vapply(as.data.frame(X), is.numeric, logical(1))))
    stop("standardization_preflight: X must be a numeric matrix or data frame", call. = FALSE)
  X <- as.matrix(X); storage.mode(X) <- "double"
  if (nrow(X) < 2L || ncol(X) < 1L || any(!is.finite(X)))
    stop("standardization_preflight: invalid dimensions or nonfinite X", call. = FALSE)
  center <- colMeans(X)
  scale <- apply(X, 2L, stats::sd)
  bad <- which(!is.finite(scale) | scale <= 0)
  if (length(bad)) stop("standardization_preflight: nonpositive/nonfinite scale in columns ",
                        paste(bad, collapse = ","), call. = FALSE)
  Z <- sweep(sweep(X, 2L, center, "-"), 2L, scale, "/")
  stopifnot(all(is.finite(Z)), identical(dim(Z), dim(X)), identical(dimnames(Z), dimnames(X)))
  list(Z = Z, center = center, scale = scale, raw_hash = raw_hash,
       standardized_hash = sha256_object(Z))
}
init_standardized_runtime <- function(project_or_review) {
  path <- normalizePath(project_or_review, mustWork = TRUE)
  project <- if (dir.exists(file.path(path, "planning/reviews"))) path else normalizePath(file.path(path, "../../.."))
  reviews <- file.path(project, "planning/reviews")
  kit <- file.path(reviews, "2026-09-16_meeting_followup/execution_20260919")
  prior <- file.path(reviews, "2026-09-26_original_pme_local")
  for (pkg in c("digest", "jsonlite", "sensitivity", "FNN", "RANN", "gtools"))
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing existing dependency: ", pkg)
  frozen_path <- file.path(kit, "FROZEN.json")
  stopifnot(identical(sha256_file(frozen_path), "a986bba9f2ec5877f426283e2284096bf391e9003f716a30e9f8265f6e4c6047"))
  frozen <- jsonlite::fromJSON(frozen_path, simplifyVector = FALSE)
  frozen_paths <- file.path(kit, names(frozen$files))
  stopifnot(identical(unname(vapply(frozen_paths, sha256_file, character(1))), unname(unlist(frozen$files))))
  source_path <- file.path(prior, "evidence/original_pme_validation/source/pme.R")
  source_sha <- "67304c8f9126d879c380a5e76ba59aa4b08fce2a6fddf8f70ae3476eeb9b1f1d"
  stopifnot(identical(sha256_file(source_path), source_sha))
  official <- new.env(parent = asNamespace("sensitivity"))
  sys.source(source_path, envir = official, keep.source = FALSE)
  fn_hash <- function(env) {
    nms <- ls(env, all.names = TRUE)
    nms <- nms[vapply(nms, function(n) is.function(env[[n]]), logical(1))]
    setNames(vapply(nms, function(n) sha256_object(list(formals(env[[n]]), body(env[[n]]))), character(1)), nms)
  }
  stopifnot(all(vapply(ls(official), function(n) identical(environment(official[[n]]), official), logical(1))))
  official_hashes <- fn_hash(official)
  lockEnvironment(official, bindings = TRUE)
  adapter <- new.env(parent = baseenv())
  sys.source(file.path(kit, "first_adapter.R"), envir = adapter, keep.source = FALSE)
  first <- adapter$load_first_reference(file.path(kit, "vendor/first.R"))
  first_hashes <- fn_hash(first)
  lockEnvironment(first, bindings = TRUE)
  lockEnvironment(adapter, bindings = TRUE)
  .libPaths(c(file.path(kit, "library"), .libPaths()))
  stopifnot(normalizePath(find.package("prime")) == normalizePath(file.path(kit, "library/prime")),
            as.character(utils::packageVersion("prime")) == "0.5.2.9101")
  ns <- asNamespace("prime")
  stopifnot(normalizePath(getNamespaceInfo(ns, "path")) == normalizePath(file.path(kit, "library/prime")))
  prime_hashes <- fn_hash(ns)
  original_paths <- c(frozen_path, frozen_paths, source_path)
  original_hashes <- setNames(vapply(original_paths, sha256_file, character(1)), original_paths)
  assert_frozen <- function() {
    stopifnot(identical(vapply(original_paths, sha256_file, character(1)), original_hashes),
              identical(fn_hash(official), official_hashes), identical(fn_hash(first), first_hashes),
              identical(fn_hash(ns), prime_hashes))
    invisible(TRUE)
  }
  settings <- list(method = "knn", n.knn = 6, noise = TRUE, rescale = FALSE,
                   marg = TRUE, tol = NULL, n.limit = 2001, nboot = NULL, parl = NULL)
  runtime <- list(project = project, kit = kit, official = official, first = first,
                  source_path = source_path, source_sha256 = source_sha,
                  first_source_sha256 = first$source_sha256,
                  frozen_signature = frozen$signature, official_settings = settings,
                  assert_frozen = assert_frozen, prime_function_hashes = prime_hashes,
                  first_function_hashes = first_hashes, official_function_hashes = official_hashes)
  .standardized_state$runtime <- runtime
  runtime
}
.standardized_runtime <- function() {
  if (!exists("runtime", .standardized_state, inherits = FALSE)) stop("Call init_standardized_runtime first")
  .standardized_state$runtime
}
fit_prime_standardized <- function(Z, y, runtime = .standardized_runtime()) {
  force(runtime)
  prime::prime(Z, y, rescale = FALSE)
}
fit_pme_standardized <- function(Z, y, runtime = .standardized_runtime()) {
  runtime$official$pme_knn(model = function(X) y, X = Z, method = "knn", n.knn = 6,
                          noise = TRUE, rescale = FALSE, marg = TRUE, tol = NULL,
                          n.limit = 2001, nboot = NULL, parl = NULL)
}
fit_first_standardized <- function(Z, y, runtime = .standardized_runtime()) {
  run_first_standardized(runtime$first, Z, y)
}

# Exact frozen run_first_shap orchestration, with its two public rescale flags FALSE.
run_first_standardized <- function(env, X, y) {
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
    rescale = FALSE, n.forward = 2, parl = NULL, verbose = FALSE)
  selected <- screen$selected
  if (!length(selected)) unavailable("empty_selection")
  if (anyNA(selected) || any(selected != as.integer(selected)) ||
      any(selected < 1L | selected > p) || anyDuplicated(selected))
    unavailable("invalid_selected_set")
  selected <- as.integer(selected)
  allocation <- env$.shap_capture(
    X = X[, selected, drop = FALSE], y = y, noise = TRUE, n.knn = 2,
    n.mc = nrow(X), twin.mc = FALSE, rescale = FALSE, parl = NULL)
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
