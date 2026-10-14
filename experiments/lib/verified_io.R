# Run scripts inside the verified workflow's source/protocol-pinned namespace.
verified_require_context <- function() {
  root <- Sys.getenv("PRIME_VERIFIED_WORK_ROOT", "")
  signature <- Sys.getenv("PRIME_VERIFIED_SIGNATURE", "")
  marker <- file.path(root, "VERIFIED_CONTEXT.dcf")
  if (!nzchar(root) || !nzchar(signature) || !file.exists(marker) ||
      !identical(normalizePath(getwd()), normalizePath(root))) {
    stop("Use the verified Makefile workflow; direct historical-root runs are disabled.",
         call. = FALSE)
  }
  context <- read.dcf(marker)
  if (!identical(unname(context[1L, "Signature"]), signature) ||
      !identical(unname(context[1L, "WorkRoot"]), normalizePath(root))) {
    stop("Verified run context does not match this process.", call. = FALSE)
  }
  invisible(TRUE)
}

# Existing atomic writers pass their temporary pathname and the final path's
# compression flag, so a .gz.tmp-PID file is still valid gzip before rename.
verified_write_csv <- function(x, file, compressed = endsWith(file, ".gz"), ...) {
  verified_require_context()
  connection <- if (compressed) gzfile(file, "wt") else base::file(file, "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.csv(x, connection, ...)
}
