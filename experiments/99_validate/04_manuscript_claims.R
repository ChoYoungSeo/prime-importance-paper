# Validate live values and complete vectors; historical constants are not a gate.
source("experiments/lib/verified_io.R")
verified_require_context()
canonical <- Sys.getenv("PRIME_CANONICAL_ROOT", "")
if (!nzchar(canonical)) stop("Missing verified canonical-root context.")
context <- jsonlite::fromJSON("VERIFIED_CONTEXT.json", simplifyVector = FALSE)
publication_check <- identical(context$pin$namespace, "release-validation")
output <- file.path(getwd(), "results/validation/numeric_claims")
arguments <- c(
  "experiments/99_validate/numeric_claims.py", "--data-root", getwd(),
  "--baseline-root", file.path(getwd(), "verification_baseline"),
  "--fixed-root", file.path(getwd(), "verification_baseline"),
  "--candidate", file.path(canonical, "manuscript"), "--out", output
)
if (publication_check) arguments <- c(arguments, "--validate-prose",
  "--finalized-inputs", "--observed-current", "--timings-current")
status <- system2("python3", vapply(arguments, shQuote, character(1)))
if (status != 0L) stop("Vector/data/prose validation failed; see ", output)
if (publication_check) {
  arguments <- c("experiments/99_validate/validate_numeric_tables.py",
    "--snapshot", file.path(output, "numeric_claim_snapshot.json"),
    "--candidate", file.path(canonical, "manuscript"),
    "--out", file.path(output, "table_validation.json"))
  status <- system2("python3", vapply(arguments, shQuote, character(1)))
  if (status != 0L) stop("Published table entries do not match verified inputs.")
}
cat(if (publication_check) "Verified-release data, prose, and tables: PASS\n" else
  "Run-data/vector integrity: PASS; manuscript alignment requires make validate-release.\n")
