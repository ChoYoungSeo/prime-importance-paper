root <- normalizePath(dirname(sub(
  "^--file=", "", grep("^--file=", commandArgs(), value = TRUE)
)))
invisible(capture.output(source(file.path(root, "run_central_ablation.R"))))

# The allocation records retain Shapley for estimator-sensitivity analysis.
# This timing run replaces only that post-allocation calculation by a constant
# vector so the measured interval ends with projected-game PME, matching the
# scope of the pure Full LL-PME timing audit.
shapley_exact <- function(game, p) rep(1 / p, p)

plan <- normalizePath(file.path(root, ".."))
input_dir <- file.path(
  plan, "execution_20260919/results/4f29d1c54aff9649/checkpoints"
)
inputs <- sort(list.files(
  input_dir,
  pattern = "^friedman_n200_rho0\\.5_r[0-9]{3}\\.input\\.rds$",
  full.names = TRUE
))
stopifnot(length(inputs) == 30L)

# One unrecorded warm-up.
warm_input <- readRDS(inputs[1L])
gc()
warm_started <- proc.time()[["elapsed"]]
warm_result <- matched_k_local_constant_game(warm_input$X, warm_input$y)
warm_elapsed <- proc.time()[["elapsed"]] - warm_started

rows <- lapply(inputs, function(input_path) {
  input <- readRDS(input_path)
  id <- input$task$id
  reference <- readRDS(file.path(
    root, "checkpoints", paste0(id, ".matched_k_lc.rds")
  ))
  gc()
  started <- proc.time()[["elapsed"]]
  result <- matched_k_local_constant_game(input$X, input$y)
  elapsed <- proc.time()[["elapsed"]] - started
  gap <- max(abs(result$pme - reference$details$pme))
  stopifnot(gap < 1e-12)
  data.frame(
    task_id = id,
    replication = as.integer(input$task$replication),
    elapsed_sec = elapsed,
    pme_max_abs_gap = gap,
    stringsAsFactors = FALSE
  )
})
rows <- do.call(rbind, rows)
rows <- rows[order(rows$replication), ]
utils::write.csv(rows, file.path(root, "pure_pme_timing.csv"), row.names = FALSE)

receipt <- list(
  status = "PASS",
  scope = paste(
    "Matched-K full-coordinate local-constant game construction,",
    "monotone projection, and PME on 30 frozen Friedman datasets",
    "at n=200 and rho=0.5"
  ),
  excluded = c(
    "package and source loading", "file input/output", "hashing",
    "pre-run garbage collection", "Shapley allocation"
  ),
  timer = "proc.time()[['elapsed']]",
  schedule = "one unrecorded warm-up, then one serial timed run per dataset",
  r_workers = 1L,
  native_threads = 1L,
  warmup_elapsed_sec = warm_elapsed,
  complete = nrow(rows),
  median_sec = stats::median(rows$elapsed_sec),
  q1_sec = unname(stats::quantile(rows$elapsed_sec, 0.25)),
  q3_sec = unname(stats::quantile(rows$elapsed_sec, 0.75)),
  mean_sec = mean(rows$elapsed_sec),
  min_sec = min(rows$elapsed_sec),
  max_sec = max(rows$elapsed_sec),
  maximum_pme_abs_gap = max(rows$pme_max_abs_gap),
  prime_version = as.character(utils::packageVersion("prime")),
  script_md5 = unname(tools::md5sum(file.path(root, "time_central_pme.R")))
)
jsonlite::write_json(
  receipt, file.path(root, "PURE_PME_TIMING_RECEIPT.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = 16
)
print(receipt)
