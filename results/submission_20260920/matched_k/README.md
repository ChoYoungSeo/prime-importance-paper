# Matched-\(K\) local-constant allocation ablation

This experiment evaluates whether the improvement of the full-coordinate
local-linear estimator can be attributed only to its larger neighborhoods.

For every nonempty Friedman subset \(U\), the ablation uses the same
\(K_n(|U|)=\lceil\max\{1.5\sqrt n,3|U|\}\rceil\) rule as PRIME and the full
local-linear estimator. It selects the exact \(K_n(|U|)\) nearest other
observations in the projected coordinates, computes their sample response
variance with divisor \(K_n(|U|)-1\), uses the full-set risk as the noise floor,
and applies the same positive-part normalization, monotone projection, PME, and
Shapley calculations. The only intended change relative to Full LL is the
local fit order.

The run reuses the 270 frozen Friedman datasets: 30 replications for every
\(n\in\{100,200,500\}\) and \(\rho\in\{0,0.5,0.9\}\). All 270 allocations
completed.

Files:

- run_full_ablation.R: complete resumable computation.
- checkpoints/: per-dataset JSON and RDS records.
- full_summary.csv: replication-level error, proxy mass, and timing.
- cell_means.csv: nine-cell averages.
- FULL_RECEIPT.json: scope, \(K\) values, completion, and aggregate checks.
- run_central_ablation.R, central_summary.csv, and RECEIPT.json: serial
  central-condition timing run used in the computation table.
- time_central_pme.R, pure_pme_timing.csv, and PURE_PME_TIMING_RECEIPT.json:
  scope-matched central timing that excludes Shapley allocation.

Full LL-PME has lower mean PME error than matched-\(K\) local-constant PME in
all nine cells. In the central condition, their mean errors are 0.09379 and
0.44029. The matched-\(K\) result is also worse than the 2NN baseline, showing
that increasing a local-constant neighborhood does not reproduce the
local-linear improvement.
