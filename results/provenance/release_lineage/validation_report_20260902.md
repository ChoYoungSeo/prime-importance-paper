# Validation

The release was validated on 2026-09-02 under R 4.6.1 on Apple silicon.

- `R CMD check --no-manual` completed with `Status: OK`.
- All package unit, regression, and public-API tests passed.
- The embedded package and canonical development implementation
  agree exactly on selected sets and to numerical precision on forward and
  backward diagnostics, subset values, monotone projection, and allocation.
- Direct R and compiled-batch local-linear calculations agree on every parity
  case.
- PME-kNN subset values match the recorded reference implementation.
- Every main-profile validator passes expected-row, unique-key, seed,
  method-set, allocation-sum, feasibility, and reference-sum checks.
- The 36 aggregate result files remain numerically identical under the PRIME
  label migration.

The numerical-claim audit confirms:

- selector-to-PME performance in 6/6 cells;
- low-dimensional PME performance in 18/18 cells;
- correlated-pipeline performance in 6/6 cells;
- ambient-dimension performance in 25/25 cells;
- local-linear subset-value improvement in 28/28 matched cells;
- 600 observed-data bootstrap samples and 75 outer cross-validation folds.

The validated profiles comprise 2,250 coordinate-selection datasets, 900
selector-to-PME rows, 540 low-dimensional PME datasets, 720 correlated
pipeline rows, 4,500 ambient-scaling rows, 240 Data-Twinning runs, and 600
observed-data bootstrap samples. Detailed reference and lineage records are
stored in `results/provenance/`.
