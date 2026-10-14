# Validated reuse of paired outputs

The benchmark generators, predictor laws, task seeds, and local-constant versus
local-linear conditional-variance diagnostic are inherited from the
professor-coauthored project. PME reference objects are recomputed from stored
coalition values or analytic moments.

The BIC-gated ridge start and paired querywise calibration change only the
PRIME selector and its downstream allocation. Therefore:

- every PRIME coordinate-recovery row is rerun;
- every downstream row supplied by the PRIME selected set is rerun;
- FIRST, LASSO, FOCI, and FORD rows are reused only when task keys and seeds
  match exactly;
- full local-linear PME, PME-kNN, NCMD-PME, and oracle-coordinate controls are
  reused only when their data, target, and seed are unchanged.

The Ackley-3 competitor rows are paired reruns from the locked replacement
benchmark. The correlated and ambient-scaling LASSO, FOCI, and FORD pipelines
are paired baselines on the same task seeds.

Main-profile assembly checks all expected keys and seeds before combining
reused rows and records `result_origin`.

The matched-neighborhood conditional-variance summary is unaffected by the
selector. Observed-data and Data-Twinning diagnostics were rerun under the
validated `prime 0.5.2` implementation. PRIME rows in the current summaries
come from these reruns; reused competitor rows retain the same task keys and
seeds.
