# Validation

## PyPI final release candidate: 0.2.1 — 2026-10-12

This local candidate updates `0.2.1.dev1` to the final version `0.2.1` in package
metadata and the `__version__` constant. Estimator logic, API aliases, tests,
examples, figures, and benchmark results remain unchanged. The PyPI description
uses ordinary installation and has no local/private/unpublished review text.

The wheel and source distribution rebuild successfully and pass strict Twine.
A fresh installation in a separate directory reports `0.2.1` and passes both
README examples and the standalone quickstart. The prior 91-test result below
is carried forward on unchanged implementation/tests; that suite and the CI
matrix are not repeated for this metadata-only release preparation. No package
has been uploaded and no account or publisher configuration has been changed.

## Standalone repository review basis — 2026-10-12

This basis retains the existing Python repository
`https://github.com/ChoYoungSeo/prime-importance`. Installation uses its root;
the paper and R package remain in their existing repositories. Youngseo Cho's
current public contact is `josun@gatech.edu`.

Source, tests, examples, figures, and benchmark results are byte-identical to the
19 protected files in the 91-test validation basis below. This standalone review
rebuilds the source distribution and wheel, checks both with strict Twine, and
installs the wheel in a separate target directory to run both README examples
and the standalone quickstart. These packaging and example checks pass; the
91-test suite is carried forward from the matching code and is not rerun here.

The checks use cached tools and numerical dependencies on macOS arm64 with
Python 3.13.3, NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2. No repository
update or PyPI publication has occurred. The operating-system/Python-version
matrix is not repeated.

## Historical consolidated-layout validation — 2026-10-12

The repository URL and installation examples now use
`https://github.com/ChoYoungSeo/prime-importance-paper`, with the Python package
in `packages/python`. Youngseo Cho's current contact is `josun@gatech.edu`.
The distribution remains `prime-importance` version `0.2.1.dev1`, and the import
namespace remains `prime`. Implementation, tests, examples, figures, and benchmark
results are byte-identical to the earlier review basis.

The local source distribution and wheel build without network access; the wheel
is built from the source distribution. Strict Twine metadata checks pass. A new
installation into a separate target directory passes both README examples, the
standalone quickstart, and all **91 existing tests**. The README examples select
`X1`, `X2`, and `X5`, with zero assigned to the remaining features.

These checks use cached build tools and numerical dependencies on macOS arm64
with Python 3.13.3, NumPy 2.5.3, SciPy 1.18.1, and Matplotlib 3.11.2. They do not
repeat the operating-system/Python-version matrix. The repository layout has
not been published, and nothing has been uploaded to PyPI.

## Historical candidate: 0.2.1.dev1 — 2026-10-11

All **91 tests pass** against a freshly built wheel installed in a separate
virtual environment on macOS arm64 with Python 3.13.3, NumPy 2.5.4, and
SciPy 1.18.1. Both README examples and the standalone quickstart pass.
The isolated sdist-to-wheel build, Ruff, strict Twine checks, and `pip check`
also pass.

This version makes BIC and paired calibration invariant to response units,
aligns ridge predictor scaling with the R implementation, and uses normalized
games and log-domain potentials for stable proportional allocation. It also
validates query indices and feature names before conversion and supports
NumPy's copy-aware array protocol. Regression tests cover these fixes; seven
shared R/Python allocation fixtures agree to a maximum relative difference of
`2.13e-16`, with a separate R ridge-reference regression test.

The sections below are historical validation records for earlier candidates.

## Historical local submission candidate: 2026-10-01

The candidate version remains `0.2.1.dev0`. Its final allocation preserves small
positive effects: `zero_tolerance`
continues to classify projected game values, but does not truncate the final
allocation vector. A regression test checks efficiency at the default tolerance
and at `0.1`, and verifies that genuinely null game values still yield zero.

The 36 Python tests pass on macOS arm64 with Python 3.13.3. Ruff passes, an
isolated sdist-to-wheel build succeeds, and Twine's strict metadata check passes.
An independent virtual environment installs the wheel and runs all 36 tests and
the complete quickstart. These local checks do not substitute for rerunning the
GitHub operating-system/Python-version matrix on this new candidate.

The README's 2026-09-13 benchmark table below is historical evidence and was not
rerun during this packaging audit. No package was uploaded.

## Neighbor backend update: 2026-09-13

Development version `0.2.1.dev0` adds streamed projected-tree queries,
scalar subset-risk caching, bounded querywise-detail caching, and explicit
`neighbor_eps`/`query_block_size` arguments. The exact default retains all
observations and the existing neighborhood and calibration rules.

- All **35 tests pass**. New tests compare streamed exact neighborhoods and
  querywise risks with the prior dense calculation, including block sizes
  1/64/256/1024, reversed query subsets, duplicate rows, and equidistant ties.
  Exact contribution arrays, screening traces, subset values, and effects
  agree exactly on these fixtures.
- Tests verify one tree per subset, bounded query batches, the complete
  neighbor pool, bounded detail caching, scalar-only allocation caching,
  the approximate neighbor-distance bound, and public argument validation.
- Source, tests, examples, and benchmark runner pass Ruff; `git diff --check`
  passes. The source distribution and wheel build successfully.

The paired benchmark uses one fresh process per case/backend, fixed seed
20260913, all queries, and one BLAS thread. Wall time excludes imports and
input generation; peak RSS covers the entire worker lifetime. Each cell is a
single development measurement, not a timing confidence interval.

| Case | Prior exact seconds / MiB | Streamed exact seconds / MiB | eps=0.2 seconds / MiB |
|---|---:|---:|---:|
| One subset, n=15,000, d=6 | 1.060 / 129.1 | 1.015 / 69.5 | 0.876 / 70.7 |
| All 31 nonempty subsets, n=3,000, p=5 | 2.063 / 90.8 | 1.846 / 65.5 | 1.803 / 65.3 |
| Complete fit, n=1,500, p=10 | 0.622 / 69.9 | 0.643 / 63.1 | 0.595 / 63.6 |

All exact risks and effects equal the prior results in the recorded cases.
For eps=0.2, the complete-fit selected set remains `(0, 1, 2)` and the largest
absolute effect change is `0.0003612173`. Across the 31 subset risks, the largest
absolute ANN change is `0.0037245094`. These examples do not establish an
accuracy guarantee or universal speedup for approximate search. Exact
streaming lowers peak memory in each measured case; the full-fit elapsed time
is slightly higher in this one measurement.

The subset-cache case retains 744,000 contribution bytes in the prior code
and zero in the new scalar cache. Temporary querywise vectors still require
O(n_queries) memory, and exact PME enumeration remains exponential in the
selected dimension.

Reproduce with
`python benchmarks/neighbor_runtime.py --output /tmp/prime-neighbors.json`.
The [final receipt](benchmarks/results/neighbors_20260913_final.json) records
the immutable baseline revision, source SHA-256 hashes, every numerical
payload, and exact/ANN comparisons. The earlier `neighbors_20260913.json`
receipt predates reuse of the forward leader's querywise contributions and
is retained only as a development measurement.

## Previous cross-language validation: 2026-09-12

- Thirteen unit and integration tests cover proportional allocation, monotone
  projection, BIC-gated ridge initialization, paired querywise forward and
  backward calibration, noise-adjusted subset risks, zero padding, explicit
  query indices, subset budgets, input validation, fitted-object plotting, and
  raw forward ranking below the zero-signal boundary.
- A fixed 320-by-8 correlated dataset produces the same selected coordinates,
  noise floor, raw and projected subset values, and PRIME effects as the R
  implementation. The maximum cross-language numerical difference is
  `7.77e-16`.
- Twenty random three-coordinate set functions, including nonmonotone inputs
  requiring projection, agree with the R implementation with maximum absolute
  allocation difference `6.66e-16`.
- The isolated source distribution and pure-Python wheel build successfully.
- The built wheel installs into a fresh environment and exposes version
  `0.2.0` and the public API.

GitHub Actions runs Python 3.10, 3.12, and 3.13 on Linux and Python 3.12 on
macOS and Windows. It also lints the source and builds both distribution
formats.
