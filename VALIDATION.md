# Validation

`make validate-release` checks the public scientific inputs, saved estimates
and current four-panel comparisons without fitting models. It runs two
validators:

- `results/numerical_20261001`: scientific input and retained-file hashes,
  complete allocation vectors, paired inputs and standardization, target
  normalization and ratios, and the earlier PRIME, PME-NA and FIRST-SHAP
  results. The nine-condition, three-method Table 2 numeric comparison is checked
  against the existing full-precision
  `friedman/earlier_four_method_table.csv` (`PRIME`, `Full LL` as PRIME-Full,
  and `PME-kNN` as PME-NA), covering all 27 displayed values.
- `results/numerical_20261011`: paired input files, saved SHAP and
  cross-validation results, timing records, four-panel figure data and box
  statistics, and retained approved PNG/PDF research graphic assets.

Paper prose, full paper PDFs, TeX/Bib sources and typeset table source are
`NOT_IN_PUBLIC_EXPORT`. Their checks are outside the public validation scope
and are not marked as passed. The numeric comparison remains fully checked;
removing paper files does not remove its scientific input or saved-result
checks.

These saved-result checks require Python 3.9 or newer. `make artifacts`
rebuilds research comparison PDF and PNG graphics only in
`build/current_numerical_20261011/` from saved data, using R with ggplot2,
patchwork and jsonlite. It does not refit models or compile a paper.

The recorded R calculations use the frozen PRIME 0.5.2.9101 source.
From the repository root, run:

```sh
make deps
make bootstrap
make verify-workflow
make validate-release
```

Bootstrap compiles the pinned package into a local library, runs a bounded
smoke check, and records installed files and dependency versions in
`.verified/local_install.json`. Workflow checks verify the source and input
manifests and the installation receipt. Historical manifests are preserved;
reviewed documentation and orchestration changes have separate workflow
manifests.

The CI workflow includes package checks, workflow tests, bootstrap and public
saved-result validators. The separately documented frozen-input replay
interfaces in the [numerical evidence guide](docs/numerical-evidence.md#reproduction-and-provenance)
are not invoked by `make validate-release` or `make artifacts`. Timing records
describe their original hardware and runtime. No new scientific experiments or
complete model refits were run during preparation of this export.
