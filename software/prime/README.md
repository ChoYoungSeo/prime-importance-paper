# PRIME: A New Measure for Variable Importance

This directory contains PRIME 0.5.2.9101, the R implementation used for the
paper's experiments. Its scientific source is pinned by the release manifests.
For new analyses, use the [maintained R and Python packages](../../README.md#install).

Run these commands from the **paper repository root**:

```sh
make deps
make bootstrap
make verify-workflow
make validate-release
```

Bootstrap builds the pinned package in the repository's local R library and
checks its installed files. The [numerical evidence guide](../../docs/numerical-evidence.md#reproduction-and-provenance)
describes the saved-result checks and archived replay source.

- [R reference manual](../../docs/prime.pdf)
- [Quickstart source](inst/examples/quickstart.R)
- [Method specification](../../METHOD_SPEC.md)
- [Citation](../../CITATION.cff)
- [MIT license](../../LICENSE)
