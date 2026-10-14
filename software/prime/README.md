# PRIME: A New Measure for Variable Importance

This directory contains PRIME 0.5.2.9101, the R implementation used for the
paper's experiments. Its scientific source is pinned by the release manifests.
The maintained [R package](https://github.com/ChoYoungSeo/prime-importance-r)
and [Python package](https://github.com/ChoYoungSeo/prime-importance) provide
PRIME for new analyses.

Run these commands from the **paper repository root**:

```sh
make deps
make bootstrap
make verify-workflow
make validate-release
```

Bootstrap builds the pinned package in the repository's local R library and
checks its installed files. The [reproduction guide](../../docs/research-materials.md)
describes the saved-result checks and estimator replay.

- [R reference manual](../../docs/prime.pdf)
- [Quickstart source](inst/examples/quickstart.R)
- [Method specification](../../METHOD_SPEC.md)
- [Citation](../../CITATION.cff)
- [MIT license](../../LICENSE)
