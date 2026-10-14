# Current validation scope

The public release pins PRIME 0.5.2.9101, all package build sources, and corrected
numerical inputs through `results/current_release.json`. Original scientific
input hashes and historical validation reports remain in the lineage records.

The paper engine includes the reviewed row-identity, singleton, feature-name,
and empty-selection contract repairs. The separate optimized R development
engine can choose different neighbors at exact distance ties; universal
cross-engine parity is not claimed. This repository pins the paper engine.

For current checks, run `make bootstrap`, `make verify-workflow`, and
`make validate-release`. Bootstrap compiles pinned source and checks the loaded
installation with a bounded smoke. Workflow tests reject altered source,
transport, decoded data, receipts, and installed files. Strict release validation
recalculates coordinate-level losses and verifies 12 numerical prose blocks and
56 table rows. CI also runs `R CMD check --no-manual`.

The public packaging changes preserve all 63 scientific input hashes and add one
fixed Ackley reference for the experiment validators. Full scientific grids are
not rerun for these transport and installation changes. Historical package and
simulation claims are retained under `release_lineage/validation_report_20260902.md`.
