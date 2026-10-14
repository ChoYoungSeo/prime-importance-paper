# Concrete observed-data release evidence

The active main-paper figure uses 100 paired samples of 824 distinct rows
drawn without replacement from the 1,030-row Concrete dataset. Water is
replaced by water/(cement + slag + fly ash), leaving eight coordinates. The
five methods are PRIME, Full LL-PME, PME-kNN, Shapley-kNN on the same
projected 2NN game, and FIRST-SHAP. The Concrete data have no known
population allocation, so paired half-L1 values describe agreement between
estimated vectors rather than error against truth.

`distinct_824/sample_inputs.csv` contains the ordered row indices and input
hash for each replication. `paired_phi.csv` contains all 4,000 coordinate
allocations, `runs.csv` records all 500 method fits, and the coordinate,
pair, and timing summaries are computed from them. `VALIDATION.json` records
the original exact seed and input-hash replays. The release validator checks
index distinctness, all 500 statuses, every vector sum, all 40 coordinate
means/medians/5th and 95th percentiles/selection frequencies, and all 1,000
paired half-L1 distances and their summaries. The main figure is generated
from `distinct_824/coordinate_summary.csv` by
`experiments/13_concrete_submission/build_concrete_p8_figure.py`;
`FIGURE_RECEIPT.json` pins its CSV and PDF hashes.

`bootstrap_1030/` preserves the earlier 200 paired, size-1,030 bootstrap
summary and full-data fit as a sensitivity comparison. Both sample size and
replacement differ from the main protocol. `distinct_824/neighbor_full_findings.json`
records the duplicate-neighbor audit: the full-coordinate 2NN zero-local-variance
fraction is 66.81% in the bootstrap and 13.54% in the distinct-row samples.
The main figure uses only the latter. The separately archived
`distinct_824/FILE_HASHES.csv` names upstream checkpoints and logs that are
not included in this compact release; it is retained as the original hash
receipt, and `PORTABLE_SOURCES.json` maps the relevant archived source
hashes to paths in this checkout.

The Concrete data and PME helper in this repository exactly match the source
hashes in the study. The local FIRST reporting adapter is included at
`experiments/13_concrete_submission/first_adapter.R`. FIRST-SHAP used the
official supplementary `code/lib/first.R` of Huang and Joseph,
[*Factor Importance Ranking and Selection Using Total Indices*](https://doi.org/10.1080/00401706.2025.2483531)
(2025), SHA-256
`759c110ef16e4cce1fcebf17e1d1ed03035a90c92878e02715ff038f1afe1ee4`.
That third-party file is not redistributed here; obtain the article's
supplementary `code.zip` and verify the hash before a full method rerun.
The copied source receipt retains the original absolute archive paths for
historical identification; `PORTABLE_SOURCES.json` supplies the checkout
paths used by the release validator.

Run `make validate-release` at the repository root to check the packaged
evidence and compile both manuscript cases. To redraw the figure, install
ReportLab 4.4.9, then run `make figure-concrete`. The builder uses Arial when
available, matching the pinned release PDF hash; its Liberation Sans and
DejaVu Sans fallbacks retain the exact plotted values with different PDF
bytes.
