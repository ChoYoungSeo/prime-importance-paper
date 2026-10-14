# Published numerical inputs

This directory contains the corrected numerical inputs selected by
`results/current_release.json`, with immutable historical verification baselines
kept separately. `REVISION_LINEAGE.json` describes the scientific repair lineage.

The original 63 logical input SHA-256 values are retained. The fixed Ackley
reference required by the experiment validator is added as input 64. The large
`scaling_phi.csv` is distributed as deterministic `.csv.gz`; the verified workflow
checks both transport integrity and the original uncompressed hash, then hydrates
the original CSV name inside an isolated work directory.

From the repository root, run `make deps`, `make bootstrap`, and
`make validate-release`. See the root README for build prerequisites and the
separate, optional scientific rerun commands. Historical installed-binary hashes
remain lineage evidence; each fresh source build receives its own local receipt.
