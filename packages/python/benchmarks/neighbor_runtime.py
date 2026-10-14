"""Paired, single-process-at-a-time old/exact/ANN timing and peak-RSS audit.

Run from the repository with a Python environment containing NumPy and SciPy:
    python benchmarks/neighbor_runtime.py --output /tmp/prime-neighbors.json

The baseline is read from the named git revision into a temporary source tree.
Wall time excludes imports and input generation; peak RSS covers the fresh
worker's lifetime, including imports/input generation. No time limit is a
subsampling instruction: unsuccessful workers fail the run explicitly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def worker(args):
    import resource
    import time
    from itertools import combinations

    sys.path.insert(0, str(Path(args.source_root) / "src"))
    import numpy as np

    from prime import prime
    from prime._core import local_risk, neighborhood_size, prepare_data

    sizes = {"single_subset": (15000, 6), "subset_cache": (3000, 5), "fit": (1500, 10)}
    n, p = sizes[args.case]
    rng = np.random.default_rng(20260913)
    X = rng.normal(size=(n, p))
    y = X[:, 0] + 0.7 * X[:, 1] ** 2 + np.sin(X[:, 2]) + rng.normal(0, 0.1, n)
    options = {} if args.mode == "old" else {
        "neighbor_eps": 0.2 if args.mode == "ann" else 0.0,
        "query_block_size": 256,
    }
    start = time.perf_counter()
    if args.case == "fit":
        result = prime(X, y, **options)
        payload = {"selected": result.selected, "effects": dict(result.effects)}
    else:
        data = prepare_data(X, y, feature_names=None, rescale=True,
                            query_indices=None, lambda_=1e-6, **options)
        subsets = [tuple(range(p))] if args.case == "single_subset" else [
            subset for d in range(1, p + 1) for subset in combinations(range(p), d)
        ]
        risks = [local_risk(data, subset, neighborhood_size(n, len(subset)))
                 for subset in subsets]
        scalar_entries = sum(isinstance(value, float) for value in data.cache.values())
        retained = sum(getattr(value, "contribution", np.empty(0)).nbytes
                       for value in data.cache.values())
        retained += sum(value.contribution.nbytes
                        for value in getattr(data, "detail_cache", {}).values())
        payload = {"risks": risks, "scalar_cache_entries": scalar_entries,
                   "retained_contribution_bytes": retained}
    elapsed = time.perf_counter() - start
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        peak_rss *= 1024
    print(json.dumps({"case": args.case, "mode": args.mode, "n": n, "p": p,
                      "wall_seconds": elapsed, "peak_rss_bytes": peak_rss,
                      "numpy": np.__version__, "payload": payload}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--baseline", default="HEAD")
    parser.add_argument("--total-seconds", type=float, default=60)
    parser.add_argument("--source-root")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--case", choices=["single_subset", "subset_cache", "fit"])
    parser.add_argument("--mode", choices=["old", "exact", "ann"])
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.output is None:
        parser.error("--output is required")
    root = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(
        ["git", "rev-parse", args.baseline], cwd=root, text=True
    ).strip()
    environment = os.environ.copy()
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS"]:
        environment[key] = "1"
    rows = []
    deadline = time.perf_counter() + args.total_seconds
    with tempfile.TemporaryDirectory(prefix="prime-neighbor-baseline-") as temporary:
        baseline = Path(temporary)
        paths = subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", revision, "src"], cwd=root, text=True
        ).splitlines()
        for relative in paths:
            target = baseline / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(subprocess.check_output(
                ["git", "show", f"{revision}:{relative}"], cwd=root
            ))
        for case in ["single_subset", "subset_cache", "fit"]:
            for mode in ["old", "exact", "ann"]:
                source = baseline if mode == "old" else root
                completed = subprocess.run(
                    [sys.executable, str(Path(__file__).resolve()), "--worker",
                     "--source-root", str(source), "--case", case, "--mode", mode],
                    env=environment, capture_output=True, text=True, check=True,
                    timeout=max(0.01, min(120, deadline - time.perf_counter())),
                )
                row = json.loads(completed.stdout)
                rows.append(row)
                print(case, mode, f'{row["wall_seconds"]:.3f}s',
                      f'{row["peak_rss_bytes"] / 2**20:.1f} MiB', flush=True)
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(
                    {"baseline_revision": revision, "seed": 20260913,
                     "threads": 1, "replications": 1, "python": sys.version,
                     "platform": sys.platform,
                     "source_sha256": {
                         str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in sorted((root / "src").rglob("*.py"))
                     }, "rows": rows}, indent=2
                ) + "\n")
    comparisons = []
    for case in ["single_subset", "subset_cache", "fit"]:
        modes = {row["mode"]: row for row in rows if row["case"] == case}
        key = "effects" if case == "fit" else "risks"
        original = modes["old"]["payload"][key]
        exact = modes["exact"]["payload"][key]
        approximate = modes["ann"]["payload"][key]
        assert original == exact, f"Exact numerical mismatch for {case}"
        if isinstance(original, dict):
            differences = [abs(original[name] - approximate[name]) for name in original]
        else:
            differences = [abs(x - y) for x, y in zip(original, approximate, strict=True)]
        comparison = {"case": case, "exact_equal": True,
                      "ann_max_absolute_difference": max(differences)}
        if case == "fit":
            comparison["ann_selected_equal"] = (
                modes["old"]["payload"]["selected"] == modes["ann"]["payload"]["selected"]
            )
        comparisons.append(comparison)
    receipt = json.loads(args.output.read_text())
    receipt["comparisons"] = comparisons
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
