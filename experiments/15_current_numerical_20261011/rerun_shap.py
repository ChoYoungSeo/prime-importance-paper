#!/usr/bin/env python3
r"""Optional fresh full-sample TreeSHAP replay; the default is read-only checking.

Run from any directory (paths are resolved relative to this repository):
  python experiments/15_current_numerical_20261011/rerun_shap.py --check-only
  python experiments/15_current_numerical_20261011/rerun_shap.py \
      --run --case-id friedman_n100_rho0.9_r001 --output-dir build/shap_example

Only --run imports the pinned numerical dependencies or fits models. Each replay
uses 18 fresh inner-CV fits and one fresh full-n fit, then explains those same n
rows. No background, normalization, preprocessing, outer holdout, or model reuse
is introduced. Replay timings are diagnostics, not a reproduction of the source
benchmark's warmups, deterministic repeats, or timing session.

Payload schema: cases.json and results/shap_cases.json are lists indexed by case_id.
Input paths in cases.json are relative to results/numerical_20261011. CSV hashes
refer to decompressed original bytes. Optional saved inner_train_indices and
inner_validation_indices may be present in either list; these are zero-based.
When only validation indices are saved, training indices are their sorted complement.
"""

from __future__ import annotations

import argparse
import csv
import gc
import gzip
import hashlib
import importlib.metadata
import io
import json
import math
import os
from pathlib import Path
import re
import statistics
import sys
import time
from types import SimpleNamespace
import warnings


REPO = Path(__file__).resolve().parents[2]
PAYLOAD = REPO / "results" / "numerical_20261011"
DEFAULT_OUTPUT = Path("build/shap_replay_20261011")
CANDIDATES = [
    {"n_estimators": n, "max_depth": d}
    for n in (200, 500, 1000) for d in (2, 4)
]
FIXED = {
    "objective": "reg:squarederror", "learning_rate": 0.05, "reg_lambda": 1,
    "min_child_weight": 1, "subsample": 1, "colsample_bytree": 1,
    "colsample_bylevel": 1, "colsample_bynode": 1, "tree_method": "hist",
    "n_jobs": 1, "verbosity": 0,
}
VERSIONS = {
    "numpy": "1.26.4", "pandas": "2.2.3", "scipy": "1.14.1",
    "scikit-learn": "1.5.2", "xgboost": "2.1.4", "shap": "0.46.0",
    "lime": "0.2.0.1", "threadpoolctl": "3.5.0",
}
SHAP = {
    "data": None, "feature_perturbation": "tree_path_dependent",
    "model_output": "raw", "approximate": False,
    "aggregation": "np.abs(values).mean(axis=0)", "normalization": "none",
}
PHASES = ("selection_seconds", "final_fit_seconds", "explanation_seconds", "total_seconds")


def require(condition, message):
    """Unlike assert, input checks remain enabled under python -O."""
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def indexed_rows(path):
    rows = read_json(path)
    require(isinstance(rows, list) and rows, f"{path.name}: expected a nonempty list")
    indexed = {}
    for row in rows:
        require(isinstance(row, dict), f"{path.name}: expected case objects")
        case_id = row.get("case_id")
        require(isinstance(case_id, str) and re.fullmatch(r"[A-Za-z0-9_.-]+", case_id),
                f"{path.name}: invalid case_id")
        require(case_id not in indexed, f"{path.name}: duplicate {case_id}")
        indexed[case_id] = row
    return indexed


def contained(path, root):
    return path == root or root in path.parents


def resolve_input(relative):
    require(isinstance(relative, str), "input_path must be a relative string")
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts,
            "input_path must remain inside the bundled payload")
    path = (PAYLOAD / path).resolve()
    require(contained(path, PAYLOAD.resolve()), "input_path escapes the payload")
    require(path.name.endswith(".csv.gz"), "Expected an original CSV compressed as .csv.gz")
    return path


def resolve_output(value):
    path = Path(value).expanduser()
    path = (REPO / path if not path.is_absolute() else path).resolve()
    repo = REPO.resolve()
    require(not contained(repo, path), "Output cannot be the repository or its ancestor")
    require(not contained(path, PAYLOAD.resolve()), "Output cannot be inside the source payload")
    require(not contained(path, (REPO / "manuscript").resolve()),
            "Output cannot be inside the manuscript")
    if contained(path, repo):
        require(contained(path, (REPO / "build").resolve()),
                "Repository-local output must be under build/")
    return path


def expected_seeds(case_id, p):
    seed = (20261008 + int(sha256(case_id.encode("utf-8"))[:7], 16)) % 2**32
    return {
        "case_seed": seed,
        "inner_fold_seed": (seed + 50021) % 2**32,
        "inner_fit_seeds": [(seed + 1000003 * (j + 1) + 7919) % 2**32 for j in range(3)],
        "fit_seed": (seed + ((1 << p) - 1) * 7919) % 2**32,
    }


def select_one_se(rows):
    """Apply the original tie-break and 1e-14 eligibility tolerance."""
    key = lambda row: (row["params"]["n_estimators"], row["params"]["max_depth"])
    best = min(rows, key=lambda row: (-row["mean_inner_r2"], *key(row)))
    cutoff = best["mean_inner_r2"] - best["se_inner_r2"]
    selected = min((row for row in rows if row["mean_inner_r2"] >= cutoff - 1e-14), key=key)
    return selected, best, cutoff


def check_saved_folds(row, n, case_id):
    train = row.get("inner_train_indices")
    valid = row.get("inner_validation_indices")
    if train is None and valid is None:
        return None
    require(isinstance(valid, list) and len(valid) == 3 and
            all(isinstance(va, list) and all(type(i) is int for i in va) for va in valid),
            f"{case_id}: saved validation folds must contain three integer lists")
    if train is None:
        train = [sorted(set(range(n)) - set(va)) for va in valid]
    require(isinstance(train, list) and isinstance(valid, list) and len(train) == len(valid) == 3,
            f"{case_id}: saved folds must contain three train/validation lists")
    universe = set(range(n))
    validation_rows = []
    for j, (tr, va) in enumerate(zip(train, valid)):
        require(isinstance(tr, list) and isinstance(va, list), f"{case_id}: malformed fold {j}")
        require(all(type(i) is int for i in tr + va), f"{case_id}: fold indices must be integers")
        require(len(set(tr)) == len(tr) and len(set(va)) == len(va),
                f"{case_id}: duplicate fold indices")
        require(set(tr).isdisjoint(va) and set(tr) | set(va) == universe,
                f"{case_id}: invalid train/validation partition")
        require(len(va) == n // 3 + (j < n % 3), f"{case_id}: wrong KFold validation size")
        validation_rows.extend(va)
    require(sorted(validation_rows) == list(range(n)), f"{case_id}: validation coverage mismatch")
    return list(zip(train, valid))


def check_reference(reference, case):
    case_id, n, p = case["case_id"], case["n"], case["p"]
    phi = reference.get("phi")
    require(isinstance(phi, list) and len(phi) == p and
            all(isinstance(x, (float, int)) and math.isfinite(x) and x >= 0 for x in phi),
            f"{case_id}: invalid saved mean-absolute SHAP vector")
    for field in ("n", "p", "feature_names", "source_input_hash", *expected_seeds(case_id, p)):
        if field in reference:
            require(reference[field] == case[field], f"{case_id}: reference {field} mismatch")
    if "input_sha256" in reference:
        require(reference["input_sha256"] == case["csv_sha256"], f"{case_id}: reference CSV hash mismatch")
    for field, expected in {
        "n_train": n, "n_evaluation": n, "in_sample": True,
        "feature_perturbation": "tree_path_dependent", "model_output": "raw",
        "metric": "response_units", "normalization": "none", "n_background": 0,
    }.items():
        if field in reference:
            require(reference[field] == expected, f"{case_id}: reference {field} mismatch")
    for field in ("train_indices", "evaluation_indices"):
        if field in reference:
            require(reference[field] == list(range(n)), f"{case_id}: saved {field} mismatch")
    rows = reference.get("candidate_scores")
    require(isinstance(rows, list) and len(rows) == 6, f"{case_id}: expected six candidate scores")
    require([row["params"] for row in rows] == CANDIDATES, f"{case_id}: candidate grid/order mismatch")
    for row in rows:
        scores = row["inner_fold_r2"]
        require(len(scores) == 3 and all(math.isfinite(v) for v in scores),
                f"{case_id}: invalid three-fold R2 scores")
        for field, value in (("mean_inner_r2", statistics.mean(scores)),
                             ("se_inner_r2", statistics.stdev(scores) / math.sqrt(3))):
            require(math.isfinite(row[field]) and math.isclose(row[field], value, rel_tol=1e-12, abs_tol=1e-12),
                    f"{case_id}: saved {field} is inconsistent")
    selected, best, cutoff = select_one_se(rows)
    require(reference["selected_params"] == selected["params"], f"{case_id}: saved one-SE selection mismatch")
    require(math.isclose(reference["one_se_cutoff"], cutoff, rel_tol=1e-12, abs_tol=1e-12),
            f"{case_id}: saved one-SE cutoff mismatch")
    if "best_mean_candidate" in reference:
        require(reference["best_mean_candidate"] == best["params"], f"{case_id}: saved best candidate mismatch")


def check_case(case, reference):
    case_id = case["case_id"]
    n, p = case.get("n"), case.get("p")
    require(type(n) is int and n >= 6 and type(p) is int and 1 <= p <= 63,
            f"{case_id}: invalid n or p")
    names = case.get("feature_names")
    require(isinstance(names, list) and len(names) == p and
            all(isinstance(name, str) and name for name in names) and len(set(names)) == p,
            f"{case_id}: invalid feature_names")
    for field in ("csv_sha256", "source_input_hash"):
        require(isinstance(case.get(field), str) and re.fullmatch(r"[a-f0-9]{64}", case[field]),
                f"{case_id}: invalid {field}")
    for field, expected in expected_seeds(case_id, p).items():
        require(case.get(field) == expected, f"{case_id}: incorrect {field}")
    path = resolve_input(case["input_path"])
    data = gzip.decompress(path.read_bytes())
    require(sha256(data) == case["csv_sha256"], f"{case_id}: decompressed CSV SHA-256 mismatch")
    rows = list(csv.reader(io.StringIO(data.decode("utf-8"), newline="")))
    require(rows and rows[0] == [f"X{i + 1}" for i in range(p)] + ["Y"],
            f"{case_id}: expected canonical X1..Xp,Y CSV columns")
    require(len(rows) == n + 1 and all(len(row) == p + 1 for row in rows[1:]),
            f"{case_id}: CSV row/column count mismatch")
    numeric_rows = [[float(v) for v in row] for row in rows[1:]]
    require(all(math.isfinite(v) for row in numeric_rows for v in row), f"{case_id}: nonfinite CSV value")
    require(len({row[-1] for row in numeric_rows}) > 1, f"{case_id}: constant response")
    check_reference(reference, case)
    folds = [saved for row in (case, reference)
             if (saved := check_saved_folds(row, n, case_id)) is not None]
    require(not folds or all(saved == folds[0] for saved in folds), f"{case_id}: conflicting saved folds")
    return {"case": case, "reference": reference, "csv_bytes": data,
            "saved_folds": folds[0] if folds else None}


def check_protocol(protocol):
    for field, expected in (("candidates", CANDIDATES), ("fixed_parameters", FIXED),
                             ("package_versions", VERSIONS), ("shap", SHAP)):
        require(protocol.get(field) == expected, f"protocol.json: {field} differs from the frozen SHAP protocol")
    if "native_threads" in protocol:
        require(protocol["native_threads"] == 1, "Only one native thread is permitted")


def load_runtime():
    """This is deliberately unreachable in the default/check-only execution."""
    mismatches = []
    for package, expected in VERSIONS.items():
        try:
            actual = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            actual = "not installed"
        if actual != expected:
            mismatches.append(f"{package}: expected {expected}, found {actual}")
    require(not mismatches, "Pinned dependencies required for --run; install requirements-shap.txt:\n" + "\n".join(mismatches))
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
        os.environ[name] = "1"
    sys.dont_write_bytecode = True
    import numpy as np
    import pandas as pd
    import shap
    from sklearn.metrics import r2_score
    from sklearn.model_selection import KFold
    from threadpoolctl import threadpool_limits
    from xgboost import XGBRegressor
    return SimpleNamespace(np=np, pd=pd, shap=shap, r2_score=r2_score,
                           KFold=KFold, threadpool_limits=threadpool_limits,
                           XGBRegressor=XGBRegressor)


def prepare_runtime_context(context, runtime):
    case = context["case"]
    frame = runtime.pd.read_csv(io.BytesIO(context["csv_bytes"]), float_precision="round_trip")
    X, y = frame.drop(columns="Y").to_numpy(float), frame.Y.to_numpy(float)
    require(X.shape == (case["n"], case["p"]) and len(y) == case["n"], "Parsed dimension mismatch")
    require(runtime.np.isfinite(X).all() and runtime.np.isfinite(y).all() and runtime.np.var(y) > 0,
            "Parsed input must be finite with nonconstant response")
    folds = list(runtime.KFold(3, shuffle=True, random_state=case["inner_fold_seed"]).split(runtime.np.zeros(len(y))))
    if context["saved_folds"] is not None:
        actual = [(tr.tolist(), va.tolist()) for tr, va in folds]
        require(actual == context["saved_folds"], f"{case['case_id']}: seeded KFold differs from saved folds")
    context.update(X=X, y=y, folds=folds)


def pipeline(context, runtime):
    np, case = runtime.np, context["case"]
    start = time.perf_counter()
    candidates = []
    for params in CANDIDATES:
        scores = []
        for j, (train, valid) in enumerate(context["folds"]):
            model = runtime.XGBRegressor(**FIXED, **params, random_state=case["inner_fit_seeds"][j])
            model.fit(context["X"][train], context["y"][train])
            pred = np.asarray(model.predict(context["X"][valid]), float)
            scores.append(float(runtime.r2_score(context["y"][valid], pred)))
        a = np.asarray(scores)
        candidates.append({"params": params, "inner_fold_r2": scores,
                           "mean_inner_r2": float(a.mean()), "se_inner_r2": float(a.std(ddof=1) / np.sqrt(3))})
    selected, best, cutoff = select_one_se(candidates)
    selection_end = time.perf_counter()
    model = runtime.XGBRegressor(**FIXED, **selected["params"], random_state=case["fit_seed"])
    model.fit(context["X"], context["y"])
    fit_end = time.perf_counter()
    explainer = runtime.shap.TreeExplainer(model, data=None, feature_perturbation="tree_path_dependent", model_output="raw")
    values = np.asarray(explainer.shap_values(context["X"], approximate=False, check_additivity=False), float)
    base = float(explainer.expected_value)
    phi = np.abs(values).mean(axis=0)
    end = time.perf_counter()
    return {"model": model, "values": values, "base": base, "phi": phi,
            "candidate_scores": candidates, "selected_params": selected["params"],
            "best_mean_candidate": best["params"], "one_se_cutoff": cutoff,
            "selection_seconds": selection_end - start, "final_fit_seconds": fit_end - selection_end,
            "explanation_seconds": end - fit_end, "total_seconds": end - start}


def save_json(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def save_replay(context, result, output, runtime, caught):
    np, case, reference = runtime.np, context["case"], context["reference"]
    case_id, n, p = case["case_id"], case["n"], case["p"]
    pred = np.asarray(result["model"].predict(context["X"]), float)
    values, phi, base = result["values"], result["phi"], result["base"]
    require(pred.shape == (n,) and values.shape == (n, p) and phi.shape == (p,), "Output shape mismatch")
    require(np.isfinite(pred).all() and np.isfinite(values).all() and np.isfinite(phi).all() and np.isfinite(base),
            "Nonfinite prediction or SHAP output")
    require(np.array_equal(phi, np.abs(values).mean(axis=0)) and np.all(phi >= 0), "Mean-absolute aggregation mismatch")
    gap = float(np.max(np.abs(base + values.sum(axis=1) - pred)))
    require(gap < 1e-4 * max(1., float(np.max(np.abs(pred)))), "SHAP additivity check failed")
    model_path = output / "models" / f"{case_id}.ubj"
    result["model"].save_model(model_path)
    model_hash = sha256(model_path.read_bytes())
    loaded = runtime.XGBRegressor(**FIXED, **result["selected_params"], random_state=case["fit_seed"])
    loaded.load_model(model_path)
    loaded_explainer = runtime.shap.TreeExplainer(loaded, data=None, feature_perturbation="tree_path_dependent", model_output="raw")
    loaded_values = np.asarray(loaded_explainer.shap_values(context["X"], approximate=False, check_additivity=True), float)
    require(np.array_equal(pred, np.asarray(loaded.predict(context["X"]), float)) and
            np.array_equal(values, loaded_values) and base == float(loaded_explainer.expected_value),
            "Serialized model prediction/SHAP reload parity failed")
    require(sha256(model_path.read_bytes()) == model_hash, "Serialized model changed during reload")
    indices = np.arange(n, dtype=np.int64)
    values_path = output / "values" / f"{case_id}.npz"
    np.savez(values_path, values=values, base=base, predictions=pred, train_indices=indices,
             evaluation_indices=indices, background_indices=np.array([], dtype=np.int64))
    ref_phi = np.asarray(reference["phi"], float)
    scores = np.asarray([row["inner_fold_r2"] for row in result["candidate_scores"]])
    require(scores.shape == (6, 3) and np.isfinite(scores).all(), "Nonfinite inner-CV scores")
    ref_scores = np.asarray([row["inner_fold_r2"] for row in reference["candidate_scores"]])
    comparison = {
        "selected_params_match": result["selected_params"] == reference["selected_params"],
        "max_absolute_phi_difference": float(np.max(np.abs(phi - ref_phi))),
        "max_absolute_inner_r2_difference": float(np.max(np.abs(scores - ref_scores))),
        "phi_close": bool(np.allclose(phi, ref_phi, rtol=1e-8, atol=1e-10)),
        "phi_rtol": 1e-8, "phi_atol": 1e-10,
        "note": "Cross-platform numerical differences are reported, not hidden or treated as bitwise reproducibility.",
    }
    row = {
        "status": "COMPLETE", "case_id": case_id, "n": n, "p": p, "feature_names": case["feature_names"],
        "phi": phi.tolist(), "signed_mean": values.mean(axis=0).tolist(), "expected_value": base,
        "metric": "response_units", "normalization": "none", "in_sample": True,
        "training_scope": "full_sample", "evaluation_scope": "same_full_sample_in_sample",
        "n_train": n, "n_evaluation": n, "train_indices": indices.tolist(), "evaluation_indices": indices.tolist(),
        "n_background": 0, "background_indices": [], "feature_perturbation": "tree_path_dependent",
        "model_output": "raw", "approximate": False, "model_reused": False,
        "input_sha256": case["csv_sha256"], "source_input_hash": case["source_input_hash"],
        **expected_seeds(case_id, p),
        "inner_train_indices": [tr.tolist() for tr, _ in context["folds"]],
        "inner_validation_indices": [va.tolist() for _, va in context["folds"]],
        "candidate_scores": result["candidate_scores"], "selected_params": result["selected_params"],
        "best_mean_candidate": result["best_mean_candidate"], "one_se_cutoff": result["one_se_cutoff"],
        "model_path": model_path.relative_to(output).as_posix(), "model_sha256": model_hash,
        "values_path": values_path.relative_to(output).as_posix(), "values_sha256": sha256(values_path.read_bytes()),
        "max_additivity_error": gap, "reloaded_predictions_and_shap_bitwise_equal": True,
        "reference_comparison": comparison, "warnings": [str(w.message) for w in caught],
        **{phase: result[phase] for phase in PHASES},
    }
    save_json(output / "cases" / f"{case_id}.json", row)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-only", action="store_true", help="Read-only validation (the default); no numerical packages needed")
    mode.add_argument("--run", action="store_true", help="Explicitly enable fresh model fitting with pinned dependencies")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--case-id", action="append", help="Case to check/replay; repeat for multiple cases")
    selection.add_argument("--all", action="store_true", help="Select every bundled SHAP case in saved order")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT), help="New output directory; relative paths use the repository root")
    args = parser.parse_args(argv)
    if args.run and not (args.case_id or args.all):
        parser.error("--run requires --case-id or --all; omit --run for safe checking")
    output = resolve_output(args.output_dir)
    protocol = read_json(PAYLOAD / "protocol.json")
    check_protocol(protocol)
    cases = indexed_rows(PAYLOAD / "cases.json")
    references = indexed_rows(PAYLOAD / "results" / "shap_cases.json")
    ids = args.case_id or list(references)
    require(len(ids) == len(set(ids)), "Duplicate --case-id arguments are not allowed")
    require(all(case_id in cases and case_id in references for case_id in ids), "Unknown or non-SHAP case_id")
    if "case_ids" in protocol:
        require(list(references) == protocol["case_ids"], "Saved SHAP case list/order differs from protocol")
    contexts = [check_case(cases[case_id], references[case_id]) for case_id in ids]
    if not args.run:
        print(json.dumps({"status": "CHECK_ONLY_PASS", "case_count": len(contexts), "model_fits": 0,
                          "dependency_checks": "not required in check-only mode",
                          "checks": ["original CSV hashes", "finite input and dimensions", "row counts",
                                     "seed formulas", "saved one-SE selection", "optional saved fold partitions"],
                          "saved_folds_cases": sum(ctx["saved_folds"] is not None for ctx in contexts),
                          "note": "Exact seeded sklearn KFold reconstruction is additionally checked before --run fitting."}, indent=2))
        return 0
    runtime = load_runtime()
    for context in contexts:
        prepare_runtime_context(context, runtime)
    require(not output.exists(), "Output directory already exists; choose a new directory to avoid overwrites")
    output.mkdir(parents=True, exist_ok=False)
    for folder in ("cases", "models", "values"):
        (output / folder).mkdir()
    save_json(output / "replay_protocol.json", {
        "case_ids": ids, "candidates": CANDIDATES, "fixed_parameters": FIXED, "package_versions": VERSIONS,
        "shap": SHAP, "native_threads": 1, "processes": 1, "python": sys.version,
        "bundled_protocol_sha256": sha256((PAYLOAD / "protocol.json").read_bytes()),
        "bundled_cases_sha256": sha256((PAYLOAD / "cases.json").read_bytes()),
        "bundled_shap_cases_sha256": sha256((PAYLOAD / "results" / "shap_cases.json").read_bytes()),
        "script_sha256": sha256(Path(__file__).read_bytes()),
        "scope": "Fresh three-fold CV, one-SE selection, full-n fit, and same-full-n mean absolute TreeSHAP.",
        "timing_scope": "Diagnostic single passes; no benchmark warmups or independent deterministic repeat fits.",
    })
    summaries = []
    with runtime.threadpool_limits(limits=1):
        for context in contexts:
            gc.collect()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = pipeline(context, runtime)
            row = save_replay(context, result, output, runtime, caught)
            summaries.append({"case_id": row["case_id"], "reference_comparison": row["reference_comparison"],
                              **{phase: row[phase] for phase in PHASES}})
            print(json.dumps({"status": "CASE_COMPLETE", "completed": len(summaries), "total": len(ids),
                              "case_id": row["case_id"], "reference_comparison": row["reference_comparison"]}), flush=True)
            del result
    summary = {"status": "COMPLETE", "case_count": len(summaries), "model_fits": 19 * len(summaries),
               "all_selected_params_match": all(row["reference_comparison"]["selected_params_match"] for row in summaries),
               "all_phi_close": all(row["reference_comparison"]["phi_close"] for row in summaries), "cases": summaries}
    save_json(output / "summary.json", summary)
    print(json.dumps({key: value for key, value in summary.items() if key != "cases"}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
