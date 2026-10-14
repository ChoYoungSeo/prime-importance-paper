#!/usr/bin/env python3
"""Validate the September 2026 manuscript and its reviewed evidence."""

import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT = ROOT / "manuscript"
EVIDENCE = ROOT / "results/submission_20260920"
CONCRETE = EVIDENCE / "concrete"
CONCRETE_DISTINCT = CONCRETE / "distinct_824"
CONCRETE_BOOTSTRAP = CONCRETE / "bootstrap_1030"
TIMING_PARITY = EVIDENCE / "friedman/timing_parity"
METHODS = {"PRIME", "Full LL-PME", "PME-kNN", "Shapley-kNN", "FIRST-SHAP"}
VARIABLES = (
    "cement", "slag", "ash", "water_binder_ratio", "superplastic",
    "coarseagg", "fineagg", "age",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_sha256(path):
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def finite(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def quantile(values, probability):
    """R's default type-7 empirical quantile."""
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[min(lower + 1, len(ordered) - 1)] * fraction


def close(actual, expected, message, tolerance=1e-11):
    require(abs(float(actual) - float(expected)) <= tolerance, message)


def selected_case():
    selector = (MANUSCRIPT / "sections/observed_case_config.tex").read_text()
    settings = re.findall(r"(?m)^\\ConcreteMain(true|false)\s*$", selector)
    require(len(settings) == 1, "Set exactly one observed-data case")
    for path in (MANUSCRIPT / "main.tex", MANUSCRIPT / "supplement.tex"):
        require(r"\input{sections/observed_case_config}" in path.read_text(),
                f"Missing common case switch in {path.name}")
    for relative, concrete, auto in (
        ("sections/06_experiments.tex", "new_concrete", "new_auto_mpg"),
        ("sections/supplement_experiments.tex", "supplement_concrete_case", "supplement_auto_mpg_case"),
    ):
        source = (MANUSCRIPT / relative).read_text()
        require(r"\ifConcreteMain" in source and r"\else" in source and r"\fi" in source,
                f"Missing case branch in {relative}")
        require(r"\input{sections/" + concrete + "}" in source and
                r"\input{sections/" + auto + "}" in source,
                f"Incomplete case branch in {relative}")
    return "concrete" if settings[0] == "true" else "auto_mpg"


def validate_legacy_inputs():
    release = json.loads((ROOT / "results/current_release.json").read_text())
    require(len(release["files"]) == 64, "Legacy release must retain all 64 inputs")
    data_root = ROOT / release["data_root"]
    checked = 0
    for logical, record in release["files"].items():
        if "transport" in record:
            transport = record["transport"]
            source = data_root / transport["path"]
            require(source.is_file(), f"Missing legacy transport: {logical}")
            require(sha256(source) == transport["sha256"],
                    f"Legacy transport hash mismatch: {logical}")
            require(gzip_sha256(source) == record["sha256"],
                    f"Legacy uncompressed hash mismatch: {logical}")
        else:
            source = data_root / record["path"]
            require(source.is_file(), f"Missing legacy input: {logical}")
            require(sha256(source) == record["sha256"],
                    f"Legacy input hash mismatch: {logical}")
        checked += 1
    return checked


def validate_matched_k():
    full = rows(EVIDENCE / "matched_k/full_summary.csv")
    require(len(full) == 270, "Matched-K summary must contain 270 datasets")
    keys = [(r["task_id"], r["n"], r["rho"], r["replication"]) for r in full]
    require(len(set(keys)) == 270, "Matched-K task keys are not unique")
    counts = Counter((int(r["n"]), float(r["rho"])) for r in full)
    require(set(counts.values()) == {30} and len(counts) == 9,
            "Matched-K cells must contain 30 datasets")
    metrics = [
        "pme_half_l1", "shapley_half_l1",
        "pme_proxy_mass", "shapley_proxy_mass",
    ]
    require(all(finite(r[m]) for r in full for m in metrics),
            "Matched-K metrics must be finite")

    cells = rows(EVIDENCE / "matched_k/cell_means.csv")
    require(len(cells) == 9, "Matched-K cell summary must contain nine cells")
    for cell in cells:
        subset = [
            r for r in full
            if int(r["n"]) == int(cell["n"])
            and abs(float(r["rho"]) - float(cell["rho"])) < 1e-12
        ]
        require(len(subset) == 30, "Matched-K cell membership mismatch")
        for metric in metrics:
            mean = sum(float(r[metric]) for r in subset) / len(subset)
            require(abs(mean - float(cell[metric])) < 1e-12,
                    f"Matched-K aggregate mismatch: {cell['n']} {cell['rho']} {metric}")

    combined = rows(EVIDENCE / "matched_k/pme_accuracy_four_methods.csv")
    require(len(combined) == 36, "Four-method PME table source must have 36 rows")
    by_cell = {}
    for row in combined:
        by_cell.setdefault((int(row["n"]), float(row["rho"])), {})[
            row["method"]
        ] = float(row["mean"])
    require(len(by_cell) == 9, "Four-method PME table must have nine cells")
    for key, values in by_cell.items():
        require(set(values) == {
            "PRIME", "Full LL-PME", "Full LC-K_n", "PME-kNN"
        }, f"Missing method in PME cell {key}")
        require(values["Full LL-PME"] < values["Full LC-K_n"],
                f"Matched-K ordering failed in cell {key}")
        require(values["PRIME"] == min(values.values()),
                f"PRIME is not the strict cell winner in {key}")

    receipt = json.loads((EVIDENCE / "matched_k/FULL_RECEIPT.json").read_text())
    require(receipt["status"] == "PASS" and receipt["datasets"] == 270,
            "Matched-K receipt failed")
    timing = json.loads(
        (EVIDENCE / "matched_k/PURE_PME_TIMING_RECEIPT.json").read_text()
    )
    require(timing["status"] == "PASS" and timing["complete"] == 30,
            "Matched-K timing receipt failed")
    return {
        "datasets": 270,
        "cells": 9,
        "full_ll_better_cells": 9,
        "central_pme_error": next(
            float(r["pme_half_l1"])
            for r in cells if int(r["n"]) == 200 and float(r["rho"]) == 0.5
        ),
        "archived_pure_pme_median_time": float(timing["median_sec"]),
    }


def validate_auto_mpg():
    full = rows(EVIDENCE / "auto_mpg/full_data_allocations.csv")
    bars = rows(EVIDENCE / "auto_mpg/figure8_fullbars_stability_values.csv")
    summary = rows(EVIDENCE / "auto_mpg/method_variable_summary.csv")
    require(len(full) == len(bars) == len(summary) == 32,
            "Auto MPG summaries must contain 32 method-coordinate rows")
    full_map = {
        (r["method"], r["variable_name"]): float(r["estimate"]) for r in full
    }
    for row in bars:
        key = (row["method"], row["variable_name"])
        require(key in full_map, f"Missing Auto MPG full-data value: {key}")
        require(abs(float(row["full_data_estimate"]) - full_map[key]) < 1e-12,
                f"Auto MPG bar mismatch: {key}")
        require(float(row["empirical_q05"]) <= float(row["empirical_q95"]),
                f"Auto MPG interval reversal: {key}")
        require(int(row["repeated_count"]) == 100,
                f"Auto MPG repeated count mismatch: {key}")
    selection = {
        (r["method"], r["variable_name"]): float(r["selected_frequency"])
        for r in summary
    }
    for variable, expected in {
        "displacement": 1.0,
        "model_year": 1.0,
        "horsepower": 0.93,
        "weight": 0.89,
    }.items():
        require(abs(selection[("PRIME", variable)] - expected) < 1e-12,
                f"Auto MPG selection frequency mismatch: {variable}")
    receipt = json.loads(
        (EVIDENCE / "auto_mpg/FULLBARS_BUILD_RECEIPT.json").read_text()
    )
    require(receipt["pass_build_checks"] is True,
            "Auto MPG figure receipt failed")
    return {"records": 404, "full_data_rows": 32, "subsamples": 100}


def validate_concrete():
    provenance = json.loads((CONCRETE / "PORTABLE_SOURCES.json").read_text())
    for relative, expected in provenance["sha256"].items():
        path = ROOT / relative
        require(path.is_file() and sha256(path) == expected,
                f"Concrete source or summary hash mismatch: {relative}")

    archived = json.loads((CONCRETE_DISTINCT / "SOURCE_MANIFEST.json").read_text())
    for name, relative in provenance["archived_source_paths"].items():
        require(archived["sources"][name]["sha256"] == provenance["sha256"][relative],
                f"Concrete archive source differs from portable copy: {name}")
    require(archived["sources"]["official_first"]["sha256"] ==
            provenance["external_sources"]["official_first"]["sha256"],
            "Pinned official FIRST supplement hash changed")
    require(sha256(CONCRETE_DISTINCT / "FILE_HASHES.csv") ==
            archived["file_hashes_csv_sha256"], "Concrete archive hash receipt changed")
    figure = json.loads((CONCRETE / "FIGURE_RECEIPT.json").read_text())
    require(figure["source_csv"] ==
            "results/submission_20260920/concrete/distinct_824/coordinate_summary.csv" and
            figure["source_sha256"] == sha256(CONCRETE_DISTINCT / "coordinate_summary.csv")
            and figure["output_pdf"] == "manuscript/figures/concrete_p8_distinct_subsamples.pdf"
            and figure["output_sha256"] == sha256(ROOT / figure["output_pdf"])
            and figure["median_markers"] == figure["percentile_ranges"] == 40,
            "Concrete figure receipt or rendered figure changed")

    study = json.loads((CONCRETE_DISTINCT / "VALIDATION.json").read_text())
    require(study["complete_calls"] == study["expected_calls"] == 500 and
            study["error_calls"] == study["unavailable_calls"] == 0,
            "Concrete five-method completion receipt failed")
    require(study["source_rows"] == 1030 and study["rows_per_replication"] == 824
            and study["p"] == 8 and study["all_sample_indices_distinct"]
            and study["exact_seed_index_replays"] == 100
            and study["exact_input_hash_replays"] == 100
            and study["distinct_input_hashes"] == 100,
            "Concrete input replay receipt failed")
    for key, relative in (
        ("concrete_sha256", "data/real_data_cache/concrete.csv"),
        ("pme_helper_sha256", "experiments/03_pme/00_common.R"),
        ("first_adapter_sha256", "experiments/13_concrete_submission/first_adapter.R"),
    ):
        require(study[key] == sha256(ROOT / relative), f"Concrete pinned source changed: {key}")
    require(study["official_first_sha256"] ==
            provenance["external_sources"]["official_first"]["sha256"],
            "Concrete FIRST source receipt changed")

    inputs = rows(CONCRETE_DISTINCT / "sample_inputs.csv")
    require(len(inputs) == 100, "Concrete must have 100 sample-input records")
    input_hashes = set()
    for record in inputs:
        replication = int(record["replication"])
        indices = [int(value) for value in record["indices"].split()]
        require(1 <= replication <= 100 and int(record["seed"]) == 202609220 + replication,
                f"Concrete seed mismatch: {replication}")
        require(int(record["n_source"]) == 1030 and int(record["n"]) == 824 and
                int(record["p"]) == 8 and len(indices) == len(set(indices)) == 824
                and min(indices) >= 1 and max(indices) <= 1030,
                f"Concrete sample indices invalid: {replication}")
        require(re.fullmatch(r"[a-f0-9]{64}", record["input_hash"]) is not None,
                f"Concrete input hash invalid: {replication}")
        input_hashes.add(record["input_hash"])
    require({int(r["replication"]) for r in inputs} == set(range(1, 101)) and
            len(input_hashes) == 100, "Concrete input records are not unique")

    runs = rows(CONCRETE_DISTINCT / "runs.csv")
    require(len(runs) == 500, "Concrete must have 500 method records")
    run_keys = {(int(r["replication"]), r["method"]) for r in runs}
    require(run_keys == {(r, m) for r in range(1, 101) for m in METHODS},
            "Concrete method records are incomplete")
    require(all(r["status"] == "complete" and
                int(r["seed"]) == 202609220 + int(r["replication"]) and
                finite(r["elapsed_sec"]) for r in runs),
            "Concrete method status, seed, or timing invalid")

    phi = rows(CONCRETE_DISTINCT / "paired_phi.csv")
    require(len(phi) == 4000, "Concrete must have 4,000 coordinate records")
    vectors = {}
    selected = {}
    for row in phi:
        replication = int(row["replication"])
        method = row["method"]
        variable = row["variable_name"]
        require(row["status"] == "complete" and method in METHODS and
                variable in VARIABLES and int(row["seed"]) == 202609220 + replication and
                int(row["variable"]) == VARIABLES.index(variable) + 1,
                "Concrete coordinate identity mismatch")
        value = float(row["estimate"])
        require(math.isfinite(value) and -1e-12 <= value <= 1 + 1e-12,
                "Concrete allocation outside [0,1]")
        require(row["selected"] in {"TRUE", "FALSE"}, "Concrete selection flag invalid")
        key = (replication, method)
        vectors.setdefault(key, {})[variable] = value
        selected[(replication, method, variable)] = row["selected"] == "TRUE"
    require(set(vectors) == run_keys and all(len(v) == 8 for v in vectors.values()),
            "Concrete coordinate vectors are incomplete")
    require(all(abs(sum(v.values()) - 1) < 1e-10 for v in vectors.values()),
            "Concrete allocation vectors do not sum to one")

    coordinates = rows(CONCRETE_DISTINCT / "coordinate_summary.csv")
    require(len(coordinates) == 40, "Concrete coordinate summary must have 40 rows")
    coordinate_map = {(r["method"], r["variable_name"]): r for r in coordinates}
    require(set(coordinate_map) == {(m, v) for m in METHODS for v in VARIABLES},
            "Concrete coordinate summary has missing or duplicate rows")
    for (method, variable), row in coordinate_map.items():
        values = [vectors[(r, method)][variable] for r in range(1, 101)]
        require(int(row["n_complete"]) == 100, "Concrete summary replication count mismatch")
        for key, actual in (
            ("mean", statistics.mean(values)), ("median", quantile(values, 0.5)),
            ("q05", quantile(values, 0.05)), ("q95", quantile(values, 0.95)),
            ("selection_frequency", sum(selected[(r, method, variable)]
                                        for r in range(1, 101)) / 100),
        ):
            close(row[key], actual, f"Concrete summary mismatch: {method} {variable} {key}")

    distances = rows(CONCRETE_DISTINCT / "paired_half_l1.csv")
    require(len(distances) == 1000, "Concrete must have ten paired distances per sample")
    distance_map = {}
    for row in distances:
        key = (int(row["replication"]), row["method_a"], row["method_b"])
        a, b = vectors[(key[0], key[1])], vectors[(key[0], key[2])]
        expected = sum(abs(a[v] - b[v]) for v in VARIABLES) / 2
        close(row["half_l1"], expected, f"Concrete paired distance mismatch: {key}")
        require(key not in distance_map, f"Duplicate Concrete paired distance: {key}")
        distance_map[key] = expected
    require(len(distance_map) == 1000, "Concrete paired distances are not unique")
    pair_summary = rows(CONCRETE_DISTINCT / "pair_summary.csv")
    require(len(pair_summary) == 10, "Concrete pair summary must have ten rows")
    for row in pair_summary:
        values = [distance_map[(r, row["method_a"], row["method_b"])]
                  for r in range(1, 101)]
        require(int(row["n"]) == 100, "Concrete pair replication count mismatch")
        for key, probability in (("median", 0.5), ("q05", 0.05), ("q95", 0.95)):
            close(row[key], quantile(values, probability),
                  f"Concrete pair summary mismatch: {row['method_a']} {row['method_b']} {key}")

    bootstrap = json.loads((CONCRETE_BOOTSTRAP / "VALIDATION.json").read_text())
    require(bootstrap["repetitions"] == 200 and
            bootstrap["complete_method_calls"] == bootstrap["expected_method_calls"] == 1000 and
            bootstrap["all_method_statuses_complete"],
            "Concrete historical bootstrap sensitivity receipt failed")
    old_coordinates = rows(CONCRETE_BOOTSTRAP / "coordinate_summary.csv")
    require(len(old_coordinates) == 40 and
            all(int(r["n"]) == 200 for r in old_coordinates),
            "Concrete bootstrap coordinate summary invalid")
    full_data = rows(CONCRETE_BOOTSTRAP / "full_data_phi.csv")
    require(len(full_data) == 40, "Concrete full-data reference must have 40 rows")
    neighbor = json.loads((CONCRETE_DISTINCT / "neighbor_full_findings.json").read_text())
    fractions = neighbor["results"]["full_8"]["metrics"]["nearest_response_tie_rate"]
    close(fractions["bootstrap_mean"], 0.6681, "Concrete bootstrap 2NN diagnostic changed")
    close(fractions["distinct_80pct_mean"], 0.1354,
          "Concrete distinct-row 2NN diagnostic changed")
    return {
        "source_rows": 1030, "coordinates": 8, "distinct_rows_per_sample": 824,
        "subsamples": 100, "complete_method_fits": 500,
        "historical_bootstraps": 200, "historical_complete_method_fits": 1000,
    }


def validate_timing_parity():
    receipt = json.loads((TIMING_PARITY / "RECEIPT.json").read_text())
    require(receipt["status"] == "PASS" and receipt["complete"] == 150 and
            receipt["replications"] == 30 and receipt["maximum_effect_abs_gap"] == 0 and
            receipt["selected_mismatches"] == receipt["total_warnings"] == 0,
            "Friedman timing parity receipt failed")
    require(sha256(TIMING_PARITY / "timings.csv") == receipt["timings_sha256"] and
            sha256(TIMING_PARITY / "source_manifest.csv") == receipt["source_manifest_sha256"],
            "Friedman timing parity source hashes changed")
    timings = rows(TIMING_PARITY / "timings.csv")
    summary = rows(TIMING_PARITY / "summary.csv")
    require(len(timings) == 150 and len(summary) == 5,
            "Friedman timing parity rows incomplete")
    for row in summary:
        subset = [float(t["elapsed_sec"]) for t in timings if t["method"] == row["method"]]
        require(len(subset) == int(row["n"]) == 30,
                f"Friedman timing parity method incomplete: {row['method']}")
        for key, probability in (("median_sec", 0.5), ("q1_sec", 0.25), ("q3_sec", 0.75)):
            close(row[key], quantile(subset, probability),
                  f"Friedman timing parity summary mismatch: {row['method']} {key}")
    return {"datasets": 30, "complete_method_fits": 150, "serial_workers": 1}


def validate_manuscript_values(case):
    main_text = "\n".join(
        path.read_text() for path in [
            MANUSCRIPT / "sections/new_friedman.tex",
            MANUSCRIPT / f"sections/new_{case}.tex",
            MANUSCRIPT / "sections/07_conclusion.tex",
            MANUSCRIPT / "tables/pme_accuracy.tex",
        ]
    )
    for forbidden in [
        "common game", "common-game", "exponential stage",
        "The numerical results separate these roles",
        "engineering-response screens",
        "allocation game are held fixed",
    ]:
        require(forbidden not in main_text, f"Stale manuscript phrase: {forbidden}")
    for required in [
        "Full LC-\\(K_n\\)",
        "local sample variances",
        "nearest other observation",
        "FIRST-SHAP targets Shapley effects",
    ]:
        require(required in main_text, f"Missing manuscript definition: {required}")

    combined = rows(EVIDENCE / "matched_k/pme_accuracy_four_methods.csv")
    table = (MANUSCRIPT / "tables/pme_accuracy.tex").read_text()
    for row in combined:
        token = f"{float(row['mean']):.3f}"
        require(token in table, f"Missing rounded PME table value: {token}")
    central = rows(EVIDENCE / "matched_k/central_estimator_sensitivity.csv")
    supplement = (MANUSCRIPT / "sections/supplement_results.tex").read_text()
    for row in central:
        token = f"{float(row['mean_error']):.5f}"
        require(token in supplement,
                f"Missing estimator-sensitivity value: {token}")
    friedman = (MANUSCRIPT / "sections/new_friedman.tex").read_text()
    for row in rows(TIMING_PARITY / "summary.csv"):
        display = {
            "PRIME": "PRIME", "Full_LL_PME": "Full LL-PME",
            "Full_LC_Kn": r"Full LC-\(K_n\)", "PME_kNN": "PME-kNN",
            "FIRST_SHAP": "FIRST-SHAP",
        }[row["method"]]
        line = next((line for line in friedman.splitlines()
                     if line.startswith(display + " &")), None)
        require(line is not None, f"Friedman timing row missing: {display}")
        timing_matches = re.findall(r" & ([0-9.]+) \[([0-9.]+),([0-9.]+)\] & ", line)
        require(timing_matches, f"Friedman timing format invalid: {display}")
        for actual, key in zip(timing_matches[-1], ("median_sec", "q1_sec", "q3_sec")):
            close(actual, row[key], f"Friedman timing value mismatch: {display} {key}",
                  tolerance=0.00051)
    if case == "concrete":
        concrete = (MANUSCRIPT / "sections/new_concrete.tex").read_text()
        supplement_case = (MANUSCRIPT / "sections/supplement_concrete_case.tex").read_text()
        summary = {(r["method"], r["variable_name"]): r for r in
                   rows(CONCRETE_DISTINCT / "coordinate_summary.csv")}
        for method, variable in (
            ("PRIME", "water_binder_ratio"),
            ("Full LL-PME", "water_binder_ratio"),
            ("PME-kNN", "water_binder_ratio"),
            ("PRIME", "cement"), ("Full LL-PME", "cement"),
        ):
            token = f"{float(summary[(method, variable)]['median']):.3f}"
            require(token in concrete, f"Missing Concrete rounded allocation: {method} {variable}")
        pair = next(r for r in rows(CONCRETE_DISTINCT / "pair_summary.csv")
                    if r["method_a"] == "PRIME" and r["method_b"] == "Full LL-PME")
        require(f"{float(pair['median']):.3f}" in concrete,
                "Missing Concrete paired-distance claim")
        require("100 subsamples of 824 distinct rows" in concrete and
                "evaluations is 16" in supplement_case and
                "versus 256" in supplement_case and
                "200 size-1,030 bootstrap samples" in supplement_case and
                "66.8\\%" in supplement_case and "13.5\\%" in supplement_case,
                "Concrete design or sensitivity claim missing")
    else:
        auto_text = (MANUSCRIPT / "sections/new_auto_mpg.tex").read_text()
        require("392 complete cars" in auto_text and
                "auto_mpg_p8_fullbars_stability.pdf" in auto_text,
                "Auto MPG main-case claims missing")


def release_files():
    relative = [
        "README.md", "Makefile", "METHOD_SPEC.md",
        "manuscript/main.tex", "manuscript/supplement.tex",
        "manuscript/references.bib", "manuscript/main_labels.aux",
        "manuscript/sections/00_abstract.tex",
        "manuscript/sections/01_introduction.tex",
        "manuscript/sections/02_allocation.tex",
        "manuscript/sections/03_estimation.tex",
        "manuscript/sections/04_first.tex",
        "manuscript/sections/05_framework.tex",
        "manuscript/sections/06_experiments.tex",
        "manuscript/sections/07_conclusion.tex",
        "manuscript/sections/after_selection.tex",
        "manuscript/sections/observed_case_config.tex",
        "manuscript/sections/new_auto_mpg.tex",
        "manuscript/sections/new_concrete.tex",
        "manuscript/sections/new_friedman.tex",
        "manuscript/sections/supplement_auto_mpg_case.tex",
        "manuscript/sections/supplement_concrete_case.tex",
        "manuscript/sections/supplement_experiments.tex",
        "manuscript/sections/supplement_proofs.tex",
        "manuscript/sections/supplement_results.tex",
        "manuscript/tables/observed_application_main.tex",
        "manuscript/tables/pme_accuracy.tex",
        "manuscript/tables/pme_correlated_pipeline.tex",
        "manuscript/tables/pme_scaling_pmax.tex",
        "manuscript/tables/projection_sensitivity.tex",
        "manuscript/tables/ridge_seed_diagnostic.tex",
        "manuscript/tables/twinning_diagnostic.tex",
    ]
    relative += [
        "manuscript/figures/" + name for name in [
            "auto_mpg_data_context.pdf",
            "auto_mpg_p8_fullbars_stability.pdf",
            "concrete_p8_distinct_subsamples.pdf",
            "conditional_variance_rmse.pdf",
            "exact_recovery_supplement.pdf",
            "friedman_ll_allocation_rules.pdf",
            "friedman_ll_dependence_path.pdf",
            "friedman_methods.pdf",
            "friedman_selection.pdf",
            "local_linear_drift.pdf",
            "observed_application_main.pdf",
            "other_functions_selection.pdf",
            "pme_estimator_calibration.pdf",
            "pme_scaling_main.pdf",
            "selection_to_pme_combined.pdf",
            "shapley_pme_proxy_population.pdf",
        ]
    ]
    return relative


def active_files(case):
    inactive = {
        "concrete": {
            "manuscript/sections/new_auto_mpg.tex",
            "manuscript/sections/supplement_auto_mpg_case.tex",
        },
        "auto_mpg": {
            "manuscript/sections/new_concrete.tex",
            "manuscript/sections/supplement_concrete_case.tex",
        },
    }[case]
    return [path for path in release_files() if path not in inactive]


def validate_tex_closure(case):
    for relative in release_files():
        require((ROOT / relative).is_file(), f"Missing active file: {relative}")
    tex_files = [ROOT / path for path in active_files(case) if path.endswith(".tex")]
    text = "\n".join(path.read_text() for path in tex_files)
    labels = re.findall(r"\\label\{([^}]+)\}", text)
    refs = re.findall(r"\\(?:eqref|ref|autoref)\{([^}]+)\}", text)
    duplicates = [key for key, count in Counter(labels).items() if count > 1]
    require(not duplicates, f"Duplicate labels: {duplicates}")
    require(not (set(refs) - set(labels)),
            f"Missing labels: {sorted(set(refs) - set(labels))}")
    return len(labels)


def build_pdfs(case):
    build_env = os.environ.copy()
    build_env["SOURCE_DATE_EPOCH"] = "1790035200"  # 2026-09-22 00:00 UTC
    build_env["FORCE_SOURCE_DATE"] = "1"
    subprocess.run(
        ["latexmk", "-g", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
         "main.tex"],
        cwd=MANUSCRIPT, env=build_env, check=True, stdout=subprocess.DEVNULL
    )
    aux = (MANUSCRIPT / "main.aux").read_text().splitlines()
    labels = ["\\relax"] + [line for line in aux if line.startswith("\\newlabel")]
    (MANUSCRIPT / "main_labels.aux").write_text("\n".join(labels) + "\n")
    subprocess.run(
        ["latexmk", "-g", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
         "supplement.tex"],
        cwd=MANUSCRIPT, env=build_env, check=True, stdout=subprocess.DEVNULL
    )
    bad = re.compile(
        r"^(LaTeX Warning|Package .* Warning|pdfTeX warning|Overfull|Underfull|!)"
    )
    for name in ["main", "supplement"]:
        log = (MANUSCRIPT / f"{name}.log").read_text().splitlines()
        warnings = [line for line in log if bad.match(line)]
        require(not warnings, f"{name} build warnings: {warnings[:5]}")
    pages = {}
    expected_pages = {
        "concrete": {"main": 28, "supplement": 28},
        "auto_mpg": {"main": 30, "supplement": 27},
    }[case]
    for name, expected in expected_pages.items():
        output = subprocess.check_output(
            ["pdfinfo", str(MANUSCRIPT / f"{name}.pdf")], text=True
        )
        match = re.search(r"^Pages:\s+(\d+)", output, re.MULTILINE)
        require(match and int(match.group(1)) == expected,
                f"Unexpected {name} page count")
        pages[name] = int(match.group(1))
    return pages


def validate_case_variants(active_case):
    """Compile both cases and leave the checked-out case and PDFs active."""
    validate_manuscript_values(active_case)
    active_labels = validate_tex_closure(active_case)
    active_pages = build_pdfs(active_case)
    alternative = "auto_mpg" if active_case == "concrete" else "concrete"
    selector_path = MANUSCRIPT / "sections/observed_case_config.tex"
    original = selector_path.read_text()
    replacement = r"\ConcreteMainfalse" if alternative == "auto_mpg" else r"\ConcreteMaintrue"
    switched = re.sub(r"(?m)^\\ConcreteMain(?:true|false)$", replacement.replace("\\", r"\\"),
                      original)
    require(switched != original, "Case selector could not be toggled")
    try:
        selector_path.write_text(switched)
        require(selected_case() == alternative, "Alternate case selector failed")
        validate_manuscript_values(alternative)
        alternative_labels = validate_tex_closure(alternative)
        alternative_pages = build_pdfs(alternative)
    finally:
        selector_path.write_text(original)
        build_pdfs(active_case)
    require(selected_case() == active_case, "Active case was not restored")
    return {
        active_case: {"pages": active_pages, "labels": active_labels},
        alternative: {"pages": alternative_pages, "labels": alternative_labels},
    }


def write_manifest(validation):
    evidence_files = sorted(
        path for path in EVIDENCE.rglob("*")
        if path.is_file() and path not in {
            EVIDENCE / "manifest.json", EVIDENCE / "VALIDATION.json"
        }
    )
    source_files = [
        ROOT / "experiments/12_submission_update/01_run_matched_k.R",
        ROOT / "experiments/12_submission_update/matched_k_helpers.R",
        ROOT / "experiments/12_submission_update/02_validate_submission.py",
        ROOT / "experiments/13_concrete_submission/build_concrete_p8_figure.py",
        ROOT / "experiments/13_concrete_submission/first_adapter.R",
    ]
    files = [ROOT / path for path in release_files()] + evidence_files + source_files
    manifest = {
        "release_id": "20260920_submission",
        "status": "validated_submission_overlay",
        "active_case": validation["active_case"],
        "files": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in sorted(files)
        },
    }
    (EVIDENCE / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    (EVIDENCE / "VALIDATION.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )


def main():
    case = selected_case()
    legacy = validate_legacy_inputs()
    matched = validate_matched_k()
    auto = validate_auto_mpg()
    concrete = validate_concrete()
    timing = validate_timing_parity()
    variants = validate_case_variants(case)
    labels = variants[case]["labels"]
    pages = variants[case]["pages"]
    validation = {
        "status": "PASS",
        "active_case": case,
        "legacy_inputs_verified": legacy,
        "matched_k": matched,
        "auto_mpg": auto,
        "concrete": concrete,
        "friedman_timing_parity": timing,
        "pages": pages,
        "case_variants": variants,
        "labels": labels,
        "latex_warnings": 0,
    }
    write_manifest(validation)
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
