#!/usr/bin/env python3
"""Run existing paper scripts against a verified release in an isolated cwd."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tempfile
import time

PROTOCOL = "prime-verified-workflow-v2"
BOOTSTRAP_RECEIPT = ".verified/local_install.json"
BOOTSTRAP_DEPENDENCIES = ["FNN", "MASS", "Rcpp", "RcppArmadillo", "testthat",
                          "twinning", "ggplot2", "jsonlite"]
BASELINES = ["selection_competitors", "selection_to_pme_competitors",
             "selection_to_pme_competitors_phi", "pme_controls", "pme_controls_phi",
             "correlated_competitors", "correlated_competitors_phi",
             "scaling_competitors", "scaling_competitors_phi"]
ARTIFACT_SCRIPTS = [f"experiments/90_artifacts/{name}" for name in
    ["01_build_main.R", "02_build_supplement.R", "03_build_concrete.R",
     "04_build_observed.R", "05_build_seed_diagnostic.R"]]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def relative(value):
    p = PurePosixPath(value)
    if not value or p.is_absolute() or ".." in p.parts or "\\" in value:
        raise ValueError(f"Expected a safe relative path: {value!r}")
    return Path(*p.parts)


def confined(root, value):
    path = root / relative(value)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes its declared root: {value}")
    return path


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def logical_input(source_root, record):
    """Open logical bytes; transport hashes are separate from scientific hashes."""
    transport = record.get("transport")
    if transport is None:
        return confined(source_root, record["path"]).open("rb")
    if transport.get("encoding") != "gzip":
        raise ValueError("Unsupported release input transport")
    path = confined(source_root, transport["path"])
    if digest(path) != transport["sha256"]:
        raise ValueError(f"Compressed transport hash mismatch: {transport['path']}")
    return gzip.open(path, "rb")


def verify_logical_input(source_root, logical, record):
    value = hashlib.sha256()
    size = 0
    with logical_input(source_root, record) as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
            size += len(chunk)
    if value.hexdigest() != record["sha256"]:
        raise ValueError(f"Release input hash mismatch: {logical}")
    expected_size = record.get("transport", {}).get("uncompressed_bytes")
    if expected_size is not None and size != expected_size:
        raise ValueError(f"Release input size mismatch: {logical}")


def package_source_inventory(root):
    """Every package build input is pinned; stale compiled objects are excluded."""
    source = root / "software/prime"
    return {str(path.relative_to(root)) for path in source.rglob("*")
            if path.is_file() and path.suffix not in (".o", ".so", ".dll", ".a")
            and path.name != ".DS_Store"}


MAINTAINED_PACKAGE_ROOTS = {"r": "packages/r", "python": "packages/python"}
MAINTAINED_TRANSIENT_DIRECTORIES = frozenset({
    "__pycache__", ".pytest_cache", ".ruff_cache", "build", "dist",
    "library", "tmp", "check"})
MAINTAINED_COMPILED_SUFFIXES = frozenset({".o", ".so", ".dll", ".a"})


def maintained_relative(value):
    """New inventory paths have one canonical spelling; old pins are unchanged."""
    if not isinstance(value, str):
        raise ValueError("Maintained package paths must be strings")
    path = relative(value)
    if path.as_posix() != value or value == ".":
        raise ValueError(f"Noncanonical maintained package path: {value!r}")
    return path


def maintained_package_source_inventory(root, package_root):
    """Pin every maintained source/build input, excluding named transient artifacts.

    Paths are repository-relative. Source symlinks and special files are refused,
    including symlinked package ancestors. Ignored directories are not traversed.
    This is deliberately separate from the frozen software/prime inventory.
    """
    root = Path(root)
    relative_root = maintained_relative(package_root)
    if package_root not in MAINTAINED_PACKAGE_ROOTS.values():
        raise ValueError(f"Unsupported maintained package root: {package_root}")
    source = confined(root, package_root)
    ancestor = root
    for component in relative_root.parts:
        ancestor = ancestor / component
        if ancestor.is_symlink():
            raise ValueError(f"Maintained package source symlink: {ancestor.relative_to(root)}")
    if not source.is_dir():
        raise ValueError(f"Missing maintained package root: {package_root}")
    inventory = set()

    def walk_error(error):
        raise ValueError(f"Cannot inspect maintained package source: {error}") from error

    for directory, directories, filenames in os.walk(source, followlinks=False, onerror=walk_error):
        directory = Path(directory)
        kept_directories = []
        for name in sorted(directories):
            path = directory / name
            if path.is_symlink():
                raise ValueError(f"Maintained package source symlink: {path.relative_to(root)}")
            if (name in MAINTAINED_TRANSIENT_DIRECTORIES or
                    name.endswith((".egg-info", ".Rcheck"))):
                continue
            kept_directories.append(name)
        directories[:] = kept_directories
        for name in sorted(filenames):
            path = directory / name
            if path.is_symlink():
                raise ValueError(f"Maintained package source symlink: {path.relative_to(root)}")
            if (path.suffix in MAINTAINED_COMPILED_SUFFIXES or name == ".DS_Store" or
                    name.endswith((".egg-info", ".Rcheck"))):
                continue
            if not path.is_file():
                raise ValueError(f"Maintained package source is not a regular file: {path.relative_to(root)}")
            inventory.add(path.relative_to(root).as_posix())
    return inventory


def maintained_sha256(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("Maintained package SHA256 must be 64 lowercase hexadecimal characters")
    return value


def verify_maintained_package_version(source, language, expected):
    """Read static source metadata without importing or installing either package."""
    if not isinstance(expected, str) or not expected:
        raise ValueError("Maintained package version must be a nonempty string")
    if language == "r":
        fields = {name: [] for name in ("Package", "Version")}
        current_field = None
        for line in (source / "DESCRIPTION").read_text().splitlines():
            if line.startswith((" ", "\t")) and current_field in fields:
                raise ValueError("Maintained R package name/version must be single-line fields")
            current_field = line.split(":", 1)[0] if ":" in line and not line[0].isspace() else None
            for name in fields:
                if line.startswith(name + ":"):
                    fields[name].append(line.split(":", 1)[1].strip())
        if fields["Package"] != ["prime"] or fields["Version"] != [expected]:
            raise ValueError("Maintained R package name/version mismatch")
    else:
        # Existing pyproject metadata uses static, single-line name/version
        # strings. Keep this guard dependency-free on Python 3.10 as well.
        text = (source / "pyproject.toml").read_text()
        # Refuse forms the bounded Python 3.10 parser cannot interpret: a fake
        # table inside a multiline string must never stand in for metadata.
        if "\"\"\"" in text or "'''" in text:
            raise ValueError("Maintained Python package metadata must use single-line strings")
        for header in re.finditer(r"(?m)^[ \t]*\[(.*)\][ \t]*(?:#.*)?$", text):
            if re.fullmatch(r"[A-Za-z0-9_-]+(?:\.[A-Za-z0-9_-]+)*", header.group(1)) is None:
                raise ValueError("Maintained Python package metadata must use bare table names")
        sections = list(re.finditer(r"(?m)^[ \t]*\[([^\]\r\n]+)\][ \t]*(?:#.*)?$", text))
        projects = [i for i, section in enumerate(sections) if section.group(1) == "project"]
        if len(projects) != 1:
            raise ValueError("Maintained Python package requires one static [project] table")
        index = projects[0]
        end = sections[index + 1].start() if index + 1 < len(sections) else len(text)
        project = text[sections[index].end():end]
        fields = {name: [] for name in ("name", "version")}
        for line in project.splitlines():
            if re.match(r"[ \t]*dynamic[ \t]*=", line):
                raise ValueError("Maintained Python package metadata must be static")
            if re.match(r"[ \t]*(['\"])[^'\"\r\n]+\1[ \t]*(?:=|\.)", line):
                raise ValueError("Maintained Python package metadata must use unquoted keys")
            match = re.fullmatch(r"[ \t]*(name|version)[ \t]*=[ \t]*(['\"])([^'\"\\\r\n]+)\2[ \t]*(?:#.*)?", line)
            if match:
                fields[match.group(1)].append(match.group(3))
            elif re.match(r"[ \t]*(name|version)[ \t]*(?:=|\.)", line):
                raise ValueError("Maintained Python package name/version must be static strings")
        if fields["name"] != ["prime-importance"] or fields["version"] != [expected]:
            raise ValueError("Maintained Python package name/version mismatch")
        # The standard parser additionally checks the complete TOML on Python
        # 3.11+. The bounded static format above remains usable on Python 3.10.
        try:
            import tomllib
        except ModuleNotFoundError:
            pass
        else:
            try:
                parsed = tomllib.loads(text)
            except tomllib.TOMLDecodeError as error:
                raise ValueError("Invalid maintained Python package metadata") from error
            if (parsed.get("project", {}).get("name") != "prime-importance" or
                    parsed.get("project", {}).get("version") != expected):
                raise ValueError("Maintained Python package name/version mismatch")


def maintained_json_object(pairs):
    """Refuse ambiguous duplicate keys in the new independently pinned manifest."""
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate maintained package manifest key: {key}")
        value[key] = item
    return value


def verify_maintained_packages(root, reference):
    """Verify optional maintained packages without changing the scientific manifest."""
    if not isinstance(reference, dict) or set(reference) != {"manifest", "sha256"}:
        raise ValueError("Unsupported maintained package manifest reference")
    maintained_relative(reference["manifest"])
    manifest_path = confined(root, reference["manifest"])
    expected_digest = maintained_sha256(reference["sha256"])
    manifest_bytes = manifest_path.read_bytes()
    if hashlib.sha256(manifest_bytes).hexdigest() != expected_digest:
        raise ValueError("Maintained package manifest hash mismatch")
    manifest = json.loads(manifest_bytes, object_pairs_hook=maintained_json_object)
    if (not isinstance(manifest, dict) or set(manifest) != {"schema_version", "packages"} or
            type(manifest["schema_version"]) is not int or manifest["schema_version"] != 1 or
            not isinstance(manifest["packages"], dict) or
            set(manifest["packages"]) != set(MAINTAINED_PACKAGE_ROOTS)):
        raise ValueError("Unsupported maintained package manifest schema")
    for language, expected_root in MAINTAINED_PACKAGE_ROOTS.items():
        package = manifest["packages"][language]
        if (not isinstance(package, dict) or set(package) != {"root", "version", "source_files"} or
                package["root"] != expected_root):
            raise ValueError(f"Unsupported maintained {language} package record/root")
        source_files = package["source_files"]
        if not isinstance(source_files, dict) or not source_files:
            raise ValueError("Maintained package must pin source files")
        for logical, expected in source_files.items():
            logical_path = maintained_relative(logical)
            if not logical_path.is_relative_to(maintained_relative(expected_root)):
                raise ValueError(f"Maintained source path escapes its package root: {logical}")
            maintained_sha256(expected)
        if maintained_package_source_inventory(root, expected_root) != set(source_files):
            raise ValueError(f"Maintained package source inventory differs from its pins: {expected_root}")
        for logical, expected in source_files.items():
            if digest(confined(root, logical)) != expected:
                raise ValueError(f"Maintained package source hash mismatch: {logical}")
        verify_maintained_package_version(confined(root, expected_root), language, package["version"])


def verify_release(root):
    manifest_path = root / "results/current_release.json"
    historical = json.loads(manifest_path.read_text())
    pointer_path = root / "results/current_workflow.json"
    maintained_reference = None
    maintained_present = False
    if pointer_path.exists():
        pointer = json.loads(pointer_path.read_text())
        if pointer.get("schema_version") != 1:
            raise ValueError("Unsupported workflow pointer schema")
        if pointer.get("historical_release_sha256") != digest(manifest_path):
            raise ValueError("Historical release manifest changed")
        revised_path = confined(root, pointer["manifest"])
        if digest(revised_path) != pointer["sha256"]:
            raise ValueError("Workflow manifest hash mismatch")
        revised = json.loads(revised_path.read_text())
        # This maintenance boundary changes orchestration/documentation only.
        # Scientific inputs, implementation, and experiment source stay pinned.
        # Only the named documentation and workflow files may be maintained.
        audited_paths = {"Makefile", "README.md", "CITATION.cff", "VALIDATION.md", ".gitignore",
                         "software/prime/README.md",
                         "software/prime/_pkgdown.yml",
                         "software/prime/man/prime-package.Rd",
                         "software/prime/man/prime.Rd",
                         "software/prime/man/prime_allocate.Rd",
                         "software/prime/man/prime_screen.Rd",
                         "software/prime/man/summary.prime_fit.Rd",
                         "software/prime/man/summary.prime_screen.Rd",
                         ".github/workflows/R-CMD-check.yaml",
                         "experiments/99_validate/verified_workflow.py",
                         "experiments/99_validate/test_verified_workflow.py"}
        if set(revised) != set(historical):
            raise ValueError("Workflow revision changed manifest fields")
        if any(revised[k] != historical[k] for k in historical if k != "source_files"):
            raise ValueError("Workflow revision changed scientific release metadata")
        if set(revised["source_files"]) != set(historical["source_files"]):
            raise ValueError("Workflow revision changed source inventory")
        changed = {p for p in historical["source_files"]
                   if revised["source_files"][p] != historical["source_files"][p]}
        if not changed <= audited_paths:
            raise ValueError("Workflow revision changed scientific source pins")
        if "maintained_packages" in pointer:
            maintained_present = True
            maintained_reference = pointer["maintained_packages"]
        manifest_path = revised_path
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 2:
        raise ValueError("Unsupported current_release schema")
    if manifest.get("status") != "validated_public_release":
        raise ValueError("Release inputs have not passed publication validation")
    if not manifest.get("source_files") or not manifest.get("files"):
        raise ValueError("Release must pin source files and input files")
    data_root = confined(root, manifest["data_root"])
    relative(manifest["release_id"])
    for logical, record in manifest["files"].items():
        relative(logical)
        if logical.startswith("results/raw/") and logical.lower().endswith(".rds"):
            raise ValueError("Historical raw checkpoints must not be hydrated")
        verify_logical_input(data_root, logical, record)
    for logical, expected in manifest["source_files"].items():
        if logical in manifest["files"]:
            raise ValueError(f"Source/input destination collision: {logical}")
        if digest(confined(root, logical)) != expected:
            raise ValueError(f"Canonical source hash mismatch: {logical}")
    package_pins = {p for p in manifest["source_files"] if p.startswith("software/prime/")}
    if package_source_inventory(root) != package_pins:
        raise ValueError("Package build-source inventory differs from the release pins")
    lineage = manifest["historical_lineage"]
    if digest(confined(root, lineage["manifest"])) != lineage["sha256"]:
        raise ValueError("Historical lineage manifest changed")
    for name in BASELINES:
        logical = f"results/provenance/paired_baselines/{name}.csv.gz"
        if logical not in manifest["files"]:
            raise ValueError(f"Missing corrected baseline: {logical}")
    if maintained_present:
        verify_maintained_packages(root, maintained_reference)
    return manifest, digest(manifest_path)


def installed_snapshot(package, version):
    description = (package / "DESCRIPTION").read_text()
    actual = next((line.split(":", 1)[1].strip() for line in description.splitlines()
                   if line.startswith("Version:")), None)
    if actual != version:
        raise ValueError(f"Installed prime version {actual!r} differs from {version!r}")
    files = {str(p.relative_to(package)): digest(p)
             for p in sorted(package.rglob("*")) if p.is_file()}
    if not files:
        raise ValueError("Installed prime package is empty")
    return package, files


def r_environment(root, library=None):
    environment = os.environ.copy()
    # R_LIBS precedes R_LIBS_USER in R's search path. Do not let an inherited
    # library silently select a different PRIME than the attested installation.
    environment.pop("R_LIBS", None)
    environment["R_LIBS_USER"] = os.pathsep.join(dict.fromkeys(filter(None, [
        str(library or root / "r-library"), str(root / "r-library"),
        environment.get("R_LIBS_USER", ""), str(Path.home() / "R/library")])))
    return environment


def runtime_snapshot(root):
    names = ",".join(json.dumps(x) for x in BOOTSTRAP_DEPENDENCIES)
    expression = f'''packages <- c({names})
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop(paste("Missing R dependencies:", paste(missing,collapse=", "), "; run make deps"))
versions <- setNames(vapply(packages,function(p)as.character(packageVersion(p)),character(1)),packages)
cat(jsonlite::toJSON(list(R_version=R.version.string, platform=R.version$platform,
  dependencies=as.list(versions)),auto_unbox=TRUE))'''
    result = subprocess.run(["Rscript", "--vanilla", "-e", expression],
                            env=r_environment(root), capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or "Cannot inspect the R runtime")
    return json.loads(result.stdout)


def package_snapshot(root, manifest, manifest_hash):
    receipt_path = root / BOOTSTRAP_RECEIPT
    if not receipt_path.is_file():
        raise ValueError("No verified local installation; run make bootstrap first")
    receipt = json.loads(receipt_path.read_text())
    if (receipt.get("schema_version") != 1 or
            receipt.get("release_sha256") != manifest_hash or
            receipt.get("source_files") != manifest["source_files"] or
            receipt.get("input_sha256") != {k: v["sha256"] for k, v in manifest["files"].items()}):
        raise ValueError("Local bootstrap receipt is stale; verify the release and run make bootstrap")
    if receipt.get("runtime") != runtime_snapshot(root):
        raise ValueError("R/dependency versions changed since bootstrap; run make bootstrap")
    package, files = installed_snapshot(confined(root, receipt["package_path"]),
                                        manifest["package_version"])
    if files != receipt.get("installed_files") or receipt.get("smoke_validation") != "passed":
        raise ValueError("Installed PRIME files differ from the verified local bootstrap")
    return package, files


def bootstrap(root):
    """Explicitly compile reviewed source; never update release/source/data pins."""
    manifest, manifest_hash = verify_release(root)
    runtime = runtime_snapshot(root)
    work_parent = root / "tmp/bootstrap"
    work_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="fresh-", dir=work_parent))
    source = work / "source/prime"
    library = work / "library"
    library.mkdir()
    package_pins = {p: h for p, h in manifest["source_files"].items()
                    if p.startswith("software/prime/")}
    for logical, expected in package_pins.items():
        target = source / relative(logical).relative_to("software/prime")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(confined(root, logical), target)
        if digest(target) != expected:
            raise ValueError(f"Source changed during bootstrap copy: {logical}")
    command = ["R", "CMD", "INSTALL", "--preclean", "--no-multiarch",
               "-l", str(library), str(source)]
    with (work / "install.log").open("w") as log:
        result = subprocess.run(command, cwd=work, env=r_environment(root, library),
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Fresh package installation failed; see {work / 'install.log'}")
    _, files = installed_snapshot(library / "prime", manifest["package_version"])
    _, after_hash = verify_release(root)
    if after_hash != manifest_hash:
        raise ValueError("Release changed during bootstrap")
    destination = root / "r-library/prime"
    destination.parent.mkdir(parents=True, exist_ok=True)
    previous = work / "previous-prime"
    if destination.exists():
        destination.rename(previous)
    (library / "prime").rename(destination)
    try:
        # The bounded source-pinned smoke runs from this fresh build directory,
        # so its generated summary cannot overwrite reviewed paper outputs.
        smoke = "experiments/99_validate/01_smoke.R"
        if smoke not in manifest["source_files"]:
            raise ValueError("The bootstrap smoke must be source-pinned")
        smoke_copy = work / "smoke.R"
        shutil.copy2(confined(root, smoke), smoke_copy)
        expression = (f'stopifnot(normalizePath(find.package("prime")) == '
                      f'normalizePath({json.dumps(str(destination))})); source("smoke.R")')
        with (work / "smoke.log").open("w") as log:
            result = subprocess.run(["Rscript", "--vanilla", "-e", expression],
                                    cwd=work, env=r_environment(root),
                                    stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Fresh-install smoke failed; see {work / 'smoke.log'}")
        _, checked_files = installed_snapshot(destination, manifest["package_version"])
        _, checked_hash = verify_release(root)
        if checked_files != files or checked_hash != manifest_hash or runtime_snapshot(root) != runtime:
            raise ValueError("Source, inputs, runtime, or installation changed during bootstrap")
    except Exception:
        destination.rename(work / "failed-prime")
        if previous.exists():
            previous.rename(destination)
        raise
    receipt = {"schema_version": 1, "kind": "explicit_fresh_build_from_pinned_source",
               "release_sha256": manifest_hash, "source_files": manifest["source_files"],
               "input_sha256": {k: v["sha256"] for k, v in manifest["files"].items()},
               "package_version": manifest["package_version"], "package_path": "r-library/prime",
               "installed_files": files, "runtime": runtime, "smoke_validation": "passed",
               "created_unix": time.time(), "build_directory": str(work.relative_to(root)),
               "install_log_sha256": digest(work / "install.log"),
               "smoke_log_sha256": digest(work / "smoke.log"),
               "build_environment": {k: os.environ.get(k, "") for k in
                   ["R_MAKEVARS_USER", "CC", "CXX", "CXXFLAGS", "CPPFLAGS", "LDFLAGS"]}}
    write_json(root / BOOTSTRAP_RECEIPT, receipt)
    print(json.dumps({"bootstrap": "verified", "receipt": BOOTSTRAP_RECEIPT,
                      "release_sha256": manifest_hash, "package_version": manifest["package_version"]}))


def check_baseline_vectors(work):
    """Require a complete, finite normalized vector for every feasible baseline."""
    directory = work / "results/provenance/paired_baselines"
    for stem in ["selection_to_pme_competitors", "pme_controls",
                 "correlated_competitors", "scaling_competitors"]:
        with gzip.open(directory / f"{stem}.csv.gz", "rt", newline="") as handle:
            rows = list(csv.DictReader(handle))
        with gzip.open(directory / f"{stem}_phi.csv.gz", "rt", newline="") as handle:
            vectors = list(csv.DictReader(handle))
        if not rows or not vectors:
            raise ValueError(f"Empty baseline: {stem}")
        keys = [k for k in ["benchmark", "design", "n", "p", "rho", "replication", "method"]
                if k in rows[0] and k in vectors[0]]
        if "method" not in keys or "replication" not in keys:
            raise ValueError(f"Missing baseline pairing fields: {stem}")
        groups = {}
        for row in vectors:
            groups.setdefault(tuple(row[k] for k in keys), []).append(row)
        seen = set()
        for row in rows:
            key = tuple(row[k] for k in keys)
            if key in seen:
                raise ValueError(f"Duplicate baseline run key: {stem} {key}")
            seen.add(key)
            if "feasible" in row and row["feasible"] not in ("TRUE", "FALSE"):
                raise ValueError(f"Unknown baseline feasibility: {stem} {key}")
            if row.get("feasible", "TRUE") != "TRUE":
                continue
            vector = groups.get(key, [])
            try:
                values = [float(v["estimate"]) for v in vector]
            except (KeyError, ValueError):
                raise ValueError(f"Missing finite baseline vector: {stem} {key}") from None
            import math
            try:
                dimension = float(row["p"])
                labels = [float(v["variable"]) for v in vector]
                valid_coordinates = (math.isfinite(dimension) and dimension > 0
                    and dimension.is_integer() and len(labels) == int(dimension)
                    and all(math.isfinite(v) and v.is_integer() for v in labels)
                    and sorted(int(v) for v in labels) == list(range(1, int(dimension) + 1)))
            except (KeyError, ValueError, OverflowError):
                valid_coordinates = False
            if not valid_coordinates:
                raise ValueError(f"Baseline vector must contain exactly coordinates 1:p: {stem} {key}")
            if (not values or not all(math.isfinite(x) and x >= -1e-8 for x in values)
                    or abs(sum(values) - 1) > 1e-6):
                raise ValueError(f"Invalid complete baseline vector: {stem} {key}")
        if set(groups) - seen:
            raise ValueError(f"Orphan baseline vector keys: {stem}")


def verify_hydrated_inputs(work, manifest):
    for logical, record in manifest["files"].items():
        if digest(confined(work, logical)) != record["sha256"]:
            raise ValueError(f"Hydrated release input hash mismatch: {logical}")


def prepare(root, namespace):
    namespace = str(relative(namespace))
    if "/" in namespace:
        raise ValueError("Namespace must be a single component")
    manifest, manifest_hash = verify_release(root)
    required_sources = ["experiments/99_validate/verified_workflow.py",
                        "experiments/lib/verified_io.R"]
    if not all(path in manifest["source_files"] for path in required_sources):
        raise ValueError("Release must pin the workflow and IO helpers")
    package, package_files = package_snapshot(root, manifest, manifest_hash)
    pin = {"protocol": PROTOCOL, "namespace": namespace,
           "release_sha256": manifest_hash, "source_files": manifest["source_files"],
           "package_version": manifest["package_version"], "package_files": package_files,
           "bootstrap_receipt_sha256": digest(root / BOOTSTRAP_RECEIPT)}
    signature = hashlib.sha256(json.dumps(pin, sort_keys=True).encode()).hexdigest()
    work = root / "tmp/verified_runs" / manifest["release_id"] / f"{namespace}-{signature}"
    context = work / "VERIFIED_CONTEXT.json"
    if work.exists():
        prior = json.loads(context.read_text())
        if prior["pin"] != pin or prior["signature"] != signature:
            raise ValueError("Existing run context does not match source/protocol pin")
        for logical, expected in manifest["source_files"].items():
            if digest(confined(work, logical)) != expected:
                raise ValueError(f"Scratch source changed: {logical}")
        for logical, expected in package_files.items():
            if digest(confined(work / "r-library/prime", logical)) != expected:
                raise ValueError(f"Scratch installed package changed: {logical}")
        if namespace == "release-validation":
            verify_hydrated_inputs(work, manifest)
        return work, manifest, pin, signature
    work.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="hydrate-", dir=work.parent))
    try:
        for logical, record in manifest["files"].items():
            destination = confined(staging, logical)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with logical_input(confined(root, manifest["data_root"]), record) as source:
                with destination.open("wb") as target:
                    shutil.copyfileobj(source, target, 1024 * 1024)
        for logical in manifest["source_files"]:
            destination = confined(staging, logical)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(confined(root, logical), destination)
        shutil.copytree(package, staging / "r-library/prime")
        for logical, expected in manifest["source_files"].items():
            if digest(confined(staging, logical)) != expected:
                raise ValueError(f"Source changed during hydration: {logical}")
        for logical, expected in package_files.items():
            if digest(confined(staging / "r-library/prime", logical)) != expected:
                raise ValueError(f"Installed package changed during hydration: {logical}")
        for logical, record in manifest["files"].items():
            if digest(confined(staging, logical)) != record["sha256"]:
                raise ValueError(f"Input changed during hydration: {logical}")
        check_baseline_vectors(staging)
        (staging / "manuscript/figures").mkdir(parents=True, exist_ok=True)
        (staging / "manuscript/tables").mkdir(parents=True, exist_ok=True)
        write_json(staging / "VERIFIED_CONTEXT.json", {"pin": pin, "signature": signature,
                   "work_root": str(work), "release_status": manifest["status"]})
        (staging / "VERIFIED_CONTEXT.dcf").write_text(
            f"WorkRoot: {work}\nSignature: {signature}\nReleaseSHA256: {manifest_hash}\n")
        staging.rename(work)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return work, manifest, pin, signature


def invoke(root, work, manifest, script, signature):
    if script not in manifest["source_files"]:
        raise ValueError(f"Invoked script is not source-pinned: {script}")
    current, current_hash = verify_release(root)
    context = json.loads((work / "VERIFIED_CONTEXT.json").read_text())
    if current_hash != context["pin"]["release_sha256"]:
        raise ValueError("Release changed before script invocation")
    _, installed = package_snapshot(root, current, current_hash)
    if (installed != context["pin"]["package_files"] or
            digest(root / BOOTSTRAP_RECEIPT) != context["pin"]["bootstrap_receipt_sha256"]):
        raise ValueError("Local installation changed before script invocation")
    for logical, expected in current["source_files"].items():
        if digest(confined(work, logical)) != expected:
            raise ValueError(f"Scratch source changed: {logical}")
    for logical, expected in installed.items():
        if digest(confined(work / "r-library/prime", logical)) != expected:
            raise ValueError(f"Scratch installed package changed: {logical}")
    environment = r_environment(root, work / "r-library")
    environment["PRIME_VERIFIED_WORK_ROOT"] = str(work)
    environment["PRIME_VERIFIED_SIGNATURE"] = signature
    environment["PRIME_CANONICAL_ROOT"] = str(root)
    cache = work / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    environment["XDG_CACHE_HOME"] = str(cache)
    expression = (f'stopifnot(normalizePath(find.package("prime")) == '
                  f'normalizePath({json.dumps(str(work / "r-library/prime"))})); '
                  f'source({json.dumps(script)})')
    command = ["Rscript", "--vanilla", "-e", expression]
    record = {"command": command, "signature": signature, "started_unix": time.time(),
              "environment": {k: v for k, v in environment.items()
                              if k.startswith("PRIME_") or k in ("R_LIBS_USER", "XDG_CACHE_HOME")}}
    log = work / "verified_logs" / (Path(script).parent.name + "-" + Path(script).name + ".json")
    # Concurrent Make invocations must not update the same checkpoint tree.
    import fcntl
    with (work / "VERIFIED_RUN.lock").open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"Another script is using this verified namespace: {work}") from None
        write_json(log, record | {"status": "running"})
        context = json.loads((work / "VERIFIED_CONTEXT.json").read_text())
        immutable_inputs = context["pin"]["namespace"] == "release-validation"
        try:
            if immutable_inputs:
                verify_hydrated_inputs(work, manifest)
            result = subprocess.run(command, cwd=work, env=environment, check=False)
            if immutable_inputs:
                verify_hydrated_inputs(work, manifest)
        except (OSError, ValueError) as error:
            write_json(log, record | {"status": "failed", "error": str(error),
                                     "finished_unix": time.time()})
            raise
    write_json(log, record | {"status": "complete" if result.returncode == 0 else "failed",
                            "returncode": result.returncode, "finished_unix": time.time()})
    if result.returncode:
        raise RuntimeError(f"Verified script failed ({result.returncode}): {script}")


def artifacts(root):
    # Artifact generation always starts from the immutable verified CSVs, never
    # from partially rerun experiment outputs in the separate main namespace.
    work, manifest, pin, signature = prepare(root, f"artifacts-{time.time_ns()}")
    for script in ARTIFACT_SCRIPTS:
        invoke(root, work, manifest, script, signature)
    # Recheck immutable inputs/source before promoting generated display files.
    current, current_hash = verify_release(root)
    if current_hash != pin["release_sha256"]:
        raise ValueError("Current release changed while artifacts were building")
    package_snapshot(root, current, current_hash)
    outputs = {}
    for directory in ["manuscript/figures", "manuscript/tables", "results/summary"]:
        for source in sorted((work / directory).rglob("*")):
            if not source.is_file():
                continue
            logical = str(source.relative_to(work))
            if directory.startswith("manuscript/"):
                # Legacy temporary layout is retained; never restore document sources to the public checkout.
                if source.suffix.lower() in {".tex", ".bib", ".sty", ".cls", ".aux", ".bbl", ".bst"}:
                    continue
                logical = "build/historical_research_graphics/" + str(source.relative_to(work / "manuscript"))
            destination = confined(root, logical)
            destination.parent.mkdir(parents=True, exist_ok=True)
            backup = work / "prior_canonical_outputs" / logical
            if destination.exists() and not backup.exists():
                backup.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(destination, backup)
            temporary = destination.with_name(destination.name + f".verified-{os.getpid()}")
            shutil.copy2(source, temporary)
            temporary.replace(destination)
            outputs[logical] = digest(destination)
    write_json(root / "results/provenance/current_artifact_build.json",
               {"signature": signature, "work_root": str(work), "outputs": outputs,
                "release_sha256": pin["release_sha256"], "status": "built_for_local_review",
                "manuscript_source_export": "NOT_IN_PUBLIC_EXPORT"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["bootstrap", "run", "artifacts", "inspect"])
    parser.add_argument("--namespace", default="main")
    parser.add_argument("--script")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    if args.action == "bootstrap":
        bootstrap(root)
    elif args.action == "artifacts":
        artifacts(root)
    else:
        work, manifest, _, signature = prepare(root, args.namespace)
        if args.action == "run":
            if not args.script:
                parser.error("run requires --script")
            invoke(root, work, manifest, args.script, signature)
        print(json.dumps({"work_root": str(work), "signature": signature}))


if __name__ == "__main__":
    main()
