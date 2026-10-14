"""Bounded integrity contracts; no simulations or compiler required."""
import copy
import builtins
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

spec = importlib.util.spec_from_file_location('workflow', Path(__file__).with_name('verified_workflow.py'))
workflow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workflow)


class WorkflowIntegrity(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='prime workflow ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        def put(name, data):
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            return workflow.digest(path)
        self.put = put
        sources = {p: put(p, value) for p, value in {
            'software/prime/DESCRIPTION': b'Package: prime\nVersion: 0.5.2.9101\n',
            'software/prime/src/Makevars.win': b'PKG_LIBS = $(BLAS_LIBS)\n',
            'software/prime/src/native.h': b'// pinned header\n',
            'software/prime/README.md': b'Frozen package documentation\n',
            'experiments/99_validate/verified_workflow.py': b'# fixture\n',
        }.items()}
        data_root = 'results/revisions/test'
        files = {}
        for name in workflow.BASELINES:
            path = f'results/provenance/paired_baselines/{name}.csv.gz'
            files[path] = {'path': path, 'sha256': put(f'{data_root}/{path}', b'fixture')}
        self.payload = b'id,value\n1,0.125\n2,0.875\n'
        self.transport = 'results/raw/example.csv.gz'
        compressed = gzip.compress(self.payload, mtime=0)
        files['results/raw/example.csv'] = {
            'path': 'results/raw/example.csv',
            'sha256': hashlib.sha256(self.payload).hexdigest(),
            'transport': {'encoding': 'gzip', 'path': self.transport,
                          'sha256': put(f'{data_root}/{self.transport}', compressed),
                          'uncompressed_bytes': len(self.payload)}}
        lineage = 'results/provenance/lineage.json'
        self.manifest = {'schema_version': 2, 'status': 'validated_public_release',
                         'release_id': 'test', 'data_root': data_root,
                         'source_files': sources, 'files': files,
                         'package_version': '0.5.2.9101',
                         'historical_lineage': {'manifest': lineage, 'sha256': put(lineage, b'{}\n')}}
        self.save_manifest()
        self.runtime = {'R_version': 'fixture', 'platform': 'fixture', 'dependencies': {}}
        self.put('r-library/prime/DESCRIPTION', b'Package: prime\nVersion: 0.5.2.9101\n')
        self.put('r-library/prime/libs/prime.so', b'fresh compiled fixture')
        _, installed = workflow.installed_snapshot(self.root/'r-library/prime', '0.5.2.9101')
        self.receipt = {'schema_version': 1, 'release_sha256': self.manifest_hash,
                        'source_files': sources,
                        'input_sha256': {p:r['sha256'] for p,r in files.items()},
                        'package_path': 'r-library/prime', 'installed_files': installed,
                        'smoke_validation': 'passed', 'runtime': self.runtime}
        self.save_receipt()
        self.runtime_patch = mock.patch.object(workflow, 'runtime_snapshot', return_value=self.runtime)
        self.runtime_patch.start()
        self.addCleanup(self.runtime_patch.stop)

    def save_manifest(self):
        workflow.write_json(self.root/'results/current_release.json', self.manifest)
        self.manifest_hash = workflow.digest(self.root/'results/current_release.json')

    def save_receipt(self):
        workflow.write_json(self.root/workflow.BOOTSTRAP_RECEIPT, self.receipt)

    def snapshot(self):
        return workflow.package_snapshot(self.root, self.manifest, self.manifest_hash)

    def workflow_revision(self, revised):
        path = "results/workflow_revisions/test.json"
        workflow.write_json(self.root/path, revised)
        workflow.write_json(self.root/"results/current_workflow.json", {
            "schema_version": 1, "manifest": path,
            "sha256": workflow.digest(self.root/path),
            "historical_release_sha256": workflow.digest(self.root/"results/current_release.json")})

    def add_maintained_packages(self):
        self.maintained = {'schema_version': 1, 'packages': {}}
        fixtures = {
            'r': ('0.5.5', {
                'DESCRIPTION': b'Package: prime\nVersion: 0.5.5\n',
                'R/prime.R': b'prime <- function(x) x\n'}),
            'python': ('0.2.1', {
                'pyproject.toml': b'[project]\nname = "prime-importance"\nversion = "0.2.1"\n',
                'src/prime_importance/__init__.py': b'__version__ = "0.2.1"\n'})}
        for language, (version, files) in fixtures.items():
            package_root = f'packages/{language}'
            self.maintained['packages'][language] = {
                'root': package_root, 'version': version,
                'source_files': {f'{package_root}/{name}': self.put(f'{package_root}/{name}', data)
                                 for name, data in files.items()}}
        self.workflow_revision(copy.deepcopy(self.manifest))
        self.maintained_path = 'results/maintained_packages/test.json'
        self.save_maintained_manifest()

    def save_maintained_manifest(self):
        workflow.write_json(self.root/self.maintained_path, self.maintained)
        self.pin_maintained_manifest()

    def pin_maintained_manifest(self):
        pointer_path = self.root/'results/current_workflow.json'
        pointer = json.loads(pointer_path.read_text())
        pointer['maintained_packages'] = {
            'manifest': self.maintained_path,
            'sha256': workflow.digest(self.root/self.maintained_path)}
        workflow.write_json(pointer_path, pointer)

    def test_optional_maintained_reference_absent_preserves_original_guard(self):
        self.add_maintained_packages()
        pointer_path = self.root/'results/current_workflow.json'
        pointer = json.loads(pointer_path.read_text())
        del pointer['maintained_packages']
        workflow.write_json(pointer_path, pointer)
        self.put('packages/python/src/prime_importance/__init__.py', b'outside optional pins\n')
        actual, manifest_hash = workflow.verify_release(self.root)
        self.assertEqual(actual, self.manifest)
        self.assertEqual(manifest_hash, self.manifest_hash)
        self.snapshot()

    def test_valid_maintained_packages_preserve_frozen_manifest_and_receipt(self):
        self.add_maintained_packages()
        before = (self.root/'results/current_release.json').read_bytes()
        actual, manifest_hash = workflow.verify_release(self.root)
        self.assertEqual(actual, self.manifest)
        self.assertEqual(manifest_hash, self.manifest_hash)
        self.assertEqual(before, (self.root/'results/current_release.json').read_bytes())
        self.assertEqual(actual['package_version'], '0.5.2.9101')
        self.snapshot()

    def test_maintained_reference_null_is_not_treated_as_absent(self):
        self.add_maintained_packages()
        pointer_path = self.root/'results/current_workflow.json'
        pointer = json.loads(pointer_path.read_text())
        pointer['maintained_packages'] = None
        workflow.write_json(pointer_path, pointer)
        with self.assertRaisesRegex(ValueError, 'manifest reference'):
            workflow.verify_release(self.root)

    def test_maintained_code_tampering_is_rejected(self):
        self.add_maintained_packages()
        self.put('packages/r/R/prime.R', b'prime <- function(x) NULL\n')
        with self.assertRaisesRegex(ValueError, 'Maintained package source hash mismatch'):
            workflow.verify_release(self.root)

    def test_unpinned_maintained_source_is_rejected(self):
        self.add_maintained_packages()
        self.put('packages/python/src/prime_importance/new.py', b'new_algorithm = True\n')
        with self.assertRaisesRegex(ValueError, 'source inventory differs'):
            workflow.verify_release(self.root)

    def test_missing_maintained_source_is_rejected(self):
        self.add_maintained_packages()
        (self.root/'packages/r/R/prime.R').unlink()
        with self.assertRaisesRegex(ValueError, 'source inventory differs'):
            workflow.verify_release(self.root)

    def test_maintained_manifest_hash_is_verified(self):
        self.add_maintained_packages()
        (self.root/self.maintained_path).write_text('{}\n')
        with self.assertRaisesRegex(ValueError, 'Maintained package manifest hash mismatch'):
            workflow.verify_release(self.root)

    def test_maintained_manifest_path_traversal_is_rejected(self):
        self.add_maintained_packages()
        pointer_path = self.root/'results/current_workflow.json'
        pointer = json.loads(pointer_path.read_text())
        pointer['maintained_packages']['manifest'] = '../outside.json'
        workflow.write_json(pointer_path, pointer)
        with self.assertRaisesRegex(ValueError, 'safe relative path'):
            workflow.verify_release(self.root)

    def test_maintained_source_paths_are_canonical_and_package_confined(self):
        self.add_maintained_packages()
        original = copy.deepcopy(self.maintained)
        for invalid in ('../outside.py', 'packages/python/../outside.py',
                        './packages/python/src/new.py', 'packages//python/src/new.py',
                        'packages/python/src/new.py/', 'packages/r/R/prime.R',
                        'software/prime/src/native.h', '/outside.py',
                        'packages\\python\\outside.py'):
            with self.subTest(path=invalid):
                self.maintained = copy.deepcopy(original)
                self.maintained['packages']['python']['source_files'][invalid] = '0'*64
                self.save_maintained_manifest()
                with self.assertRaises(ValueError):
                    workflow.verify_release(self.root)

    def test_maintained_schema_roots_and_hashes_are_strict(self):
        self.add_maintained_packages()
        original = copy.deepcopy(self.maintained)
        invalid_manifests = []
        for schema in (True, 2, '1'):
            value = copy.deepcopy(original)
            value['schema_version'] = schema
            invalid_manifests.append(value)
        value = copy.deepcopy(original)
        value['extra'] = 'unreviewed field'
        invalid_manifests.append(value)
        value = copy.deepcopy(original)
        value['packages']['r']['root'] = 'software/prime'
        invalid_manifests.append(value)
        value = copy.deepcopy(original)
        value['packages']['r']['source_files']['packages/r/DESCRIPTION'] = 'A'*64
        invalid_manifests.append(value)
        value = copy.deepcopy(original)
        del value['packages']['python']
        invalid_manifests.append(value)
        for invalid in invalid_manifests:
            with self.subTest(manifest=invalid):
                self.maintained = invalid
                self.save_maintained_manifest()
                with self.assertRaises(ValueError):
                    workflow.verify_release(self.root)

    def test_maintained_manifest_duplicate_keys_are_rejected(self):
        self.add_maintained_packages()
        text = json.dumps(self.maintained).replace('"schema_version": 1',
                                                  '"schema_version": 1, "schema_version": 1')
        (self.root/self.maintained_path).write_text(text)
        self.pin_maintained_manifest()
        with self.assertRaisesRegex(ValueError, 'Duplicate maintained package manifest key'):
            workflow.verify_release(self.root)

    def test_rehashed_maintained_metadata_must_match_name_and_version(self):
        self.add_maintained_packages()
        original = copy.deepcopy(self.maintained)
        r_description = (self.root/'packages/r/DESCRIPTION').read_bytes()
        py_project = (self.root/'packages/python/pyproject.toml').read_bytes()
        for language, logical, data in (
                ('r', 'packages/r/DESCRIPTION', b'Package: prime\nVersion: 9.9.9\n'),
                ('r', 'packages/r/DESCRIPTION', b'Package: other\nVersion: 0.5.5\n'),
                ('r', 'packages/r/DESCRIPTION', b'Package: prime\nVersion: 0.5.5\nVersion: 0.5.5\n'),
                ('r', 'packages/r/DESCRIPTION', b'Package: prime\nVersion: 0.5.5\n continuation\n'),
                ('python', 'packages/python/pyproject.toml', py_project.replace(b'0.2.1', b'9.9.9')),
                ('python', 'packages/python/pyproject.toml', py_project.replace(b'prime-importance', b'other')),
                ('python', 'packages/python/pyproject.toml', py_project + b'version = "0.2.1"\n'),
                ('python', 'packages/python/pyproject.toml', py_project + b'"version" = "0.2.1"\n'),
                ('python', 'packages/python/pyproject.toml', py_project + b'dynamic = ["version"]\n'),
                ('python', 'packages/python/pyproject.toml', py_project + b'[project]\nversion = "0.2.1"\n')):
            with self.subTest(language=language, metadata=data):
                self.put('packages/r/DESCRIPTION', r_description)
                self.put('packages/python/pyproject.toml', py_project)
                self.maintained = copy.deepcopy(original)
                self.maintained['packages'][language]['source_files'][logical] = self.put(logical, data)
                self.save_maintained_manifest()
                with self.assertRaisesRegex(ValueError, '(name/version|metadata|project)'):
                    workflow.verify_release(self.root)

    def test_maintained_build_and_cache_artifacts_are_ignored(self):
        self.add_maintained_packages()
        ignored_directories = ('__pycache__', '.pytest_cache', '.ruff_cache', 'build', 'dist',
                               'prime_importance.egg-info', 'library', 'tmp', 'check', 'prime.Rcheck')
        for package_root in ('packages/r', 'packages/python'):
            for name in ignored_directories:
                self.put(f'{package_root}/{name}/generated.py', b'transient generated artifact\n')
            for name in ('stale.o', 'stale.so', 'stale.dll', 'stale.a', '.DS_Store',
                         'legacy.egg-info', 'temporary.Rcheck'):
                self.put(f'{package_root}/src/{name}', b'transient compiled artifact\n')
        workflow.verify_release(self.root)
        for package in self.maintained['packages'].values():
            self.assertEqual(workflow.maintained_package_source_inventory(self.root, package['root']),
                             set(package['source_files']))

    def test_python_310_static_metadata_fallback_accepts_current_format(self):
        self.add_maintained_packages()
        original_import = builtins.__import__
        def without_tomllib(name, *args, **kwargs):
            if name == 'tomllib':
                raise ModuleNotFoundError('Python 3.10 fixture')
            return original_import(name, *args, **kwargs)
        with mock.patch('builtins.__import__', side_effect=without_tomllib):
            workflow.verify_release(self.root)

    def test_python_310_metadata_fallback_rejects_fake_table_in_string(self):
        self.add_maintained_packages()
        metadata = b'''content = """\n[project]\nname = "prime-importance"\nversion = "0.2.1"\n[foo]\n"""\n["project"]\nname = "other"\nversion = "9.9.9"\n'''
        logical = 'packages/python/pyproject.toml'
        self.maintained['packages']['python']['source_files'][logical] = self.put(logical, metadata)
        self.save_maintained_manifest()
        original_import = builtins.__import__
        def without_tomllib(name, *args, **kwargs):
            if name == 'tomllib':
                raise ModuleNotFoundError('Python 3.10 fixture')
            return original_import(name, *args, **kwargs)
        with mock.patch('builtins.__import__', side_effect=without_tomllib):
            with self.assertRaisesRegex(ValueError, 'metadata'):
                workflow.verify_release(self.root)

    def test_maintained_source_symlinks_are_rejected(self):
        self.add_maintained_packages()
        destination = self.root/'outside-source.py'
        destination.write_bytes(b'external source\n')
        for logical, target in (('packages/python/src/linked.py', destination),
                                ('packages/python/src/broken.py', self.root/'missing.py'),
                                ('packages/python/build', self.root)):
            with self.subTest(path=logical):
                link = self.root/logical
                link.parent.mkdir(parents=True, exist_ok=True)
                link.symlink_to(target, target_is_directory=target.is_dir())
                try:
                    with self.assertRaisesRegex(ValueError, 'source symlink'):
                        workflow.verify_release(self.root)
                finally:
                    link.unlink()

    def test_symlinked_maintained_package_ancestor_is_rejected(self):
        self.add_maintained_packages()
        (self.root/'packages').rename(self.root/'actual-packages')
        (self.root/'packages').symlink_to(self.root/'actual-packages', target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'source symlink'):
            workflow.verify_release(self.root)

    def test_symlinked_maintained_package_root_is_rejected(self):
        self.add_maintained_packages()
        (self.root/'packages/r').rename(self.root/'packages/actual-r')
        (self.root/'packages/r').symlink_to(self.root/'packages/actual-r', target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'source symlink'):
            workflow.verify_release(self.root)

    def test_maintained_pins_do_not_relax_frozen_source_guard(self):
        self.add_maintained_packages()
        self.put('software/prime/src/native.h', b'changed frozen implementation\n')
        with self.assertRaisesRegex(ValueError, 'Canonical source hash mismatch'):
            workflow.verify_release(self.root)

    def test_maintained_pins_do_not_allow_historical_inventory_extension(self):
        self.add_maintained_packages()
        revised = copy.deepcopy(self.manifest)
        revised['source_files'].update(self.maintained['packages']['r']['source_files'])
        self.workflow_revision(revised)
        self.pin_maintained_manifest()
        with self.assertRaisesRegex(ValueError, 'Workflow revision changed source inventory'):
            workflow.verify_release(self.root)

    def test_workflow_revision_preserves_historical_manifest(self):
        before = (self.root/"results/current_release.json").read_bytes()
        revised = copy.deepcopy(self.manifest)
        name = "experiments/99_validate/verified_workflow.py"
        revised["source_files"][name] = self.put(name, b"# audited orchestration repair\n")
        self.workflow_revision(revised)
        actual, _ = workflow.verify_release(self.root)
        self.assertEqual(actual, revised)
        self.assertEqual(before, (self.root/"results/current_release.json").read_bytes())

    def test_workflow_revision_rejects_scientific_source_repin(self):
        revised = copy.deepcopy(self.manifest)
        name = "software/prime/src/native.h"
        revised["source_files"][name] = self.put(name, b"changed method\n")
        self.workflow_revision(revised)
        with self.assertRaisesRegex(ValueError, "scientific source pins"):
            workflow.verify_release(self.root)

    def test_workflow_revision_can_correct_frozen_package_documentation(self):
        before = (self.root / 'results/current_release.json').read_bytes()
        revised = copy.deepcopy(self.manifest)
        name = 'software/prime/README.md'
        revised['source_files'][name] = self.put(name, b'Corrected relative links\n')
        self.workflow_revision(revised)
        actual, _ = workflow.verify_release(self.root)
        self.assertEqual(actual, revised)
        self.assertEqual(before, (self.root / 'results/current_release.json').read_bytes())
        self.assertEqual(actual['files'], self.manifest['files'])
        self.assertEqual(actual['source_files']['software/prime/src/native.h'],
                         self.manifest['source_files']['software/prime/src/native.h'])

    def test_workflow_revision_rejects_scientific_input_repin(self):
        revised = copy.deepcopy(self.manifest)
        revised["files"]["results/raw/example.csv"]["sha256"] = "0"*64
        self.workflow_revision(revised)
        with self.assertRaisesRegex(ValueError, "scientific release metadata"):
            workflow.verify_release(self.root)

    def test_workflow_revision_rejects_unpinned_manifest(self):
        self.workflow_revision(self.manifest)
        (self.root/"results/workflow_revisions/test.json").write_text("{}\n")
        with self.assertRaisesRegex(ValueError, "manifest hash mismatch"):
            workflow.verify_release(self.root)

    def test_ambient_R_LIBS_cannot_preempt_attested_library(self):
        with mock.patch.dict(workflow.os.environ, {"R_LIBS": "/different/prime",
                                                    "R_LIBS_USER": "/dependencies"}):
            env = workflow.r_environment(self.root, self.root/'scratch/r-library')
        self.assertNotIn('R_LIBS', env)
        self.assertEqual(env['R_LIBS_USER'].split(workflow.os.pathsep)[:2],
                         [str(self.root/'scratch/r-library'), str(self.root/'r-library')])

    def test_valid_transport_preserves_original_bytes(self):
        workflow.verify_release(self.root)
        record = self.manifest['files']['results/raw/example.csv']
        with workflow.logical_input(self.root/self.manifest['data_root'], record) as handle:
            self.assertEqual(handle.read(), self.payload)
        self.snapshot()

    def test_transport_tampering_is_rejected(self):
        self.put(self.manifest['data_root']+'/'+self.transport, b'corrupt')
        with self.assertRaisesRegex(ValueError, 'transport hash'):
            workflow.verify_release(self.root)

    def test_recompressed_different_logical_bytes_are_rejected(self):
        record = self.manifest['files']['results/raw/example.csv']
        record['transport']['sha256'] = self.put(self.manifest['data_root']+'/'+self.transport,
                                                gzip.compress(b'changed data', mtime=0))
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, 'input hash'):
            workflow.verify_release(self.root)

    def test_modified_native_header_is_rejected(self):
        self.put('software/prime/src/native.h', b'changed')
        with self.assertRaisesRegex(ValueError, 'source hash'):
            workflow.verify_release(self.root)

    def test_unpinned_build_file_is_rejected(self):
        self.put('software/prime/src/Makevars', b'new flags')
        with self.assertRaisesRegex(ValueError, 'build-source inventory'):
            workflow.verify_release(self.root)

    def test_stale_compiled_objects_are_not_build_inputs(self):
        self.put('software/prime/src/stale.o', b'stale object')
        workflow.verify_release(self.root)
        self.assertNotIn('software/prime/src/stale.o', workflow.package_source_inventory(self.root))

    def test_missing_receipt_is_rejected(self):
        (self.root/workflow.BOOTSTRAP_RECEIPT).unlink()
        with self.assertRaisesRegex(ValueError, 'bootstrap first'):
            self.snapshot()

    def test_changed_source_manifest_invalidates_receipt(self):
        self.manifest['source_files'] = copy.deepcopy(self.manifest['source_files'])
        self.manifest['source_files']['software/prime/src/native.h'] = 'different'
        with self.assertRaisesRegex(ValueError, 'receipt is stale'):
            self.snapshot()

    def test_changed_release_receipt_is_rejected(self):
        self.receipt['release_sha256'] = 'old release'
        self.save_receipt()
        with self.assertRaisesRegex(ValueError, 'receipt is stale'):
            self.snapshot()

    def test_changed_installed_binary_is_rejected(self):
        self.put('r-library/prime/libs/prime.so', b'different installed code')
        with self.assertRaisesRegex(ValueError, 'Installed PRIME files'):
            self.snapshot()

    def test_changed_runtime_is_rejected(self):
        self.receipt['runtime'] = {'R_version': 'old runtime'}
        self.save_receipt()
        with self.assertRaisesRegex(ValueError, 'versions changed'):
            self.snapshot()

    def test_escaping_transport_is_rejected(self):
        self.manifest['files']['results/raw/example.csv']['transport']['path'] = '../outside.gz'
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, 'safe relative path'):
            workflow.verify_release(self.root)


if __name__ == '__main__':
    unittest.main()
