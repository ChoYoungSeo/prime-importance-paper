"""Bounded integrity contracts; no simulations or compiler required."""
import copy
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
