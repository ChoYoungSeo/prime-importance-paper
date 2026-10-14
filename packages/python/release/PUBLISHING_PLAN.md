# PyPI publishing

`prime-importance` 0.2.1 is published on [PyPI](https://pypi.org/project/prime-importance/0.2.1/).
The wheel and source distribution match the approved release files. The public
installation, dependency check, and quickstart were verified in a fresh Python
3.12 environment.

## Published release

- Distribution: `prime-importance`; import: `prime`.
- Version: `0.2.1`; Python: >=3.10.
- Dependencies: NumPy >=1.24 and SciPy >=1.10; optional plotting: Matplotlib >=3.7.
- License: MIT.
- [GitHub release and approved assets](https://github.com/ChoYoungSeo/prime-importance-paper/releases/tag/v0.2.1)
- Original PyPI publication provenance: [wheel](https://pypi.org/integrity/prime-importance/0.2.1/prime_importance-0.2.1-py3-none-any.whl/provenance) · [source distribution](https://pypi.org/integrity/prime-importance/0.2.1/prime_importance-0.2.1.tar.gz/provenance)
- [Paper](https://arxiv.org/abs/2610.17053) · [DOI](https://doi.org/10.48550/arXiv.2610.17053)

| File | Approved SHA256 |
| --- | --- |
| `prime_importance-0.2.1.tar.gz` | `2b5ac6c8a12e1c8979bdb466e3b0193ffcfb5c3b63321161ea195cbcc327220d` |
| `prime_importance-0.2.1-py3-none-any.whl` | `c04a3696495d1b9d47d44ed9b045c7032f604f848671220434e6aff0ed7f3642` |

## Publishing workflow configuration

`.github/workflows/publish.yml` is configured for the canonical repository
`ChoYoungSeo/prime-importance-paper` and its `main` branch. Publishing uses the
matching PyPI Trusted Publisher and a GitHub `pypi` environment restricted
to `main`.

The workflow runs only through `workflow_dispatch`; main pushes, tags, and
release events do not trigger an upload. The `upload` input defaults to `false`.

The verification job downloads the two approved GitHub release assets, checks
their fixed SHA256 values, and runs strict Twine checks. A separate publishing
job receives those same files, rechecks their hashes, and uses OIDC. Only that
job has `id-token: write`. No personal API token is needed. Production duplicate
uploads fail rather than silently skipping files.

The workflow is fixed to 0.2.1. After the canonical repository connection and
approved release assets are in place, `upload=false` verify
the release files without uploading them. Version 0.2.1 is already published
and should not be uploaded again. Its published files and original PyPI
publication provenance remain unchanged.

## Future versions

1. Prepare and review the new source, version metadata, README, and citation.
2. Build the source distribution, then build the wheel from that archive.
   Check both archives, strict Twine metadata, installed behavior, and the
   applicable tests for the exact new files.
3. Freeze filenames, sizes, SHA256 values, and the source commit in a release
   receipt. Attach those exact files to the corresponding GitHub release.
4. Update the workflow's version, release tag, filenames, artifact name,
   concurrency group, and both hash checks to the newly reviewed release.
5. Run `upload=false`. After verification and release review,
   explicitly run `upload=true`.
6. Verify the public PyPI file hashes and a fresh installation. Record the
   actual publishing result. Published filenames cannot be replaced; later
   source or metadata changes belong to a new version.

## Official references

- [Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
- [Packaging and builds](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
- [PyPI immutable files](https://pypi.org/help/)
