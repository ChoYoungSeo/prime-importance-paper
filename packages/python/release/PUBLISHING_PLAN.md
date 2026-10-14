# PyPI publishing

## Version 0.2.2 configuration

For `prime-importance` version `0.2.2`, configure the manual
`.github/workflows/publish.yml` workflow in `ChoYoungSeo/prime-importance-paper`
for release tag `v0.2.2` and these two files:

- `prime_importance-0.2.2-py3-none-any.whl`
- `prime_importance-0.2.2.tar.gz`

Publishing uses the matching PyPI Trusted Publisher and the GitHub `pypi`
environment restricted to the `main` branch. The workflow guards also require
the canonical repository and `main` branch. The `upload` input defaults to
`false`; main pushes, tags, and release events do not trigger an upload.

The verification job downloads the two release assets, checks their reviewed
SHA256 values, and runs strict Twine checks. A separate publishing job receives
those same files, rechecks their hashes, and uses OIDC. Only that job has
`id-token: write`. No personal API token is needed. Production duplicate
uploads fail rather than silently skipping files.

Keep the filenames, sizes, hashes, source commit, and actual publication result
in an external release receipt. The version, tag, artifact name, concurrency
group, and both hash checks must identify the same release. Run `upload=false`
for verification before invoking `upload=true`. After publication, verify the
public PyPI hashes and a fresh installation. Published filenames cannot be
replaced; later source or metadata changes belong to a new version.

## Historical published release: 0.2.1

`prime-importance` `0.2.1` was published on
[PyPI](https://pypi.org/project/prime-importance/0.2.1/). Its wheel and source
distribution match the approved original release files. Its installation,
dependency check, and quickstart were verified in a fresh Python 3.12 environment.

- Distribution: `prime-importance`; import: `prime`.
- Version: `0.2.1`; Python: >=3.10.
- Dependencies: NumPy >=1.24 and SciPy >=1.10; plotting: Matplotlib >=3.7.
- License: MIT.
- [GitHub release and original assets](https://github.com/ChoYoungSeo/prime-importance-paper/releases/tag/v0.2.1)
- Original PyPI publication provenance: [wheel](https://pypi.org/integrity/prime-importance/0.2.1/prime_importance-0.2.1-py3-none-any.whl/provenance) · [source distribution](https://pypi.org/integrity/prime-importance/0.2.1/prime_importance-0.2.1.tar.gz/provenance)

| Original file | SHA256 |
| --- | --- |
| `prime_importance-0.2.1.tar.gz` | `2b5ac6c8a12e1c8979bdb466e3b0193ffcfb5c3b63321161ea195cbcc327220d` |
| `prime_importance-0.2.1-py3-none-any.whl` | `c04a3696495d1b9d47d44ed9b045c7032f604f848671220434e6aff0ed7f3642` |

Version `0.2.1` must not be uploaded again. Its files and original publication
provenance remain unchanged.

## Official references

- [Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
- [Packaging and builds](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
- [PyPI immutable files](https://pypi.org/help/)
