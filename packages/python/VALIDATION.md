# Validation

## Metadata update: 0.2.2 — 2026-10-14

Version `0.2.2` updates project links and package documentation to the canonical
PRIME repository. The seven scientific modules are unchanged from `0.2.1`
except for the `__version__` constant. Algorithms, public API, dependencies,
tests, examples, figures, and benchmark data are preserved.

All 91 existing tests passed locally on macOS arm64 with Python 3.13.3,
NumPy 2.4.0, SciPy 1.16.3, and Matplotlib 3.10.8. These local checks did not
repeat the operating-system/Python-version CI matrix.

The source distribution was built first, and the wheel was built from that
archive. Both archives passed strict Twine metadata checks. The installed
wheel passed the basic example with plotting imports blocked, and passed the
plot example and standalone quickstart with Matplotlib. Both examples selected
`X1`, `X2`, and `X5`.

The published `0.2.1` files and their original provenance are preserved.
