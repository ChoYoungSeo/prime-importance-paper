import numpy as np
import pytest

from prime import PrimeAllocation, ProjectionDiagnostics


@pytest.fixture
def allocation() -> PrimeAllocation:
    return PrimeAllocation(
        effects=np.array([0.3, 0.7]),
        projected_values=np.array([0.0, 0.3, 0.7, 1.0]),
        projection=ProjectionDiagnostics(0, True, 0.0),
    )


def test_allocation_array_protocol_shares_or_copies_as_requested(allocation) -> None:
    assert allocation.__array__(copy=False) is allocation.effects
    assert allocation.__array__(dtype=np.float64, copy=False) is allocation.effects
    assert np.shares_memory(np.asarray(allocation), allocation.effects)

    copied = allocation.__array__(copy=True)
    np.testing.assert_array_equal(copied, allocation.effects)
    assert not np.shares_memory(copied, allocation.effects)
    converted = allocation.__array__(dtype=np.float32)
    assert converted.dtype == np.float32
    np.testing.assert_allclose(converted, allocation.effects)


def test_allocation_rejects_dtype_conversion_without_copy(allocation) -> None:
    with pytest.raises(ValueError, match="copy"):
        allocation.__array__(dtype=np.float32, copy=False)


@pytest.mark.skipif(
    np.lib.NumpyVersion(np.__version__) < "2.0.0",
    reason="Strict copy=False array conversion requires NumPy 2.",
)
def test_numpy_two_honors_allocation_copy_protocol(allocation) -> None:
    converted = np.array(allocation, copy=False)
    assert np.shares_memory(converted, allocation.effects)
    with pytest.raises(ValueError, match="copy"):
        np.array(allocation, dtype=np.float32, copy=False)
