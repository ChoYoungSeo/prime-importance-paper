import numpy as np
import pytest

from prime import prime_screen
from prime._core import prepare_data


def test_constant_predictor_is_rejected() -> None:
    X = np.column_stack((np.arange(20.0), np.ones(20)))
    y = np.arange(20.0)
    try:
        prime_screen(X, y)
    except ValueError as error:
        assert "constant" in str(error)
    else:
        raise AssertionError("Expected a constant predictor to fail.")


@pytest.mark.parametrize("names", [[], ["x"], ["x", "x"], ["x", 2], ["x", ["y"]]])
def test_invalid_names_are_rejected(names) -> None:
    rng = np.random.default_rng(1)
    X = rng.normal(size=(30, 2))
    y = rng.normal(size=30)
    with pytest.raises(ValueError, match="feature_names"):
        prime_screen(X, y, feature_names=names)


def test_numpy_feature_names_are_accepted() -> None:
    rng = np.random.default_rng(1)
    X = rng.normal(size=(30, 2))
    screen = prime_screen(X, X[:, 0], feature_names=np.array(["first", "second"]))
    assert tuple(screen.ridge_weights) == ("first", "second")


@pytest.mark.parametrize(
    "queries",
    [
        [0.5, 2],
        [False, 2],
        np.array([True, False]),
        [0, np.nan],
        [0, np.inf],
        [-1, 2],
        [0, 20],
        [1, 1.0],
        ["0", 2],
        [],
        [[0, 1]],
    ],
)
def test_invalid_query_indices_are_rejected_before_conversion(queries) -> None:
    X = np.arange(20.0).reshape(-1, 1)
    with pytest.raises(ValueError, match="query_indices"):
        prime_screen(X, X[:, 0], query_indices=queries)


def test_whole_valued_float_query_indices_preserve_rows_and_order() -> None:
    X = np.arange(20.0).reshape(-1, 1)
    data = prepare_data(
        X,
        X[:, 0],
        feature_names=None,
        rescale=True,
        query_indices=[4.0, 0.0, 19.0],
        lambda_=1e-6,
    )
    np.testing.assert_array_equal(data.query_indices, [4, 0, 19])
    assert data.query_indices.dtype == np.int64
