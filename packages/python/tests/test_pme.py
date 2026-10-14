import numpy as np
import pytest

from prime import prime_allocate


def test_two_variable_allocation_matches_reference() -> None:
    result = prime_allocate([0.0, 0.25, 0.60, 1.0])
    np.testing.assert_allclose(result.effects, [0.2941176471, 0.7058823529])
    np.testing.assert_allclose(result.projected_values, [0.0, 0.25, 0.60, 1.0])
    assert result.projection.converged


def test_three_variable_allocation_matches_r_reference() -> None:
    result = prime_allocate([0.0, 0.4, 0.2, 0.7, 0.3, 0.8, 0.6, 1.0])
    np.testing.assert_allclose(
        result.effects,
        [0.43545879, 0.22861586, 0.33592535],
        rtol=1e-7,
    )
    assert np.isclose(result.effects.sum(), 1.0)


def test_projection_repairs_nonmonotone_values() -> None:
    result = prime_allocate([0.0, 0.8, 0.2, 0.6])
    projected = result.projected_values
    assert projected[0] == 0.0
    assert projected[-1] == 0.6
    assert projected[1] <= projected[-1]
    assert projected[2] <= projected[-1]
    assert np.isclose(result.effects.sum(), 0.6)


def test_invalid_subset_values_length_is_rejected() -> None:
    try:
        prime_allocate([0.0, 0.3, 1.0])
    except ValueError as error:
        assert "2**d" in str(error)
    else:
        raise AssertionError("Expected an invalid subset-values length to fail.")


def test_positive_effect_below_subset_tolerance_preserves_efficiency() -> None:
    # The tolerance identifies null game values, not small positive allocations.
    for tolerance in (float(np.sqrt(np.finfo(float).eps)), 0.1):
        a = tolerance + tolerance**2 / 2
        game = np.array([0.0, a, 1.0, 1.0])
        allocation = prime_allocate(game, zero_tolerance=tolerance)
        expected = np.array([a, 1.0]) / (1.0 + a)
        assert 0 < allocation.effects[0] < tolerance
        np.testing.assert_allclose(allocation.effects, expected, atol=1e-14, rtol=0)
        np.testing.assert_allclose(allocation.effects.sum(), 1.0, atol=1e-14, rtol=0)
        np.testing.assert_allclose(allocation.projected_values, game, atol=1e-14, rtol=0)
        null = prime_allocate([0.0, tolerance / 2, 1.0, 1.0], zero_tolerance=tolerance)
        np.testing.assert_array_equal(null.effects, [0.0, 1.0])


@pytest.mark.parametrize("scale", [1e-200, 1e-12, 1e-9, 1.0, 1e100, 1e200])
@pytest.mark.parametrize("tolerance", [0.0, float(np.sqrt(np.finfo(float).eps))])
def test_allocation_preserves_game_units(scale: float, tolerance: float) -> None:
    game = np.array([0.0, 0.2, 0.8, 1.0])
    result = prime_allocate(scale * game, zero_tolerance=tolerance)
    assert np.all(np.isfinite(result.effects))
    np.testing.assert_allclose(result.effects / scale, [0.2, 0.8], rtol=1e-13)
    np.testing.assert_allclose(result.projected_values / scale, game, rtol=1e-13)
    assert result.projected_values[0] == 0
    assert result.projected_values[-1] == scale


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_nonmonotone_projection_is_scale_equivariant(scale: float) -> None:
    game = np.array([0.0, 0.8, 0.2, 0.6])
    result = prime_allocate(scale * game)
    np.testing.assert_allclose(result.projected_values / scale, [0.0, 0.6, 0.2, 0.6])
    np.testing.assert_allclose(result.effects / scale, [0.45, 0.15])


def test_zero_tolerance_preserves_tiny_positive_game_values() -> None:
    game = np.array([0.0, 1e-300, 1.0, 1.0])
    result = prime_allocate(game, zero_tolerance=0)
    np.testing.assert_array_equal(result.projected_values, game)
    np.testing.assert_allclose(result.effects, [1e-300, 1.0], rtol=1e-12, atol=0)


def test_log_potentials_handle_products_below_float_range() -> None:
    game = np.full(8, 1e-300)
    game[0], game[-1] = 0.0, 1.0
    result = prime_allocate(game, zero_tolerance=0)
    np.testing.assert_allclose(result.effects, np.full(3, 1 / 3), rtol=1e-12)
    np.testing.assert_array_equal(result.projected_values, game)


@pytest.mark.parametrize("tolerance", [1.0, 2.0])
def test_null_threshold_preserves_full_endpoint(tolerance: float) -> None:
    result = prime_allocate([0.0, 0.2, 0.8, 1.0], zero_tolerance=tolerance)
    np.testing.assert_array_equal(result.projected_values, [0.0, 0.0, 0.0, 1.0])
    np.testing.assert_allclose(result.effects, [0.5, 0.5])


def test_zero_full_value_projects_to_the_unique_zero_game() -> None:
    result = prime_allocate([0.0, -3.0, 4.0, 0.0], zero_tolerance=0)
    np.testing.assert_array_equal(result.projected_values, np.zeros(4))
    np.testing.assert_array_equal(result.effects, np.zeros(2))
    assert result.projection.converged


@pytest.mark.parametrize("game", [[0, 1e300, 0, 1e-300], [0, 1e-300, 0, 1e300]])
def test_unrepresentable_normalized_game_is_rejected(game: list[float]) -> None:
    with pytest.raises(ValueError, match="dynamic range"):
        prime_allocate(game, zero_tolerance=0)


@pytest.mark.parametrize("tolerance", [True, [0.1], -0.1, np.inf, np.nan])
def test_invalid_zero_tolerance_is_rejected(tolerance: object) -> None:
    with pytest.raises(ValueError, match="zero_tolerance"):
        prime_allocate([0.0, 0.2, 0.8, 1.0], zero_tolerance=tolerance)
