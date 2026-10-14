import numpy as np
import pytest

import prime as prime_package
from prime import prime, prime_screen
from prime._core import forward_score, local_risk_details, neighborhood_size, prepare_data


def make_example() -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(5)
    X = rng.uniform(-1, 1, size=(120, 8))
    y = X[:, 0] - X[:, 1] + X[:, 4] ** 2 + rng.normal(0, 0.08, size=120)
    return X, y


def test_short_api_aliases_are_stable() -> None:
    assert prime_package.fit is prime_package.prime
    assert prime_package.screen is prime_package.prime_screen
    assert prime_package.allocate is prime_package.prime_allocate


def test_screen_recovers_the_active_coordinates() -> None:
    X, y = make_example()
    screen = prime_screen(X, y)
    assert screen.selected_names == ("X1", "X2", "X5")
    assert screen.forward_trace
    assert screen.backward_trace
    assert "alpha" not in screen.settings
    assert (
        screen.settings["ridge_chosen_size"]
        == min(screen.settings["ridge_path"], key=lambda row: row["bic"])["size"]
    )
    for row in screen.forward_trace:
        expected = row["query_scale"] * np.sqrt(2 * np.log(max(2, row["dimension"])))
        assert np.isclose(row["threshold"], expected)
        assert np.isclose(
            row["calibrated_gain"],
            max(row["raw_gain"] - row["threshold"], 0.0),
        )
    for row in screen.backward_trace:
        expected = row["query_scale"] * np.sqrt(2 * np.log(max(2, row["dimension"])))
        assert np.isclose(row["threshold"], expected)
        assert np.isclose(
            row["calibrated_loss"],
            max(row["raw_loss"] - row["threshold"], 0.0),
        )
        assert row["removed"] == (row["calibrated_loss"] == 0.0)


def test_forward_ranking_uses_unclipped_risk_scores() -> None:
    rng = np.random.default_rng(40)
    X = rng.normal(size=(15, 8))
    y = rng.normal(size=15)
    data = prepare_data(
        X,
        y,
        feature_names=None,
        rescale=True,
        query_indices=None,
        lambda_=1e-6,
    )
    selected = [0]
    candidates = list(range(1, data.p))
    k = neighborhood_size(data.n, len(selected) + 1)
    assert k is not None
    scores = np.array(
        [forward_score(data, [*selected, candidate], k) for candidate in candidates]
    )
    risks = np.array(
        [local_risk_details(data, [*selected, candidate], k).risk for candidate in candidates]
    )

    assert np.all(scores < 0)
    assert int(np.argmax(scores)) == int(np.argmin(risks))


def test_prime_returns_zero_padded_effects() -> None:
    X, y = make_example()
    fit = prime(X, y)
    assert fit.selected_names == ("X1", "X2", "X5")
    assert set(fit.effects) == {f"X{index}" for index in range(1, 9)}
    assert all(fit.effects[name] == 0.0 for name in ("X3", "X4", "X6", "X7", "X8"))
    assert np.isclose(sum(fit.effects.values()), 1.0)
    assert len(fit.subset_values) == 8
    assert np.isclose(fit.risk_full, fit.noise_floor)
    for value in fit.subset_values:
        assert np.isclose(
            value.adjusted_risk,
            max(value.complement_risk - fit.noise_floor, 0.0),
        )
    np.testing.assert_allclose(
        list(fit.effects.values()),
        [
            0.438623583576,
            0.441997159276,
            0.0,
            0.0,
            0.119379257148,
            0.0,
            0.0,
            0.0,
        ],
        atol=1e-10,
    )


def test_plot_uses_the_fitted_effects() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    X, y = make_example()
    fit = prime(X, y)
    ax = prime_package.plot(fit)

    assert ax.get_xlabel() == "PRIME importance"
    assert ax.get_title(loc="left") == ""
    assert {label.get_text() for label in ax.get_yticklabels()} == set(fit.effects)
    plt.close(ax.figure)


def test_feature_names_and_query_indices_are_supported() -> None:
    X, y = make_example()
    names = [f"feature_{index}" for index in range(X.shape[1])]
    screen = prime_screen(
        X, y, feature_names=names, query_indices=range(0, X.shape[0], 2)
    )
    assert screen.selected_names == ("feature_0", "feature_1", "feature_4")
    assert screen.settings["query_count"] == 60


def test_subset_budget_is_explicit() -> None:
    X, y = make_example()
    try:
        prime(X, y, max_subsets=4)
    except RuntimeError as error:
        assert "max_subsets" in str(error)
    else:
        raise AssertionError("Expected the explicit subset budget to fail.")


@pytest.mark.parametrize("scale", [1e-100, 1e-8, 1.0, 1e8, 1e100])
def test_response_units_preserve_suppressor_selection(scale: float) -> None:
    location = np.repeat(np.arange(50, dtype=float), 2)
    signal = np.tile([-1.0, 1.0], 50)
    X = np.column_stack([location + 0.1 * signal, location - 0.1 * signal])
    fit = prime(X, scale * signal)
    assert fit.selected_names == ("X1", "X2")
    assert fit.settings["ridge_chosen_size"] == 2
    np.testing.assert_allclose(list(fit.effects.values()), [0.5, 0.5], rtol=1e-10)
    # Diagnostics retain squared response units even though BIC is dimensionless.
    np.testing.assert_allclose(fit.settings["ridge_path"][0]["rss"] / scale**2, 100.0)
    assert np.isfinite(fit.screen.backward_trace[0]["query_scale"])


@pytest.mark.parametrize("scale", [1e-100, 1e100])
def test_response_units_preserve_forward_calibration(scale: float) -> None:
    X, y = make_example()
    reference = prime(X, y)
    fit = prime(X, scale * y)
    assert fit.selected_names == reference.selected_names
    np.testing.assert_allclose(
        list(fit.effects.values()), list(reference.effects.values()), rtol=1e-10
    )
    assert all(np.isfinite(row["query_scale"]) for row in fit.screen.forward_trace)


def test_ridge_predictor_scaling_matches_r_reference() -> None:
    index = np.arange(1, 21, dtype=float)
    base = np.sin(index * np.sqrt(2))
    perturbation = 1e-4 * np.cos(index * np.sqrt(3))
    X = np.column_stack([base, base + perturbation])
    result = prime_screen(X, base + 10 * perturbation)
    # Fresh R MASS::lm.ridge reference, 2026-10-11. Near collinearity makes
    # the ddof=0 penalty convention observable; ddof=1 gives about 0.0566.
    np.testing.assert_allclose(
        list(result.ridge_weights.values()),
        [0.06647917453127333, 0.9335208254687266],
        atol=1e-7,
        rtol=0,
    )
    assert result.settings["ridge_chosen_size"] == 2
