from itertools import combinations

import numpy as np
import pytest
from scipy.spatial import cKDTree

from prime import _core, prime, prime_screen
from prime._core import local_risk, local_risk_details, neighbor_blocks, prepare_data
from prime._prime import _estimate_set_function


def make_data(n=513, p=5, **kwargs):
    rng = np.random.default_rng(20260913)
    X = rng.normal(size=(n, p))
    y = X[:, 0] + X[:, 1] ** 2 + rng.normal(0, 0.1, size=n)
    return prepare_data(
        X, y, feature_names=None, rescale=True, query_indices=None, lambda_=1e-6, **kwargs
    )


def dense_neighbors(X, queries, k):
    # Pre-streaming reference: one exact query over the full requested query set.
    _, indices = cKDTree(X).query(X[queries], k=k + 1, workers=1)
    return np.array([row[row != query][:k] for row, query in zip(indices, queries, strict=True)])


def dense_contributions(data, subset, k):
    X = data.X[:, subset]
    neighbors = dense_neighbors(X, data.query_indices, k)
    local_sse = np.empty(data.query_indices.size)
    for row, query in enumerate(data.query_indices):
        index = neighbors[row]
        design = np.column_stack((np.ones(k), X[index] - X[query]))
        cross = design.T @ design
        cross[1:, 1:] += data.lambda_ * np.eye(len(subset))
        coefficient = np.linalg.solve(cross, design.T @ data.y[index])
        residuals = data.y[index] - design @ coefficient
        local_sse[row] = float(residuals @ residuals)
    return local_sse / (k - len(subset) - 1)


@pytest.mark.parametrize("block_size", [1, 64, 256, 1024])
def test_exact_blocks_preserve_dense_neighbors_and_risks(block_size):
    data = make_data(query_block_size=block_size)
    data.query_indices = np.arange(data.n - 1, -1, -2)
    expected = dense_contributions(data, [0, 1, 3], 35)
    actual = local_risk_details(data, [0, 1, 3], 35)
    np.testing.assert_array_equal(actual.contribution, expected)
    assert actual.risk == float(np.mean(expected))
    assert local_risk(data, [0, 1, 3], 35) == actual.risk


def test_duplicate_and_equidistant_rows_keep_exact_tree_tie_behavior():
    grid = np.array(list(np.ndindex(2, 2, 2, 2)), dtype=float)
    X = np.tile(grid, (20, 1))
    queries = np.arange(len(X) - 1, -1, -1)
    expected = dense_neighbors(X, queries, 21)
    actual = np.concatenate([value for _, value in neighbor_blocks(X, queries, 21)])
    np.testing.assert_array_equal(actual, expected)
    assert np.all(actual != queries[:, None])


def test_streaming_builds_one_tree_uses_every_query_and_all_neighbor_rows(monkeypatch):
    constructed = []
    queried = []

    class ObservedTree:
        def __init__(self, X):
            constructed.append(len(X))
            self.tree = cKDTree(X)

        def query(self, X, **kwargs):
            queried.append((len(X), kwargs["eps"]))
            return self.tree.query(X, **kwargs)

    monkeypatch.setattr(_core, "cKDTree", ObservedTree)
    data = make_data(query_block_size=32, neighbor_eps=0.2)
    local_risk(data, [0, 1], 35)
    assert constructed == [513]
    assert sum(count for count, _ in queried) == 513
    assert max(count for count, _ in queried) <= 32
    assert all(eps == 0.2 for _, eps in queried)


def test_scalar_and_paired_caches_stay_separate_and_bounded():
    data = make_data(n=80, p=4)
    subsets = [subset for d in range(1, 5) for subset in combinations(range(4), d)]
    for subset in subsets:
        local_risk(data, subset, 20)
    assert len(data.cache) == 15
    assert all(isinstance(value, float) for value in data.cache.values())
    assert not data.detail_cache
    for subset in subsets:
        details = local_risk_details(data, subset, 20)
        assert len(data.detail_cache) <= 2
        assert details.risk == data.cache[(20, subset)]
    _estimate_set_function(data, [0, 1, 2, 3], max_subsets=None)
    assert not data.detail_cache
    assert all(isinstance(value, float) for value in data.cache.values())


def test_approximate_neighbors_satisfy_distance_bound_and_exclude_self():
    data = make_data()
    queries = data.query_indices
    k = 20
    exact = dense_neighbors(data.X, queries, k)
    approx = np.concatenate([
        value for _, value in neighbor_blocks(data.X, queries, k, neighbor_eps=0.5)
    ])
    exact_radius = np.linalg.norm(data.X[exact[:, -1]] - data.X, axis=1)
    approx_radius = np.linalg.norm(data.X[approx[:, -1]] - data.X, axis=1)
    assert np.all(approx != queries[:, None])
    assert np.all(approx_radius <= 1.5 * exact_radius + 1e-12)


def test_public_block_sizes_preserve_selection_and_allocation():
    data = make_data(n=90, p=4)
    small = prime(data.X, data.y, query_block_size=7)
    full = prime(data.X, data.y, query_block_size=1000)
    assert small.selected == full.selected
    assert small.effects == full.effects
    assert small.subset_values == full.subset_values
    assert small.screen.forward_trace == full.screen.forward_trace
    assert small.screen.backward_trace == full.screen.backward_trace


def test_approximation_is_explicit_in_fit_and_screen_settings():
    data = make_data(n=90, p=4)
    exact = prime_screen(data.X, data.y)
    fit = prime(data.X, data.y, neighbor_eps=0.2, query_block_size=23)
    assert exact.settings["neighbor_search"] == "exact"
    assert exact.settings["neighbor_eps"] == 0
    assert fit.settings["neighbor_search"] == "approximate"
    assert fit.settings["neighbor_eps"] == 0.2
    assert fit.settings["query_block_size"] == 23
    assert fit.settings["query_count"] == fit.settings["neighbor_pool_size"] == 90
    assert np.isclose(sum(fit.effects.values()), 1)


@pytest.mark.parametrize("value", [-0.1, np.nan, np.inf, True, [], "0.1"])
def test_invalid_neighbor_eps_is_rejected(value):
    data = make_data(n=20)
    with pytest.raises(ValueError, match="neighbor_eps"):
        prime_screen(data.X, data.y, neighbor_eps=value)


@pytest.mark.parametrize("value", [0, -1, 1.5, np.nan, True, "256"])
def test_invalid_query_block_size_is_rejected(value):
    data = make_data(n=20)
    with pytest.raises(ValueError, match="query_block_size"):
        prime(data.X, data.y, query_block_size=value)
