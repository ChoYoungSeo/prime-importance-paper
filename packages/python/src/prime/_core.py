from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class LocalRiskDetails:
    risk: float
    contribution: FloatArray


@dataclass
class PreparedData:
    X: FloatArray
    y: FloatArray
    feature_names: tuple[str, ...]
    response_variance: float
    center: FloatArray
    scale: FloatArray
    query_indices: IntArray
    lambda_: float
    neighbor_eps: float = 0.0
    query_block_size: int = 256
    ridge_lambda: float = 1e-6
    cache: dict[tuple[int, tuple[int, ...]], float] = field(default_factory=dict)
    detail_cache: OrderedDict[tuple[int, tuple[int, ...]], LocalRiskDetails] = field(
        default_factory=OrderedDict
    )

    @property
    def n(self) -> int:
        return int(self.X.shape[0])

    @property
    def p(self) -> int:
        return int(self.X.shape[1])


def prepare_data(
    X: ArrayLike,
    y: ArrayLike,
    *,
    feature_names: Sequence[str] | None,
    rescale: bool,
    query_indices: Sequence[int] | None,
    lambda_: float,
    neighbor_eps: float = 0.0,
    query_block_size: int = 256,
) -> PreparedData:
    inferred_names: Sequence[str] | None = None
    if hasattr(X, "columns") and hasattr(X, "to_numpy"):
        inferred_names = tuple(str(value) for value in X.columns)
        X_array = np.asarray(X.to_numpy(), dtype=float)
    else:
        X_array = np.asarray(X, dtype=float)
    y_array = np.asarray(y, dtype=float)

    if X_array.ndim != 2 or X_array.shape[0] < 10 or X_array.shape[1] < 1:
        raise ValueError("X must have at least 10 rows and one column.")
    if y_array.ndim != 1 or y_array.size != X_array.shape[0]:
        raise ValueError("y must be one-dimensional and match the rows of X.")
    if not np.all(np.isfinite(X_array)) or not np.all(np.isfinite(y_array)):
        raise ValueError("X and y must contain only finite values.")
    scales = np.std(X_array, axis=0, ddof=1)
    if np.any(~np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("X must not contain constant columns.")
    response_variance = float(np.var(y_array, ddof=1))
    if not np.isfinite(response_variance) or response_variance <= 0:
        raise ValueError("y must be nonconstant.")
    if not isinstance(rescale, (bool, np.bool_)):
        raise TypeError("rescale must be a boolean.")
    if not np.isfinite(lambda_) or lambda_ < 0:
        raise ValueError("lambda_ must be finite and nonnegative.")
    if (
        isinstance(neighbor_eps, (bool, np.bool_))
        or not np.isscalar(neighbor_eps)
        or not isinstance(neighbor_eps, (int, float, np.integer, np.floating))
        or not np.isfinite(neighbor_eps)
        or neighbor_eps < 0
    ):
        raise ValueError("neighbor_eps must be a finite nonnegative number.")
    if (
        isinstance(query_block_size, (bool, np.bool_))
        or not isinstance(query_block_size, (int, np.integer))
        or query_block_size < 1
    ):
        raise ValueError("query_block_size must be a positive integer.")

    if feature_names is not None:
        names = tuple(feature_names)
    elif inferred_names is not None:
        names = tuple(inferred_names)
    else:
        names = tuple(f"X{i + 1}" for i in range(X_array.shape[1]))
    if (
        len(names) != X_array.shape[1]
        or not all(isinstance(name, str) for name in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError(
            "feature_names must be unique strings and match the columns of X."
        )

    center = np.zeros(X_array.shape[1], dtype=float)
    scale = np.ones(X_array.shape[1], dtype=float)
    if rescale:
        center = np.mean(X_array, axis=0)
        scale = scales
        X_array = (X_array - center) / scale
    else:
        X_array = X_array.copy()

    if query_indices is None:
        queries = np.arange(X_array.shape[0], dtype=np.int64)
    else:
        # Preserve booleans in mixed sequences until their types are checked.
        raw_queries = np.asarray(query_indices, dtype=object)
        if raw_queries.ndim != 1 or raw_queries.size == 0:
            raise ValueError(
                "query_indices must be a nonempty one-dimensional sequence."
            )
        for query in raw_queries:
            if (
                isinstance(query, (bool, np.bool_))
                or not isinstance(query, (int, float, np.integer, np.floating))
                or (
                    isinstance(query, (float, np.floating))
                    and (not np.isfinite(query) or query != np.floor(query))
                )
                or query < 0
                or query >= X_array.shape[0]
            ):
                raise ValueError("query_indices contains an invalid row index.")
        queries = np.asarray(raw_queries, dtype=np.int64)
        if np.unique(queries).size != queries.size:
            raise ValueError("query_indices must not contain duplicates.")

    return PreparedData(
        X=X_array,
        y=y_array.copy(),
        feature_names=names,
        response_variance=response_variance,
        center=center,
        scale=scale,
        query_indices=queries,
        lambda_=float(lambda_),
        neighbor_eps=float(neighbor_eps),
        query_block_size=int(query_block_size),
    )


def dimension_limit(n: int, p: int) -> int:
    return int(min(p, np.floor(n * 0.5 / 3.0)))


def ridge_cap(n: int, p: int) -> int:
    target = max(1, int(np.ceil(1.5 * np.log(max(2, p)))))
    return min(p, target, dimension_limit(n, p))


def neighborhood_size(n: int, dimension: int) -> int | None:
    k = int(np.ceil(max(1.5 * np.sqrt(n), 3 * dimension)))
    cap = min(n - 1, int(np.floor(0.5 * n)))
    if k <= dimension + 1 or k > cap:
        return None
    return k


def neighbor_blocks(
    X: FloatArray,
    queries: IntArray,
    k: int,
    *,
    neighbor_eps: float = 0.0,
    query_block_size: int = 256,
) -> Iterator[tuple[IntArray, IntArray]]:
    """Build one projected tree and yield non-self neighbors for each query block."""
    tree = cKDTree(X)
    for start in range(0, queries.size, query_block_size):
        block = queries[start : start + query_block_size]
        _, raw = tree.query(
            X[block], k=min(k + 1, X.shape[0]), eps=neighbor_eps, workers=1
        )
        raw = np.atleast_2d(raw)
        neighbors = np.empty((block.size, k), dtype=np.int64)
        for row, query in enumerate(block):
            candidates = raw[row]
            candidates = candidates[candidates != query]
            if candidates.size < k:
                distances = np.sum((X - X[query]) ** 2, axis=1)
                distances[query] = np.inf
                candidates = np.argsort(distances, kind="stable")[:k]
            neighbors[row] = candidates[:k]
        yield block, neighbors


def neighbor_indices(X: FloatArray, queries: IntArray, k: int) -> IntArray:
    """Materialize exact neighbors for callers that explicitly need the full matrix."""
    return np.concatenate([neighbors for _, neighbors in neighbor_blocks(X, queries, k)])


def local_risk_details(
    data: PreparedData, subset: Sequence[int], k: int
) -> LocalRiskDetails:
    key_subset = tuple(sorted(set(int(value) for value in subset)))
    if not key_subset:
        return LocalRiskDetails(
            risk=data.response_variance,
            contribution=np.full(
                data.query_indices.size, data.response_variance, dtype=float
            ),
        )
    key = (int(k), key_subset)
    if key in data.detail_cache:
        data.detail_cache.move_to_end(key)
        return data.detail_cache[key]
    details = _evaluate_local_risk(data, key_subset, k)
    data.cache[key] = details.risk
    data.detail_cache[key] = details
    # Paired selection comparisons need two querywise vectors. Keep scalar risks
    # separately so neither candidate ranking nor subset enumeration retains n_Q
    # contributions for every evaluated subset.
    while len(data.detail_cache) > 2:
        data.detail_cache.popitem(last=False)
    return details


def _evaluate_local_risk(
    data: PreparedData, key_subset: tuple[int, ...], k: int
) -> LocalRiskDetails:
    dimension = len(key_subset)
    residual_df = k - dimension - 1
    if residual_df <= 0:
        raise ValueError("Neighborhood size is too small for the local-linear fit.")

    Xs = data.X[:, key_subset]
    local_sse = np.empty(data.query_indices.size, dtype=float)
    intercept = np.ones(k)
    slope_penalty = data.lambda_ * np.eye(dimension)
    offset = 0
    for queries, neighbors in neighbor_blocks(
        Xs,
        data.query_indices,
        k,
        neighbor_eps=data.neighbor_eps,
        query_block_size=data.query_block_size,
    ):
        for row, query in enumerate(queries):
            index = neighbors[row]
            centered = Xs[index] - Xs[query]
            design = np.column_stack((intercept, centered))
            cross = design.T @ design
            cross[1:, 1:] += slope_penalty
            rhs = design.T @ data.y[index]
            try:
                coefficients = np.linalg.solve(cross, rhs)
            except np.linalg.LinAlgError:
                coefficients = np.linalg.lstsq(cross, rhs, rcond=None)[0]
            residuals = data.y[index] - design @ coefficients
            local_sse[offset + row] = float(residuals @ residuals)
        offset += queries.size
    contribution = local_sse / residual_df
    return LocalRiskDetails(
        risk=float(np.mean(contribution)),
        contribution=contribution,
    )


def local_risk(data: PreparedData, subset: Sequence[int], k: int) -> float:
    key_subset = tuple(sorted(set(int(value) for value in subset)))
    if not key_subset:
        return data.response_variance
    key = (int(k), key_subset)
    if key not in data.cache:
        # The temporary vector preserves the original NumPy mean reduction.
        # It is discarded here; only paired comparisons retain contributions.
        data.cache[key] = _evaluate_local_risk(data, key_subset, k).risk
    return data.cache[key]


def forward_score(data: PreparedData, subset: Sequence[int], k: int) -> float:
    if len(subset) == 0:
        return 0.0
    risk = local_risk(data, subset, k)
    return (data.response_variance - risk) / data.response_variance


def settings_dict(data: PreparedData, *, rescale: bool) -> dict[str, Any]:
    return {
        "ridge_rule": (
            "BIC-selected size along the standardized ridge-weight ranking, "
            "searched through ceil(1.5*log(p)) coordinates"
        ),
        "ridge_lambda": data.ridge_lambda,
        "forward_rule": (
            "start from the ridge seed and extend the local-linear path while "
            "the leading query-calibrated gain is positive"
        ),
        "backward_rule": (
            "sequentially delete zero calibrated losses using paired querywise "
            "scales and the sqrt(2*log(d)) threshold"
        ),
        "neighborhood_rule": "ceil(max(1.5*sqrt(n), 3*d))",
        "query_count": int(data.query_indices.size),
        "neighbor_pool_size": data.n,
        "neighbor_search": "exact" if data.neighbor_eps == 0 else "approximate",
        "neighbor_eps": data.neighbor_eps,
        "query_block_size": data.query_block_size,
        "lambda": data.lambda_,
        "rescale": rescale,
    }
