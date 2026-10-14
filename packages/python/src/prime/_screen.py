from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from ._core import (
    PreparedData,
    dimension_limit,
    local_risk_details,
    neighborhood_size,
    prepare_data,
    ridge_cap,
    settings_dict,
)
from ._types import PrimeScreen


def _ridge_initial(
    data: PreparedData,
) -> tuple[list[int], np.ndarray, int, int, tuple[dict[str, float | int], ...]]:
    cap = ridge_cap(data.n, data.p)
    # Unit mean square agrees with the manuscript and MASS::lm.ridge scaling.
    standardized = (data.X - np.mean(data.X, axis=0)) / np.std(data.X, axis=0, ddof=0)
    centered_y = data.y - np.mean(data.y)
    standardized_y = centered_y / np.sqrt(data.response_variance)
    cross = standardized.T @ standardized + data.ridge_lambda * np.eye(data.p)
    try:
        coefficients = np.linalg.solve(cross, standardized.T @ standardized_y)
    except np.linalg.LinAlgError:
        coefficients = np.linalg.lstsq(cross, standardized.T @ standardized_y, rcond=None)[
            0
        ]
    squared = coefficients**2
    total = float(np.sum(squared))
    weights = squared / total if total > 0 else np.zeros(data.p)
    ordered = sorted(range(data.p), key=lambda index: (-weights[index], index))

    path: list[dict[str, float | int]] = []
    for size in range(cap + 1):
        if size == 0:
            residual = standardized_y
            slope_rank = 0
        else:
            design = np.column_stack((np.ones(data.n), data.X[:, ordered[:size]]))
            fitted, _, rank, _ = np.linalg.lstsq(design, standardized_y, rcond=None)
            residual = standardized_y - design @ fitted
            slope_rank = max(0, int(rank) - 1)
        standardized_rss = float(residual @ residual)
        rss = standardized_rss * data.response_variance
        # A dimensionless floor preserves BIC ordering when response units change.
        # The omitted n*log(response_variance) is common to every candidate.
        bic = data.n * np.log(max(standardized_rss / data.n, np.finfo(float).eps))
        bic += slope_rank * np.log(data.n)
        path.append({"size": size, "rank": slope_rank, "rss": rss, "bic": float(bic)})
    chosen_size = min(path, key=lambda row: float(row["bic"]))["size"]
    selected = sorted(ordered[: int(chosen_size)])
    return selected, weights, cap, int(chosen_size), tuple(path)


def _forward(data: PreparedData) -> dict[str, Any]:
    ridge, ridge_weights, cap, ridge_chosen_size, ridge_path = _ridge_initial(data)
    selected = list(ridge)
    trace: list[dict[str, Any]] = []
    limit = dimension_limit(data.n, data.p)
    stop_reason = "all coordinates selected"

    while len(selected) < data.p:
        if len(selected) >= limit:
            stop_reason = "local-fit dimension limit reached"
            break
        candidates = [index for index in range(data.p) if index not in selected]
        dimension = len(selected) + 1
        k = neighborhood_size(data.n, dimension)
        if k is None:
            stop_reason = "local-fit dimension limit reached"
            break
        current_details = local_risk_details(data, selected, k)
        current = (
            (data.response_variance - current_details.risk) / data.response_variance
            if selected else 0.0
        )
        candidate_values = []
        best_details = None
        leading_value = -np.inf
        for candidate in candidates:
            details = local_risk_details(data, [*selected, candidate], k)
            value = (data.response_variance - details.risk) / data.response_variance
            candidate_values.append(value)
            # Preserve np.argmax's first-maximum (and first-NaN) convention.
            # Holding only the current leader avoids recomputing its local fits
            # without retaining contributions for every candidate.
            if best_details is None or (
                not np.isnan(leading_value) and (np.isnan(value) or value > leading_value)
            ):
                leading_value = value
                best_details = details
        best_position = int(np.argmax(candidate_values))
        added = candidates[best_position]
        best_value = float(candidate_values[best_position])
        gain = best_value - current
        assert best_details is not None
        augmented_details = best_details
        leader_key = (k, tuple(sorted([*selected, added])))
        data.detail_cache[leader_key] = augmented_details
        data.detail_cache.move_to_end(leader_key)
        while len(data.detail_cache) > 2:
            data.detail_cache.popitem(last=False)
        paired_gain = current_details.contribution - augmented_details.contribution
        normalized_gain = paired_gain / data.response_variance
        gain_mean = float(np.mean(normalized_gain))
        gain_scale = (
            float(np.std(normalized_gain, ddof=1) / np.sqrt(paired_gain.size))
            if paired_gain.size > 1
            else 0.0
        )
        gain_threshold = gain_scale * np.sqrt(2 * np.log(max(2, dimension)))
        calibrated = max(gain_mean - gain_threshold, 0.0)
        raw_gain = gain_mean * data.response_variance
        query_scale = gain_scale * data.response_variance
        threshold = gain_threshold * data.response_variance
        calibrated_gain = calibrated * data.response_variance
        accepted = bool(np.isfinite(calibrated) and calibrated > 0)
        trace.append(
            {
                "step": len(trace) + 1,
                "added": added,
                "added_name": data.feature_names[added],
                "selected_before": tuple(selected),
                "dimension": dimension,
                "k": k,
                "query_count": int(data.query_indices.size),
                "current_value": current,
                "candidate_value": best_value,
                "gain": gain,
                "raw_gain": raw_gain,
                "query_scale": query_scale,
                "threshold": float(threshold),
                "calibrated_gain": calibrated_gain,
                "accepted": accepted,
            }
        )
        if not accepted:
            stop_reason = "largest calibrated gain was nonpositive"
            break
        selected = sorted([*selected, added])

    candidates = list(selected)
    return {
        "candidates": candidates,
        "ridge": ridge,
        "ridge_weights": ridge_weights,
        "ridge_cap": cap,
        "ridge_chosen_size": ridge_chosen_size,
        "ridge_path": ridge_path,
        "local_path": [index for index in candidates if index not in ridge],
        "trace": trace,
        "stop_reason": stop_reason,
    }


def _backward(
    data: PreparedData, candidates: Sequence[int]
) -> tuple[list[int], list[dict[str, Any]]]:
    selected = sorted(set(int(index) for index in candidates))
    trace: list[dict[str, Any]] = []
    while selected:
        dimension = len(selected)
        k = neighborhood_size(data.n, dimension)
        if k is None:
            raise RuntimeError(
                f"Backward cleaning is infeasible at dimension {dimension}."
            )
        full = local_risk_details(data, selected, k)
        diagnostics: list[dict[str, Any]] = []
        calibrated: list[float] = []
        for variable in selected:
            reduced = [index for index in selected if index != variable]
            reduced_details = local_risk_details(data, reduced, k)
            paired_loss = reduced_details.contribution - full.contribution
            normalized_loss = paired_loss / data.response_variance
            loss_mean = float(np.mean(normalized_loss))
            loss_scale = (
                float(np.std(normalized_loss, ddof=1) / np.sqrt(paired_loss.size))
                if paired_loss.size > 1
                else 0.0
            )
            loss_threshold = loss_scale * np.sqrt(2 * np.log(max(2, dimension)))
            calibrated_loss = max(loss_mean - loss_threshold, 0.0)
            calibrated.append(calibrated_loss)
            raw_loss = loss_mean * data.response_variance
            query_scale = loss_scale * data.response_variance
            threshold = loss_threshold * data.response_variance
            diagnostics.append(
                {
                    "raw_loss": raw_loss,
                    "clipped_loss": max(raw_loss, 0.0),
                    "query_scale": query_scale,
                    "threshold": float(threshold),
                    "calibrated_loss": calibrated_loss * data.response_variance,
                }
            )
        tested_position = int(np.argmin(calibrated))
        tested = selected[tested_position]
        diagnostic = diagnostics[tested_position]
        remove = bool(
            np.isfinite(calibrated[tested_position])
            and calibrated[tested_position] == 0
        )
        trace.append(
            {
                "step": len(trace) + 1,
                "tested": tested,
                "tested_name": data.feature_names[tested],
                "selected_before": tuple(selected),
                "dimension": dimension,
                "k": k,
                "query_count": int(data.query_indices.size),
                "full_risk": full.risk,
                **diagnostic,
                "removed": remove,
            }
        )
        if not remove:
            break
        selected.pop(tested_position)
    return selected, trace


def _run_screen(
    X: ArrayLike,
    y: ArrayLike,
    *,
    feature_names: Sequence[str] | None,
    rescale: bool,
    query_indices: Sequence[int] | None,
    lambda_: float,
    neighbor_eps: float = 0.0,
    query_block_size: int = 256,
) -> tuple[PrimeScreen, PreparedData]:
    data = prepare_data(
        X,
        y,
        feature_names=feature_names,
        rescale=rescale,
        query_indices=query_indices,
        lambda_=lambda_,
        neighbor_eps=neighbor_eps,
        query_block_size=query_block_size,
    )
    forward = _forward(data)
    selected, backward_trace = _backward(data, forward["candidates"])
    ridge_weights = {
        name: float(weight)
        for name, weight in zip(
            data.feature_names, forward["ridge_weights"], strict=True
        )
    }
    screen = PrimeScreen(
        candidates=tuple(forward["candidates"]),
        candidate_names=tuple(
            data.feature_names[index] for index in forward["candidates"]
        ),
        selected=tuple(selected),
        selected_names=tuple(data.feature_names[index] for index in selected),
        forward_trace=tuple(forward["trace"]),
        backward_trace=tuple(backward_trace),
        ridge_weights=ridge_weights,
        settings={
            **settings_dict(data, rescale=rescale),
            "ridge_cap": forward["ridge_cap"],
            "ridge_chosen_size": forward["ridge_chosen_size"],
            "ridge_path": forward["ridge_path"],
            "forward_stop_reason": forward["stop_reason"],
        },
    )
    return screen, data


def prime_screen(
    X: ArrayLike,
    y: ArrayLike,
    *,
    feature_names: Sequence[str] | None = None,
    rescale: bool = True,
    query_indices: Sequence[int] | None = None,
    lambda_: float = 1e-6,
    neighbor_eps: float = 0.0,
    query_block_size: int = 256,
) -> PrimeScreen:
    """Select response-relevant coordinates using PRIME screening.

    ``neighbor_eps=0`` uses exact projected nearest neighbors. A positive value
    enables SciPy's approximate tree search and can change selection. All rows
    remain eligible neighbors and the neighborhood-size rule is unchanged.
    ``query_block_size`` limits queries handled together, without subsampling.
    """

    screen, _ = _run_screen(
        X,
        y,
        feature_names=feature_names,
        rescale=rescale,
        query_indices=query_indices,
        lambda_=lambda_,
        neighbor_eps=neighbor_eps,
        query_block_size=query_block_size,
    )
    return screen
