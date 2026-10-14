from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike

from ._core import PreparedData, local_risk, neighborhood_size
from ._pme import prime_allocate
from ._screen import _run_screen
from ._types import PrimeFit, SubsetValue


def _subset(
    mask: int, selected: Sequence[int], names: Sequence[str]
) -> tuple[str, ...]:
    return tuple(
        names[selected[position]]
        for position in range(len(selected))
        if mask & (1 << position)
    )


def _estimate_set_function(
    data: PreparedData, selected: Sequence[int], *, max_subsets: int | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    dimension = len(selected)
    subset_count = 1 << dimension
    if max_subsets is not None and subset_count > max_subsets:
        raise RuntimeError(
            f"The retained set requires {subset_count} subset values, "
            f"exceeding max_subsets={max_subsets}."
        )
    # Allocation only needs scalar risks, not selection's querywise contrasts.
    data.detail_cache.clear()
    k_full = neighborhood_size(data.n, dimension)
    if k_full is None:
        raise RuntimeError(
            "The retained set is infeasible for local-linear estimation."
        )
    noise_floor = local_risk(data, selected, k_full)
    signal = max(data.response_variance - noise_floor, 0.0)
    guard = np.sqrt(np.finfo(float).eps) * data.response_variance
    if not np.isfinite(signal) or signal <= guard:
        raise RuntimeError(
            "The retained coordinates have no estimable explained variance."
        )

    full_mask = subset_count - 1
    raw = np.zeros(subset_count, dtype=float)
    complement_risk = np.empty(subset_count, dtype=float)
    adjusted_risk = np.empty(subset_count, dtype=float)
    for mask in range(subset_count):
        keep = [
            selected[position]
            for position in range(dimension)
            if mask & (1 << position) == 0
        ]
        if keep:
            k = neighborhood_size(data.n, len(keep))
            if k is None:
                raise RuntimeError(
                    "A retained subset is infeasible for local-linear estimation."
                )
            risk_keep = local_risk(data, keep, k)
        else:
            risk_keep = data.response_variance
        complement_risk[mask] = risk_keep
        adjusted_risk[mask] = max(risk_keep - noise_floor, 0.0)
        raw[mask] = adjusted_risk[mask] / signal
    raw[0] = 0.0
    raw[full_mask] = 1.0
    return raw, complement_risk, adjusted_risk, signal, noise_floor


def prime(
    X: ArrayLike,
    y: ArrayLike,
    *,
    feature_names: Sequence[str] | None = None,
    rescale: bool = True,
    query_indices: Sequence[int] | None = None,
    lambda_: float = 1e-6,
    max_subsets: int | None = None,
    neighbor_eps: float = 0.0,
    query_block_size: int = 256,
) -> PrimeFit:
    """Run PRIME selection and proportional allocation.

    ``neighbor_eps=0`` keeps exact projected nearest-neighbor search. A positive
    value enables approximate search in screening and subset-risk estimation;
    selection and effects can change. Every row remains an eligible neighbor.
    ``query_block_size`` bounds query batches without reducing their number.
    Exact subset enumeration remains exponential in the selected dimension;
    ``max_subsets`` provides an explicit computation limit.
    """

    screen, data = _run_screen(
        X,
        y,
        feature_names=feature_names,
        rescale=rescale,
        query_indices=query_indices,
        lambda_=lambda_,
        neighbor_eps=neighbor_eps,
        query_block_size=query_block_size,
    )
    full_effects = {name: 0.0 for name in data.feature_names}
    if not screen.selected:
        return PrimeFit(
            effects=full_effects,
            reduced_effects={},
            selected=(),
            selected_names=(),
            screen=screen,
            subset_values=(),
            signal=0.0,
            noise_floor=data.response_variance,
            risk_full=data.response_variance,
            projection=None,
            settings=screen.settings,
        )

    raw, complement_risk, adjusted_risk, signal, noise_floor = _estimate_set_function(
        data, screen.selected, max_subsets=max_subsets
    )
    allocation = prime_allocate(raw)
    reduced = {
        data.feature_names[index]: float(effect)
        for index, effect in zip(screen.selected, allocation.effects, strict=True)
    }
    full_effects.update(reduced)
    subset_values = tuple(
        SubsetValue(
            mask=mask,
            subset=_subset(mask, screen.selected, data.feature_names),
            raw=float(raw[mask]),
            projected=float(allocation.projected_values[mask]),
            complement_risk=float(complement_risk[mask]),
            adjusted_risk=float(adjusted_risk[mask]),
        )
        for mask in range(raw.size)
    )
    return PrimeFit(
        effects=full_effects,
        reduced_effects=reduced,
        selected=screen.selected,
        selected_names=screen.selected_names,
        screen=screen,
        subset_values=subset_values,
        signal=float(signal),
        noise_floor=float(noise_floor),
        risk_full=float(noise_floor),
        projection=allocation.projection,
        settings=screen.settings,
    )
