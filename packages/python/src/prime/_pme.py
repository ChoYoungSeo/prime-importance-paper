from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._types import PrimeAllocation, ProjectionDiagnostics

FloatArray = NDArray[np.float64]


def _mask_size(mask: int) -> int:
    return mask.bit_count()


def _hasse_edges(dimension: int) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    lower: list[int] = []
    upper: list[int] = []
    full = (1 << dimension) - 1
    for mask in range(full + 1):
        for position in range(dimension):
            bit = 1 << position
            if mask & bit == 0:
                lower.append(mask)
                upper.append(mask | bit)
    return np.asarray(lower, dtype=np.int64), np.asarray(upper, dtype=np.int64)


def _project_monotone(
    game: FloatArray,
    dimension: int,
    tolerance: float = 1e-10,
    max_iterations: int = 20_000,
) -> tuple[FloatArray, ProjectionDiagnostics]:
    full = (1 << dimension) - 1
    full_value = float(game[full])
    if not np.isfinite(full_value) or full_value < 0:
        raise ValueError("The full-set value must be finite and nonnegative.")
    if dimension == 1:
        projected = np.asarray([0.0, full_value], dtype=float)
        diagnostics = ProjectionDiagnostics(0, True, 0.0)
        return projected, diagnostics

    interior_masks = np.arange(1, full, dtype=np.int64)
    x = game[interior_masks].astype(float, copy=True)
    variable_index = np.full(full + 1, -1, dtype=np.int64)
    variable_index[interior_masks] = np.arange(interior_masks.size)
    lower, upper = _hasse_edges(dimension)
    lower_index = variable_index[lower]
    upper_index = variable_index[upper]
    norm_squared = (lower_index >= 0).astype(float) + (upper_index >= 0).astype(float)
    rhs = np.where(upper == full, -full_value, 0.0)
    correction = np.zeros(lower.size, dtype=float)

    projected = np.empty(full + 1, dtype=float)
    max_violation = math.inf
    iterations = 0
    for _ in range(max_iterations):
        iterations += 1
        max_move = 0.0
        for edge in range(lower.size):
            li = int(lower_index[edge])
            ui = int(upper_index[edge])
            dot_current = 0.0
            if ui >= 0:
                dot_current += x[ui]
            if li >= 0:
                dot_current -= x[li]
            dot_shifted = dot_current + correction[edge] * norm_squared[edge]
            delta = (
                (rhs[edge] - dot_shifted) / norm_squared[edge]
                if dot_shifted < rhs[edge]
                else 0.0
            )
            move = correction[edge] + delta
            if ui >= 0:
                x[ui] += move
            if li >= 0:
                x[li] -= move
            correction[edge] = -delta
            max_move = max(max_move, abs(move))

        projected.fill(0.0)
        projected[interior_masks] = x
        projected[full] = full_value
        violations = projected[lower] - projected[upper]
        max_violation = max(0.0, float(np.max(violations)))
        if max_move <= tolerance and max_violation <= tolerance:
            break
    else:
        raise RuntimeError("Monotone set-function projection did not converge.")

    projected = np.clip(projected, 0.0, full_value)
    projected[0] = 0.0
    projected[full] = full_value
    diagnostics = ProjectionDiagnostics(iterations, True, max_violation)
    return projected, diagnostics


def _log_ratio_potential_shifted(
    game: FloatArray, dimension: int, zero_mask: int
) -> tuple[FloatArray, int]:
    full = (1 << dimension) - 1
    remaining = full ^ zero_mask
    log_potential = np.full(full + 1, -np.inf, dtype=float)
    log_potential[0] = 0.0
    admissible = [mask for mask in range(full + 1) if mask & zero_mask == 0]
    for order in range(1, _mask_size(remaining) + 1):
        for mask in admissible:
            if _mask_size(mask) != order:
                continue
            predecessors = [mask ^ (1 << pos) for pos in range(dimension) if mask & (1 << pos)]
            log_denominator = np.logaddexp.reduce(-log_potential[predecessors])
            shifted_value = game[zero_mask | mask]
            log_potential[mask] = math.log(shifted_value) - log_denominator
    return log_potential, remaining


def _extended_proportional_value(game: FloatArray, dimension: int) -> FloatArray:
    full = (1 << dimension) - 1
    if game[full] == 0:
        return np.zeros(dimension, dtype=float)
    zero_masks = [mask for mask in range(full) if game[mask] == 0]
    largest_size = max(_mask_size(mask) for mask in zero_masks)
    largest_zero = [mask for mask in zero_masks if _mask_size(mask) == largest_size]

    log_denominator = -np.inf
    log_numerator = np.full(dimension, -np.inf, dtype=float)
    for zero_mask in largest_zero:
        log_potential, remaining = _log_ratio_potential_shifted(game, dimension, zero_mask)
        log_denominator = np.logaddexp(log_denominator, -log_potential[remaining])
        for position in range(dimension):
            bit = 1 << position
            if zero_mask & bit == 0:
                without = remaining ^ bit
                log_numerator[position] = np.logaddexp(
                    log_numerator[position], -log_potential[without]
                )
    effects = np.exp(log_numerator - log_denominator)
    # Remove rounding drift in the log-domain calculation without truncating
    # small positive effects. The game passed here has full-set value one.
    return effects / np.sum(effects)


def prime_allocate(
    subset_values: ArrayLike,
    *,
    zero_tolerance: float = math.sqrt(np.finfo(float).eps),
) -> PrimeAllocation:
    """Allocate precomputed subset values using the PRIME operator.

    Parameters
    ----------
    subset_values:
        Values ordered by binary subset masks from 0 through ``2**d - 1``.
        The returned importance values sum to the projected full-set value.
    zero_tolerance:
        Nonnegative threshold relative to the full-set value, applied only to
        proper subsets after projection. Zero preserves every positive value.
        The empty/full endpoints stay fixed; values at least one make every
        proper subset null.

    Notes
    -----
    Positive-full games are normalized internally and the result is returned
    in the input units. A zero full-set value gives an all-zero allocation.
    A game whose normalized values cannot be represented as finite, nonzero
    float64 values raises ``ValueError`` rather than changing its null sets.
    """

    values = np.asarray(subset_values, dtype=float)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ValueError("subset_values must be a finite one-dimensional array.")
    dimension = int(round(math.log2(values.size)))
    if 1 << dimension != values.size:
        raise ValueError("subset_values length must equal 2**d for an integer d >= 1.")
    if (
        isinstance(zero_tolerance, (bool, np.bool_))
        or not isinstance(zero_tolerance, (int, float, np.integer, np.floating))
        or not np.isfinite(zero_tolerance)
        or zero_tolerance < 0
    ):
        raise ValueError("zero_tolerance must be a finite nonnegative number.")

    values = values.copy()
    values[0] = 0.0
    full_value = float(values[-1])
    if full_value < 0:
        raise ValueError("The full-set value must be finite and nonnegative.")
    if full_value == 0:
        return PrimeAllocation(
            np.zeros(dimension), np.zeros(values.size), ProjectionDiagnostics(0, True, 0.0)
        )
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        normalized = values / full_value
    if np.any(~np.isfinite(normalized)) or np.any((values != 0) & (normalized == 0)):
        raise ValueError("The game has too much dynamic range for float64 normalization.")
    projected, diagnostics = _project_monotone(normalized, dimension)
    interior = projected[1:-1]
    interior[interior <= zero_tolerance] = 0.0
    effects = _extended_proportional_value(projected, dimension)
    diagnostics = ProjectionDiagnostics(
        diagnostics.iterations, diagnostics.converged, diagnostics.max_violation * full_value
    )
    return PrimeAllocation(effects * full_value, projected * full_value, diagnostics)
