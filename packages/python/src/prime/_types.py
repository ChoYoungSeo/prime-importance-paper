from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ProjectionDiagnostics:
    iterations: int
    converged: bool
    max_violation: float


@dataclass(frozen=True)
class PrimeAllocation:
    effects: FloatArray
    projected_values: FloatArray
    projection: ProjectionDiagnostics

    def __array__(
        self, dtype: Any = None, copy: bool | None = None
    ) -> NDArray[Any]:
        if copy is False:
            if dtype is not None and np.dtype(dtype) != self.effects.dtype:
                raise ValueError("A dtype conversion requires a copy of the effects.")
            return self.effects
        if copy is True:
            return np.array(self.effects, dtype=dtype, copy=True)
        return np.asarray(self.effects, dtype=dtype)


@dataclass(frozen=True)
class SubsetValue:
    mask: int
    subset: tuple[str, ...]
    raw: float
    projected: float
    complement_risk: float
    adjusted_risk: float


@dataclass(frozen=True)
class PrimeScreen:
    candidates: tuple[int, ...]
    candidate_names: tuple[str, ...]
    selected: tuple[int, ...]
    selected_names: tuple[str, ...]
    forward_trace: tuple[Mapping[str, Any], ...]
    backward_trace: tuple[Mapping[str, Any], ...]
    ridge_weights: Mapping[str, float]
    settings: Mapping[str, Any]

    def __repr__(self) -> str:
        names = ", ".join(self.selected_names) or "<empty>"
        return f"PrimeScreen(selected=[{names}])"


@dataclass(frozen=True)
class PrimeFit:
    effects: Mapping[str, float]
    reduced_effects: Mapping[str, float]
    selected: tuple[int, ...]
    selected_names: tuple[str, ...]
    screen: PrimeScreen
    subset_values: tuple[SubsetValue, ...]
    signal: float | None
    noise_floor: float | None
    risk_full: float | None
    projection: ProjectionDiagnostics | None
    settings: Mapping[str, Any]

    def __repr__(self) -> str:
        names = ", ".join(self.selected_names) or "<empty>"
        return f"PrimeFit(selected=[{names}], effects={dict(self.effects)!r})"
