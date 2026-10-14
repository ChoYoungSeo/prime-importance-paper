from __future__ import annotations

from typing import Any

import numpy as np

from ._types import PrimeFit


def plot(
    fit: PrimeFit,
    *,
    max_display: int = 20,
    ax: Any = None,
) -> Any:
    """Plot the proportional effects in a fitted PRIME result."""

    if not isinstance(fit, PrimeFit):
        raise TypeError("fit must be a PrimeFit object returned by prime.fit().")
    if not isinstance(max_display, int) or max_display < 1:
        raise ValueError("max_display must be a positive integer.")
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "Plotting requires matplotlib; install prime-importance[plot]."
        ) from error

    names = np.asarray(list(fit.effects), dtype=object)
    effects = np.asarray(list(fit.effects.values()), dtype=float)
    order = np.argsort(-effects, kind="stable")[: min(max_display, effects.size)]
    names = names[order]
    effects = effects[order]
    selected = np.isin(names, fit.selected_names)

    created_axis = ax is None
    if created_axis:
        height = max(2.5, 0.30 * len(names) + 1.0)
        _, ax = plt.subplots(figsize=(5.8, height))

    teal = "#087F8C"
    gray = "#A5ADB4"
    light_gray = "#E4E8EB"
    ink = "#20262B"
    positions = np.arange(len(names))

    ax.hlines(positions, 0, effects, color=light_gray, linewidth=1.5, zorder=1)
    ax.scatter(
        effects,
        positions,
        s=np.where(selected, 48, 30),
        color=np.where(selected, teal, gray),
        edgecolor="white",
        linewidth=0.7,
        zorder=2,
    )
    ax.set_yticks(positions, names)
    ax.invert_yaxis()
    ax.set_xlabel("PRIME importance", color=ink)

    scale = max(float(np.max(effects, initial=0)), 0.1)
    ax.set_xlim(-0.025 * scale, 1.12 * scale)
    ax.grid(axis="x", color=light_gray, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0, colors=ink)
    ax.tick_params(axis="x", colors=ink)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(gray)
    for label, is_selected in zip(ax.get_yticklabels(), selected, strict=True):
        label.set_color(ink if is_selected else gray)

    if created_axis:
        ax.figure.tight_layout()
    return ax
