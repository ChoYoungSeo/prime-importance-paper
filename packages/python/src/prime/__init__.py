"""PRIME: proportional global variable importance under dependent inputs."""

from ._plot import plot
from ._pme import prime_allocate
from ._prime import prime
from ._screen import prime_screen
from ._types import (
    PrimeAllocation,
    PrimeFit,
    PrimeScreen,
    ProjectionDiagnostics,
    SubsetValue,
)

fit = prime
screen = prime_screen
allocate = prime_allocate

__all__ = [
    "PrimeAllocation",
    "PrimeFit",
    "PrimeScreen",
    "ProjectionDiagnostics",
    "SubsetValue",
    "allocate",
    "fit",
    "prime_allocate",
    "plot",
    "prime",
    "prime_screen",
    "screen",
]

__version__ = "0.2.1"
