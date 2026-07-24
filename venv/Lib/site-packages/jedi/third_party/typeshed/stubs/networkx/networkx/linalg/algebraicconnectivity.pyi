from typing import Literal

import numpy as np
from networkx._typing import Array1D, Seed
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["algebraic_connectivity", "fiedler_vector", "spectral_ordering", "spectral_bisection"]

@_dispatchable
def algebraic_connectivity(
    G: Graph[_Node],
    weight: str | None = "weight",
    normalized: bool = False,
    tol: float = 1e-08,
    method: Literal["tracemin_pcg", "tracemin_lu", "lanczos", "lobpcg"] = "tracemin_pcg",
    seed: Seed | None = None,
) -> float: ...
@_dispatchable
def fiedler_vector(
    G: Graph[_Node],
    weight: str | None = "weight",
    normalized: bool = False,
    tol: float = 1e-08,
    method: Literal["tracemin_pcg", "tracemin_lu", "lanczos", "lobpcg"] = "tracemin_pcg",
    seed: Seed | None = None,
) -> Array1D[np.float64]: ...
@_dispatchable
def spectral_ordering(
    G: Graph[_Node],
    weight: str | None = "weight",
    normalized: bool = False,
    tol: float = 1e-08,
    method: Literal["tracemin_pcg", "tracemin_lu", "lanczos", "lobpcg"] = "tracemin_pcg",
    seed: Seed | None = None,
) -> list[_Node]: ...
@_dispatchable
def spectral_bisection(
    G: Graph[_Node],
    weight: str | None = "weight",
    normalized: bool = False,
    tol: float = 1e-08,
    method: Literal["tracemin_pcg", "tracemin_lu", "lanczos", "lobpcg"] = "tracemin_pcg",
    seed: Seed | None = None,
) -> tuple[set[_Node], set[_Node]]: ...
