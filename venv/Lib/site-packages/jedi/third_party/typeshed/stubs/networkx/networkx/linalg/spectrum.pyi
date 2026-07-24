import numpy as np
from networkx._typing import Array1D
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "laplacian_spectrum",
    "adjacency_spectrum",
    "modularity_spectrum",
    "normalized_laplacian_spectrum",
    "bethe_hessian_spectrum",
]

@_dispatchable
def laplacian_spectrum(G: Graph[_Node], weight: str | None = "weight") -> Array1D[np.float64]: ...
@_dispatchable
def normalized_laplacian_spectrum(G: Graph[_Node], weight: str | None = "weight") -> Array1D[np.float64]: ...
@_dispatchable
def adjacency_spectrum(G: Graph[_Node], weight: str | None = "weight") -> Array1D[np.complex128]: ...
@_dispatchable
def modularity_spectrum(G: Graph[_Node]) -> Array1D[np.complex128]: ...
@_dispatchable
def bethe_hessian_spectrum(G: Graph[_Node], r: float | None = None) -> Array1D[np.float64]: ...
