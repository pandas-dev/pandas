from collections.abc import Collection
from typing import Literal

import numpy as np
from networkx._typing import Array2D
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from scipy.sparse import csr_array  # type: ignore[import-untyped]  # pyright: ignore[reportMissingImports]

__all__ = [
    "laplacian_matrix",
    "normalized_laplacian_matrix",
    "directed_laplacian_matrix",
    "directed_combinatorial_laplacian_matrix",
]

@_dispatchable
def laplacian_matrix(G: Graph[_Node], nodelist: Collection[_Node] | None = None, weight: str | None = "weight") -> csr_array: ...
@_dispatchable
def normalized_laplacian_matrix(
    G: Graph[_Node], nodelist: Collection[_Node] | None = None, weight: str | None = "weight"
) -> csr_array: ...
@_dispatchable
def directed_laplacian_matrix(
    G: DiGraph[_Node],
    nodelist: Collection[_Node] | None = None,
    weight: str | None = "weight",
    walk_type: Literal["random", "lazy", "pagerank"] | None = None,
    alpha: float = 0.95,
) -> Array2D[np.float64]: ...
@_dispatchable
def directed_combinatorial_laplacian_matrix(
    G: DiGraph[_Node],
    nodelist: Collection[_Node] | None = None,
    weight: str | None = "weight",
    walk_type: Literal["random", "lazy", "pagerank"] | None = None,
    alpha: float = 0.95,
) -> Array2D[np.float64]: ...
