from collections.abc import Collection, Hashable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.typing import DTypeLike
from scipy.sparse import csc_array, csr_array  # type: ignore[import-untyped]  # pyright: ignore[reportMissingImports]

__all__ = ["incidence_matrix", "adjacency_matrix"]

@_dispatchable
def incidence_matrix(
    G: Graph[_Node],
    nodelist: Collection[_Node] | None = None,
    edgelist: (
        Collection[
            # Requiring tuples to represent an edge might be too strict as runtime does not check the type of
            # the collection. We can replace the tuples by `Collection[_Node | Hashable]` if people complain.
            tuple[_Node, _Node]  # for normal graphs, this is (u, v)
            | tuple[_Node, _Node, Hashable]  # for multigraphs, this is (u, v, key)
        ]
        | None
    ) = None,
    oriented: bool = False,
    weight: str | None = None,
    *,
    dtype: DTypeLike | None = None,
) -> csc_array: ...
@_dispatchable
def adjacency_matrix(
    G: Graph[_Node], nodelist: Collection[_Node] | None = None, dtype: DTypeLike | None = None, weight: str | None = "weight"
) -> csr_array: ...
