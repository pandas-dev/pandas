from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["biadjacency_matrix", "from_biadjacency_matrix"]

@_dispatchable
def biadjacency_matrix(
    G: Graph[_Node],
    row_order: Iterable[_Node],
    column_order: Iterable[Incomplete] | None = None,
    dtype=None,
    weight: str | None = "weight",
    format="csr",
): ...  # Return is a complex union of scipy classes depending on the format param
@_dispatchable
def from_biadjacency_matrix(
    A,
    create_using: Graph[_Node] | None = None,
    edge_attribute: str = "weight",
    *,
    row_order: Iterable[Incomplete] | None = None,
    column_order: Iterable[Incomplete] | None = None,
): ...
