from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["average_neighbor_degree"]

@_dispatchable
def average_neighbor_degree(
    G: Graph[_Node],
    source: str | None = "out",
    target: str | None = "out",
    nodes: Iterable[Incomplete] | None = None,
    weight: str | None = None,
) -> dict[Incomplete, Incomplete]: ...
