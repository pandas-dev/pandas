from _typeshed import Incomplete
from collections.abc import Mapping

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["percolation_centrality"]

@_dispatchable
def percolation_centrality(
    G: Graph[_Node],
    attribute: str | None = "percolation",
    states: Mapping[Incomplete, Incomplete] | None = None,
    weight: str | None = None,
) -> dict[Incomplete, float]: ...
