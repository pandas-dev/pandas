from collections.abc import Iterable

from networkx.classes.graph import Graph, _Edge, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["betweenness_centrality_subset", "edge_betweenness_centrality_subset"]

@_dispatchable
def betweenness_centrality_subset(
    G: Graph[_Node],
    sources: Iterable[_Node],
    targets: Iterable[_Node],
    normalized: bool | None = False,
    weight: str | None = None,
) -> dict[_Node, float]: ...
@_dispatchable
def edge_betweenness_centrality_subset(
    G: Graph[_Node],
    sources: Iterable[_Node],
    targets: Iterable[_Node],
    normalized: bool | None = False,
    weight: str | None = None,
) -> dict[_Edge[_Node], float]: ...
