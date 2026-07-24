from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["current_flow_betweenness_centrality_subset", "edge_current_flow_betweenness_centrality_subset"]

@_dispatchable
def current_flow_betweenness_centrality_subset(
    G: Graph[_Node],
    sources: Iterable[_Node],
    targets: Iterable[_Node],
    normalized: bool | None = True,
    weight: str | None = None,
    dtype: type = ...,
    solver: str = "lu",
) -> dict[Incomplete, float]: ...
@_dispatchable
def edge_current_flow_betweenness_centrality_subset(
    G: Graph[_Node],
    sources: Iterable[_Node],
    targets: Iterable[_Node],
    normalized: bool | None = True,
    weight: str | None = None,
    dtype: type = ...,
    solver: str = "lu",
) -> dict[tuple[Incomplete, Incomplete], float]: ...
