from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["edge_betweenness_partition", "edge_current_flow_betweenness_partition"]

@_dispatchable
def edge_betweenness_partition(G: Graph[_Node], number_of_sets: int, *, weight: str | None = None) -> list[Incomplete]: ...
@_dispatchable
def edge_current_flow_betweenness_partition(
    G: Graph[_Node], number_of_sets: int, *, weight: str | None = None
) -> list[Incomplete]: ...
