from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["min_weighted_dominating_set", "min_edge_dominating_set"]

@_dispatchable
def min_weighted_dominating_set(G: Graph[_Node], weight: str | None = None) -> set[Incomplete]: ...
@_dispatchable
def min_edge_dominating_set(G: Graph[_Node]) -> set[Incomplete]: ...
