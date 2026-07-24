from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["second_order_centrality"]

@_dispatchable
def second_order_centrality(G: Graph[_Node], weight: str | None = "weight") -> dict[Incomplete, float]: ...
