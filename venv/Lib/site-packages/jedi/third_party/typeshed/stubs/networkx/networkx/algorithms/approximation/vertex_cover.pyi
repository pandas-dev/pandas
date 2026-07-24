from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["min_weighted_vertex_cover"]

@_dispatchable
def min_weighted_vertex_cover(G: Graph[_Node], weight: str | None = None) -> set[Incomplete]: ...
