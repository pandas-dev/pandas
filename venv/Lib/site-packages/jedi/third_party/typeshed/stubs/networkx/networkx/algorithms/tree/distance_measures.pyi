from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["center", "centroid"]

def center(G: Graph[_Node]) -> list[Incomplete]: ...
@_dispatchable
def centroid(G: Graph[_Node]) -> list[Incomplete]: ...
