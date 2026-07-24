from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["junction_tree"]

@_dispatchable
def junction_tree(G: Graph[_Node]) -> Graph[Incomplete]: ...
