from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["min_maximal_matching"]

@_dispatchable
def min_maximal_matching(G: Graph[_Node]) -> set[Incomplete]: ...
