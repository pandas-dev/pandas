from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["tutte_polynomial", "chromatic_polynomial"]

@_dispatchable
def tutte_polynomial(G: Graph[_Node]): ...
@_dispatchable
def chromatic_polynomial(G: Graph[_Node]): ...
