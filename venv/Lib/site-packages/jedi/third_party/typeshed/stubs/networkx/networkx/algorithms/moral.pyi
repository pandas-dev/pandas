from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["moral_graph"]

@_dispatchable
def moral_graph(G: Graph[_Node]): ...
