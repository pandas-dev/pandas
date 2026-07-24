from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["line_graph", "inverse_line_graph"]

@_dispatchable
def line_graph(G: Graph[_Node], create_using=None): ...
@_dispatchable
def inverse_line_graph(G: Graph[_Node]): ...
