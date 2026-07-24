from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_perfect_graph"]

@_dispatchable
def is_perfect_graph(G: Graph[_Node]) -> bool: ...
