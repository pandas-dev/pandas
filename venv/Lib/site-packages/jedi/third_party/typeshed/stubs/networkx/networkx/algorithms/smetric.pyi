from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["s_metric"]

@_dispatchable
def s_metric(G: Graph[_Node]) -> float: ...
