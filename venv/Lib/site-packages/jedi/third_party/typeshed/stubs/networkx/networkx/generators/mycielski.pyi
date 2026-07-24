from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["mycielskian", "mycielski_graph"]

@_dispatchable
def mycielskian(G: Graph[_Node], iterations: int = 1): ...
@_dispatchable
def mycielski_graph(n): ...
