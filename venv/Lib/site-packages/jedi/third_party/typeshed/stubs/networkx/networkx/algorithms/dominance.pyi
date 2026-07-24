from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["immediate_dominators", "dominance_frontiers"]

@_dispatchable
def immediate_dominators(G: Graph[_Node], start: _Node): ...
@_dispatchable
def dominance_frontiers(G: Graph[_Node], start: _Node): ...
