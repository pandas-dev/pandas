from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_semiconnected"]

@_dispatchable
def is_semiconnected(G: DiGraph[_Node]) -> bool: ...
