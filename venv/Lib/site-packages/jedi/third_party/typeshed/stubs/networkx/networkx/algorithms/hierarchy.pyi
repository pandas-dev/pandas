from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = ["flow_hierarchy"]

@_dispatchable
def flow_hierarchy(G: DiGraph[_Node], weight: str | None = None) -> float: ...
