from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["shortest_augmenting_path"]

@_dispatchable
def shortest_augmenting_path(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    capacity: str = "capacity",
    residual: Graph[_Node] | None = None,
    value_only: bool = False,
    two_phase: bool = False,
    cutoff: float | None = None,
): ...
