from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["capacity_scaling"]

@_dispatchable
def capacity_scaling(
    G: Graph[_Node], demand: str = "demand", capacity: str = "capacity", weight: str = "weight", heap: type = ...
): ...
