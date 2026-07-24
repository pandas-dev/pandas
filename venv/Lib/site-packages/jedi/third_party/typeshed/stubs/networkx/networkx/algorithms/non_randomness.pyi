from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["non_randomness"]

@_dispatchable
def non_randomness(G: Graph[_Node], k: int | None = None, weight: str | None = "weight") -> tuple[float, float]: ...
