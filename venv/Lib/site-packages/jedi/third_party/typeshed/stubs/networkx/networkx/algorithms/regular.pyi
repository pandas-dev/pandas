from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_regular", "is_k_regular", "k_factor"]

@_dispatchable
def is_regular(G: Graph[_Node]) -> bool: ...
@_dispatchable
def is_k_regular(G: Graph[_Node], k) -> bool: ...
@_dispatchable
def k_factor(G: Graph[_Node], k, matching_weight: str | None = "weight"): ...
