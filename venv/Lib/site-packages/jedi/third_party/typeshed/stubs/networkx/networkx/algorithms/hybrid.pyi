from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["kl_connected_subgraph", "is_kl_connected"]

@_dispatchable
def kl_connected_subgraph(G: Graph[_Node], k: int, l: int, low_memory: bool = False, same_as_graph: bool = False): ...
@_dispatchable
def is_kl_connected(G: Graph[_Node], k: int, l: int, low_memory: bool = False) -> bool: ...
