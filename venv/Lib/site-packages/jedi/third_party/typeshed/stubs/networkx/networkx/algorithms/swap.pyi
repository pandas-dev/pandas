from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["double_edge_swap", "connected_double_edge_swap", "directed_edge_swap"]

@_dispatchable
def directed_edge_swap(G: DiGraph[_Node], *, nswap: int = 1, max_tries: int = 100, seed: int | RandomState | None = None): ...
@_dispatchable
def double_edge_swap(G: Graph[_Node], nswap: int = 1, max_tries: int = 100, seed: int | RandomState | None = None): ...
@_dispatchable
def connected_double_edge_swap(
    G: Graph[_Node], nswap: int = 1, _window_threshold: int = 3, seed: int | RandomState | None = None
) -> int: ...
