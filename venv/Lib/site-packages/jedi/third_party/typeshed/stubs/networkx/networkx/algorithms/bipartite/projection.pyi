from _typeshed import Incomplete
from collections.abc import Callable, Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "projected_graph",
    "weighted_projected_graph",
    "collaboration_weighted_projected_graph",
    "overlap_weighted_projected_graph",
    "generic_weighted_projected_graph",
]

@_dispatchable
def projected_graph(B: Graph[_Node], nodes: Iterable[Incomplete], multigraph: bool = False): ...
@_dispatchable
def weighted_projected_graph(B: Graph[_Node], nodes: Iterable[Incomplete], ratio: bool = False): ...
@_dispatchable
def collaboration_weighted_projected_graph(B: Graph[_Node], nodes: Iterable[Incomplete]): ...
@_dispatchable
def overlap_weighted_projected_graph(B: Graph[_Node], nodes: Iterable[Incomplete], jaccard: bool = True): ...
@_dispatchable
def generic_weighted_projected_graph(
    B: Graph[_Node], nodes: Iterable[Incomplete], weight_function: Callable[..., Incomplete] | None = None
): ...
