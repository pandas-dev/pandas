from _typeshed import Incomplete
from collections.abc import Iterable, Mapping

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["maximum_matching", "hopcroft_karp_matching", "eppstein_matching", "to_vertex_cover", "minimum_weight_full_matching"]

@_dispatchable
def hopcroft_karp_matching(G: Graph[_Node], top_nodes: Iterable[_Node] | None = None) -> dict[Incomplete, Incomplete]: ...
@_dispatchable
def eppstein_matching(G: Graph[_Node], top_nodes: Iterable[Incomplete] | None = None) -> dict[Incomplete, Incomplete]: ...
@_dispatchable
def to_vertex_cover(
    G: Graph[_Node], matching: Mapping[Incomplete, Incomplete], top_nodes: Iterable[Incomplete] | None = None
): ...

maximum_matching = hopcroft_karp_matching

@_dispatchable
def minimum_weight_full_matching(
    G: Graph[_Node], top_nodes: Iterable[Incomplete] | None = None, weight: str | None = "weight"
) -> dict[Incomplete, Incomplete]: ...
