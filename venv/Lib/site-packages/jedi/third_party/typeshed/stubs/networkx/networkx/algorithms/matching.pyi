from _typeshed import Incomplete
from collections.abc import Iterable, Mapping

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "is_matching",
    "is_maximal_matching",
    "is_perfect_matching",
    "max_weight_matching",
    "min_weight_matching",
    "maximal_matching",
]

@_dispatchable
def maximal_matching(G: Graph[_Node]) -> set[Incomplete]: ...
def matching_dict_to_set(matching: Mapping[Incomplete, Incomplete]) -> set[Incomplete]: ...
@_dispatchable
def is_matching(G: Graph[_Node], matching: dict[Incomplete, Incomplete] | Iterable[Iterable[Incomplete]]) -> bool: ...
@_dispatchable
def is_maximal_matching(G: Graph[_Node], matching: dict[Incomplete, Incomplete] | Iterable[Iterable[Incomplete]]) -> bool: ...
@_dispatchable
def is_perfect_matching(G: Graph[_Node], matching: dict[Incomplete, Incomplete] | Iterable[Iterable[Incomplete]]) -> bool: ...
@_dispatchable
def min_weight_matching(G: Graph[_Node], weight: str | None = "weight") -> set[Incomplete]: ...
@_dispatchable
def max_weight_matching(
    G: Graph[_Node], maxcardinality: bool | None = False, weight: str | None = "weight"
) -> set[Incomplete]: ...
