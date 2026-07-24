from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "bidirectional_shortest_path",
    "single_source_shortest_path",
    "single_source_shortest_path_length",
    "single_target_shortest_path",
    "single_target_shortest_path_length",
    "all_pairs_shortest_path",
    "all_pairs_shortest_path_length",
    "predecessor",
]

@_dispatchable
def single_source_shortest_path_length(G: Graph[_Node], source: _Node, cutoff: int | None = None) -> dict[Incomplete, int]: ...
@_dispatchable
def single_target_shortest_path_length(G: Graph[_Node], target: _Node, cutoff: int | None = None): ...
@_dispatchable
def all_pairs_shortest_path_length(G: Graph[_Node], cutoff: int | None = None) -> Generator[Incomplete]: ...
@_dispatchable
def bidirectional_shortest_path(G: Graph[_Node], source: _Node, target: _Node) -> list[Incomplete]: ...
@_dispatchable
def single_source_shortest_path(
    G: Graph[_Node], source: _Node, cutoff: int | None = None
) -> dict[Incomplete, list[Incomplete]]: ...
@_dispatchable
def single_target_shortest_path(
    G: Graph[_Node], target: _Node, cutoff: int | None = None
) -> dict[Incomplete, list[Incomplete]]: ...
@_dispatchable
def all_pairs_shortest_path(
    G: Graph[_Node], cutoff: int | None = None
) -> Generator[tuple[Incomplete, dict[Incomplete, list[Incomplete]]]]: ...
@_dispatchable
def predecessor(
    G: Graph[_Node], source: _Node, target: _Node | None = None, cutoff: int | None = None, return_seen: bool | None = None
): ...
