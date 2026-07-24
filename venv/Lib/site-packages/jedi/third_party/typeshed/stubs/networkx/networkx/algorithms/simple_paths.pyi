from _typeshed import Incomplete, SupportsGetItem
from collections.abc import Callable, Collection, Generator, Iterable
from typing import Any

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["all_simple_paths", "is_simple_path", "shortest_simple_paths", "all_simple_edge_paths"]

@_dispatchable
def is_simple_path(G: Graph[_Node], nodes: Collection[Incomplete]) -> bool: ...
@_dispatchable
def all_simple_paths(
    G: Graph[_Node], source: _Node, target: _Node | Iterable[_Node], cutoff: int | None = None
) -> Generator[list[_Node]]: ...
@_dispatchable
def all_simple_edge_paths(
    G: Graph[_Node], source: _Node, target: _Node | Iterable[_Node], cutoff: int | None = None
) -> Generator[list[_Node] | list[tuple[_Node, _Node]], None, list[_Node] | None]: ...
@_dispatchable
def shortest_simple_paths(
    G: Graph[_Node],
    source: _Node,
    target: _Node,
    weight: str | Callable[[Any, Any, SupportsGetItem[str, Any]], float | None] | None = None,
) -> Generator[list[_Node]]: ...

class PathBuffer:
    paths: Incomplete
    sortedpaths: Incomplete
    counter: Incomplete

    def __init__(self) -> None: ...
    def __len__(self): ...
    def push(self, cost, path) -> None: ...
    def pop(self): ...
