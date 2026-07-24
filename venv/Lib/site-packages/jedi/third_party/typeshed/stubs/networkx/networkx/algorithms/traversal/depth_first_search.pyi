from collections.abc import Callable, Generator, Iterable, Iterator
from typing import Literal

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "dfs_edges",
    "dfs_tree",
    "dfs_predecessors",
    "dfs_successors",
    "dfs_preorder_nodes",
    "dfs_postorder_nodes",
    "dfs_labeled_edges",
]

@_dispatchable
def dfs_edges(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[tuple[_Node, _Node]]: ...
@_dispatchable
def dfs_tree(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> DiGraph[_Node]: ...
@_dispatchable
def dfs_predecessors(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> dict[_Node, _Node]: ...
@_dispatchable
def dfs_successors(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> dict[_Node, list[_Node]]: ...
@_dispatchable
def dfs_postorder_nodes(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[_Node]: ...
@_dispatchable
def dfs_preorder_nodes(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[_Node]: ...
@_dispatchable
def dfs_labeled_edges(
    G: Graph[_Node],
    source: _Node | None = None,
    depth_limit: int | None = None,
    *,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[tuple[_Node, _Node, Literal["forward", "nontree", "reverse", "reverse-depth_limit"]]]: ...
