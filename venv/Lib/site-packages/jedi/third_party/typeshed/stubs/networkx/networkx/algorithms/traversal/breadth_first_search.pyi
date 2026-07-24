from collections.abc import Callable, Generator, Iterable, Iterator
from typing import Final, Literal

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "bfs_edges",
    "bfs_tree",
    "bfs_predecessors",
    "bfs_successors",
    "descendants_at_distance",
    "bfs_layers",
    "bfs_labeled_edges",
    "generic_bfs_edges",
]

@_dispatchable
def generic_bfs_edges(
    G: Graph[_Node], source: _Node, neighbors: Callable[[_Node], Iterable[_Node]] | None = None, depth_limit: int | None = None
) -> Generator[tuple[_Node, _Node]]: ...
@_dispatchable
def bfs_edges(
    G: Graph[_Node],
    source: _Node,
    reverse: bool | None = False,
    depth_limit: int | None = None,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[tuple[_Node, _Node]]: ...
@_dispatchable
def bfs_tree(
    G: Graph[_Node],
    source: _Node,
    reverse: bool | None = False,
    depth_limit: int | None = None,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> DiGraph[_Node]: ...
@_dispatchable
def bfs_predecessors(
    G: Graph[_Node],
    source: _Node,
    depth_limit: int | None = None,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[tuple[_Node, _Node]]: ...
@_dispatchable
def bfs_successors(
    G: Graph[_Node],
    source: _Node,
    depth_limit: int | None = None,
    sort_neighbors: Callable[[Iterator[_Node]], Iterable[_Node]] | None = None,
) -> Generator[tuple[_Node, list[_Node]]]: ...
@_dispatchable
def bfs_layers(G: Graph[_Node], sources: _Node | Iterable[_Node]) -> Generator[list[_Node]]: ...

REVERSE_EDGE: Final = "reverse"
TREE_EDGE: Final = "tree"
FORWARD_EDGE: Final = "forward"
LEVEL_EDGE: Final = "level"

@_dispatchable
def bfs_labeled_edges(
    G: Graph[_Node], sources: _Node | Iterable[_Node]
) -> Generator[tuple[_Node, _Node, Literal["tree", "level", "forward", "reverse"]]]: ...
@_dispatchable
def descendants_at_distance(G: Graph[_Node], source: _Node, distance: int) -> set[_Node]: ...
