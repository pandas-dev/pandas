from collections.abc import Callable, Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["bfs_beam_edges"]

@_dispatchable
def bfs_beam_edges(
    G: Graph[_Node], source: _Node, value: Callable[[_Node], float], width: int | None = None
) -> Generator[tuple[_Node, _Node]]: ...
