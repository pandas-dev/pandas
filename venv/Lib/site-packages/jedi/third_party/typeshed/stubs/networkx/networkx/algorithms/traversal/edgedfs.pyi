from _typeshed import Incomplete
from collections.abc import Generator, Iterable
from typing import Final, Literal

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["edge_dfs"]

FORWARD: Final = "forward"
REVERSE: Final = "reverse"

@_dispatchable
def edge_dfs(
    G: Graph[_Node],
    source: _Node | Iterable[_Node] | None = None,
    orientation: Literal["original", "reverse", "ignore"] | None = None,
) -> Generator[tuple[Incomplete, ...]]: ...
