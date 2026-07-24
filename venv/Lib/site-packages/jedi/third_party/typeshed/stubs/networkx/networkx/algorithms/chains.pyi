from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["chain_decomposition"]

@_dispatchable
def chain_decomposition(G: Graph[_Node], root: _Node | None = None) -> Generator[list[tuple[_Node, _Node]]]: ...
