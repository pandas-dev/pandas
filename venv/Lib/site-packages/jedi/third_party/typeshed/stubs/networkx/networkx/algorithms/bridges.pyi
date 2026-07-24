from collections.abc import Generator
from typing import overload

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["bridges", "has_bridges", "local_bridges"]

@_dispatchable
def bridges(G: Graph[_Node], root: _Node | None = None) -> Generator[_Node]: ...
@_dispatchable
def has_bridges(G: Graph[_Node], root: _Node | None = None) -> bool: ...
@overload
def local_bridges(G: Graph[_Node], with_span: bool = True, weight: str | None = None) -> Generator[tuple[_Node, _Node]]: ...
@overload
def local_bridges(G: Graph[_Node], with_span: bool = True, weight: str | None = None) -> Generator[tuple[_Node, _Node, int]]: ...
