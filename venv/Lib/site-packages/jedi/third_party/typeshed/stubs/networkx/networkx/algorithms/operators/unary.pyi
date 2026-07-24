from collections.abc import Hashable
from typing import TypeVar

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

_G = TypeVar("_G", bound=Graph[Hashable])

__all__ = ["complement", "reverse"]

@_dispatchable
def complement(G: Graph[_Node]): ...
@_dispatchable
def reverse(G: _G, copy: bool = True) -> _G: ...
