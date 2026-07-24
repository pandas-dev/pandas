from _typeshed import Incomplete
from collections.abc import Hashable
from typing import TypeVar

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

_X = TypeVar("_X", bound=Hashable)
_Y = TypeVar("_Y", bound=Hashable)

__all__ = [
    "tensor_product",
    "cartesian_product",
    "lexicographic_product",
    "strong_product",
    "power",
    "rooted_product",
    "corona_product",
    "modular_product",
]

@_dispatchable
def tensor_product(G: Graph[_X], H: Graph[_Y]) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def cartesian_product(G: Graph[_X], H: Graph[_Y]) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def lexicographic_product(G: Graph[_X], H: Graph[_Y]) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def strong_product(G: Graph[_X], H: Graph[_Y]) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def power(G: Graph[_Node], k): ...
@_dispatchable
def rooted_product(G: Graph[_X], H: Graph[_Y], root: _Y) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def corona_product(G: Graph[_X], H: Graph[_Y]) -> Graph[tuple[_X, _Y]]: ...
@_dispatchable
def modular_product(G: Graph[_Node], H) -> Graph[Incomplete]: ...
