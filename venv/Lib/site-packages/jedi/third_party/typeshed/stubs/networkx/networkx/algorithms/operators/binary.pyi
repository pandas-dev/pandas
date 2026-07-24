from _typeshed import Incomplete
from collections.abc import Hashable, Iterable
from typing import TypeVar

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["union", "compose", "disjoint_union", "intersection", "difference", "symmetric_difference", "full_join"]

@_dispatchable
def disjoint_union(G: Graph[_Node], H: Graph[_Node]): ...
@_dispatchable
def intersection(G: Graph[_Node], H: Graph[_Node]): ...
@_dispatchable
def difference(G: Graph[_Node], H: Graph[_Node]): ...
@_dispatchable
def symmetric_difference(G: Graph[_Node], H: Graph[_Node]): ...

_X_co = TypeVar("_X_co", bound=Hashable, covariant=True)
_Y_co = TypeVar("_Y_co", bound=Hashable, covariant=True)

@_dispatchable
def compose(G: Graph[_X_co], H: Graph[_Y_co]) -> DiGraph[_X_co | _Y_co]: ...
@_dispatchable
def full_join(G: Graph[_Node], H, rename=(None, None)): ...
@_dispatchable
def union(G: Graph[_X_co], H: Graph[_Y_co], rename: Iterable[Incomplete] | None = ()) -> DiGraph[_X_co | _Y_co]: ...
