from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["is_eulerian", "eulerian_circuit", "eulerize", "is_semieulerian", "has_eulerian_path", "eulerian_path"]

@_dispatchable
def is_eulerian(G: Graph[_Node]) -> bool: ...
@_dispatchable
def is_semieulerian(G: Graph[_Node]) -> bool: ...
@_dispatchable
def eulerian_circuit(G: Graph[_Node], source: _Node | None = None, keys: bool = False) -> Generator[Incomplete, Incomplete]: ...
@_dispatchable
def has_eulerian_path(G: Graph[_Node], source: _Node | None = None) -> bool: ...
@_dispatchable
def eulerian_path(G: Graph[_Node], source=None, keys: bool = False) -> Generator[Incomplete, Incomplete]: ...
@_dispatchable
def eulerize(G: Graph[_Node]): ...
