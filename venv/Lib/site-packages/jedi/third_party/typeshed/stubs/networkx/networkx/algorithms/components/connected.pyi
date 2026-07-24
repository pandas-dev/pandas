from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["number_connected_components", "connected_components", "is_connected", "node_connected_component"]

@_dispatchable
def connected_components(G: Graph[_Node]) -> Generator[set[_Node]]: ...
@_dispatchable
def number_connected_components(G: Graph[_Node]) -> int: ...
@_dispatchable
def is_connected(G: Graph[_Node]) -> bool: ...
@_dispatchable
def node_connected_component(G: Graph[_Node], n: _Node) -> set[_Node]: ...
