import sys
from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.exception import NetworkXException
from networkx.utils.backends import _dispatchable

__all__ = [
    "is_chordal",
    "find_induced_nodes",
    "chordal_graph_cliques",
    "chordal_graph_treewidth",
    "NetworkXTreewidthBoundExceeded",
    "complete_to_chordal_graph",
]

class NetworkXTreewidthBoundExceeded(NetworkXException): ...

@_dispatchable
def is_chordal(G: Graph[_Node]) -> bool: ...
@_dispatchable
def find_induced_nodes(G: Graph[_Node], s: _Node, t: _Node, treewidth_bound: float = sys.maxsize) -> set[_Node]: ...
@_dispatchable
def chordal_graph_cliques(G: Graph[_Node]) -> Generator[frozenset[_Node]]: ...
@_dispatchable
def chordal_graph_treewidth(G: Graph[_Node]) -> int: ...
@_dispatchable
def complete_to_chordal_graph(G: Graph[_Node]) -> tuple[Incomplete, dict[Incomplete, int]]: ...
