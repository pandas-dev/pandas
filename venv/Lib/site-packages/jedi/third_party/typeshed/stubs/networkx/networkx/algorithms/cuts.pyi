from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "boundary_expansion",
    "conductance",
    "cut_size",
    "edge_expansion",
    "mixing_expansion",
    "node_expansion",
    "normalized_cut_size",
    "volume",
]

@_dispatchable
def cut_size(G: Graph[_Node], S: Iterable[_Node], T: Iterable[_Node] | None = None, weight: str | None = None): ...
@_dispatchable
def volume(G: Graph[_Node], S: Iterable[_Node], weight: str | None = None): ...
@_dispatchable
def normalized_cut_size(G: Graph[_Node], S: Iterable[_Node], T: Iterable[_Node] | None = None, weight: str | None = None): ...
@_dispatchable
def conductance(G: Graph[_Node], S: Iterable[_Node], T: Iterable[_Node] | None = None, weight: str | None = None): ...
@_dispatchable
def edge_expansion(G: Graph[_Node], S: Iterable[_Node], T: Iterable[_Node] | None = None, weight: str | None = None): ...
@_dispatchable
def mixing_expansion(G: Graph[_Node], S: Iterable[_Node], T: Iterable[_Node] | None = None, weight: str | None = None): ...
@_dispatchable
def node_expansion(G: Graph[_Node], S: Iterable[_Node]): ...
@_dispatchable
def boundary_expansion(G: Graph[_Node], S: Iterable[_Node]): ...
