from _typeshed import Incomplete
from collections.abc import Callable
from typing_extensions import deprecated

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["could_be_isomorphic", "fast_could_be_isomorphic", "faster_could_be_isomorphic", "is_isomorphic"]

@_dispatchable
def could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node], *, properties: str = "dtc") -> bool: ...
@deprecated("`graph_could_be_isomorphic` is a deprecated alias for `could_be_isomorphic`. Use `could_be_isomorphic` instead.")
def graph_could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node]) -> bool: ...
@_dispatchable
def fast_could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node]) -> bool: ...
@deprecated(
    "`fast_graph_could_be_isomorphic` is a deprecated alias for `fast_could_be_isomorphic`. "
    "Use `fast_could_be_isomorphic` instead."
)
def fast_graph_could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node]) -> bool: ...
@_dispatchable
def faster_could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node]) -> bool: ...
@deprecated(
    "`faster_graph_could_be_isomorphic` is a deprecated alias for `faster_could_be_isomorphic`. "
    "Use `faster_could_be_isomorphic` instead."
)
def faster_graph_could_be_isomorphic(G1: Graph[_Node], G2: Graph[_Node]) -> bool: ...
@_dispatchable
def is_isomorphic(
    G1: Graph[_Node],
    G2: Graph[_Node],
    node_match: Callable[..., Incomplete] | None = None,
    edge_match: Callable[..., Incomplete] | None = None,
) -> bool: ...
