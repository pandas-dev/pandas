from collections.abc import Callable, Iterable

from networkx.algorithms.flow import edmonds_karp
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "average_node_connectivity",
    "local_node_connectivity",
    "node_connectivity",
    "local_edge_connectivity",
    "edge_connectivity",
    "all_pairs_node_connectivity",
]
default_flow_func = edmonds_karp

@_dispatchable
def local_node_connectivity(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
    cutoff: float | None = None,
) -> float: ...
@_dispatchable
def node_connectivity(
    G: Graph[_Node],
    s: _Node | None = None,
    t: _Node | None = None,
    flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None,
) -> float: ...
@_dispatchable
def average_node_connectivity(
    G: Graph[_Node], flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None
) -> float: ...
@_dispatchable
def all_pairs_node_connectivity(
    G: Graph[_Node],
    nbunch: Iterable[tuple[_Node, _Node]] | None = None,
    flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None,
) -> dict[_Node, dict[_Node, float]]: ...
@_dispatchable
def local_edge_connectivity(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
    cutoff: float | None = None,
) -> float: ...
@_dispatchable
def edge_connectivity(
    G: Graph[_Node],
    s: _Node | None = None,
    t: _Node | None = None,
    flow_func: Callable[[DiGraph[_Node], _Node, _Node], DiGraph[_Node]] | None = None,
    cutoff: float | None = None,
) -> float: ...
