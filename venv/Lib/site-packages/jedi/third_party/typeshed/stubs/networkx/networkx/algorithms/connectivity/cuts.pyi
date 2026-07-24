from _typeshed import Incomplete
from collections.abc import Callable
from typing_extensions import Never

from networkx.algorithms.flow import edmonds_karp
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["minimum_st_node_cut", "minimum_node_cut", "minimum_st_edge_cut", "minimum_edge_cut"]
default_flow_func = edmonds_karp

@_dispatchable
def minimum_st_edge_cut(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[..., Incomplete] | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
) -> set[tuple[Incomplete, Incomplete]]: ...
@_dispatchable
def minimum_st_node_cut(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[..., Incomplete] | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
) -> dict[Never, Never] | set[Incomplete]: ...
@_dispatchable
def minimum_node_cut(
    G: Graph[_Node], s: _Node | None = None, t: _Node | None = None, flow_func: Callable[..., Incomplete] | None = None
) -> dict[Never, Never] | set[Incomplete]: ...
@_dispatchable
def minimum_edge_cut(
    G: Graph[_Node], s: _Node | None = None, t: _Node | None = None, flow_func: Callable[..., Incomplete] | None = None
) -> set[Incomplete]: ...
