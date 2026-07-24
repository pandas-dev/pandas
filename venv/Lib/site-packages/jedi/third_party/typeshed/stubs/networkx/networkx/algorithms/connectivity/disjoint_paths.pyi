from _typeshed import Incomplete
from collections.abc import Callable, Generator

from networkx.algorithms.flow import edmonds_karp
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["edge_disjoint_paths", "node_disjoint_paths"]
default_flow_func = edmonds_karp

@_dispatchable
def edge_disjoint_paths(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[..., Incomplete] | None = None,
    cutoff: int | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
) -> Generator[Incomplete]: ...
@_dispatchable
def node_disjoint_paths(
    G: Graph[_Node],
    s: _Node,
    t: _Node,
    flow_func: Callable[..., Incomplete] | None = None,
    cutoff: int | None = None,
    auxiliary: DiGraph[_Node] | None = None,
    residual: DiGraph[_Node] | None = None,
) -> Generator[Incomplete]: ...
