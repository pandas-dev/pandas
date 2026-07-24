from _typeshed import Incomplete
from collections.abc import Callable, Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["clustering", "average_clustering", "latapy_clustering", "robins_alexander_clustering"]

def cc_dot(nu, nv) -> float: ...
def cc_max(nu, nv) -> float: ...
def cc_min(nu, nv) -> float: ...

modes: dict[str, Callable[[Incomplete, Incomplete], float]]

@_dispatchable
def latapy_clustering(
    G: Graph[_Node], nodes: Iterable[Incomplete] | None = None, mode: str = "dot"
) -> dict[Incomplete, Incomplete]: ...

clustering = latapy_clustering

@_dispatchable
def average_clustering(G: Graph[_Node], nodes: Iterable[Incomplete] | None = None, mode: str = "dot") -> float: ...
@_dispatchable
def robins_alexander_clustering(G: Graph[_Node]) -> float: ...
