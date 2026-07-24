from _typeshed import Incomplete
from collections.abc import Callable

from networkx.algorithms.flow import edmonds_karp
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["k_components"]
default_flow_func = edmonds_karp

@_dispatchable
def k_components(G: Graph[_Node], flow_func: Callable[..., Incomplete] | None = None) -> dict[Incomplete, Incomplete]: ...
def build_k_number_dict(kcomps) -> dict[Incomplete, Incomplete]: ...
