from _typeshed import Incomplete
from collections.abc import Callable, Generator

from networkx.algorithms.flow import edmonds_karp
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["all_node_cuts"]
default_flow_func = edmonds_karp

@_dispatchable
def all_node_cuts(
    G: Graph[_Node], k: int | None = None, flow_func: Callable[..., Incomplete] | None = None
) -> Generator[Incomplete]: ...
