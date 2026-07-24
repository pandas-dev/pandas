from _typeshed import Incomplete
from collections.abc import Callable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

from .edmondskarp import edmonds_karp

__all__ = ["gomory_hu_tree"]
default_flow_func = edmonds_karp

@_dispatchable
def gomory_hu_tree(G: Graph[_Node], capacity: str = "capacity", flow_func: Callable[..., Incomplete] | None = None): ...
