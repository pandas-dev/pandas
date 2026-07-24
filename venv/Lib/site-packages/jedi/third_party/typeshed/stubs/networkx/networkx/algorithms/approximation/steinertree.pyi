from _typeshed import Incomplete
from collections.abc import Iterable
from typing_extensions import deprecated

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["metric_closure", "steiner_tree"]

@_dispatchable
@deprecated(
    "`metric_closure` is deprecated and will be removed in NetworkX 3.8. Use `networkx.all_pairs_shortest_path_length` instead."
)
def metric_closure(G: Graph[_Node], weight="weight"): ...
@_dispatchable
def steiner_tree(G: Graph[_Node], terminal_nodes: Iterable[Incomplete], weight: str = "weight", method: str | None = None): ...
