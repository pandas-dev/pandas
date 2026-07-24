from _typeshed import Incomplete
from collections.abc import Callable, Mapping
from typing import Generic

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["treewidth_min_degree", "treewidth_min_fill_in"]

@_dispatchable
def treewidth_min_degree(G: Graph[_Node]) -> tuple[int, Graph[frozenset[_Node]]]: ...
@_dispatchable
def treewidth_min_fill_in(G: Graph[_Node]) -> tuple[int, Graph[frozenset[_Node]]]: ...

class MinDegreeHeuristic(Generic[_Node]):
    count: Incomplete

    def __init__(self, graph: Graph[_Node]) -> None: ...
    def best_node(self, graph: Mapping[_Node, set[_Node]]) -> _Node | None: ...

def min_fill_in_heuristic(graph_dict: Mapping[_Node, set[_Node]]) -> _Node | None: ...
@_dispatchable
def treewidth_decomp(
    G: Graph[_Node], heuristic: Callable[[dict[_Node, set[_Node]]], _Node | None] = ...
) -> tuple[int, Graph[frozenset[_Node]]]: ...
