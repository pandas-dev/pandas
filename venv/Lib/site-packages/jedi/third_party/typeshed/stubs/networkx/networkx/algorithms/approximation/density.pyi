from collections.abc import Callable, Hashable
from typing import Literal
from typing_extensions import TypeAlias

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

_Algorithm: TypeAlias = Literal["greedy++", "fista"]

__all__ = ["densest_subgraph"]

ALGORITHMS: dict[_Algorithm, Callable[[Graph[Hashable], int], tuple[float, set[int]]]]

@_dispatchable
def densest_subgraph(G: Graph[_Node], iterations: int = 1, *, method: _Algorithm = "fista") -> tuple[float, set[int]]: ...
