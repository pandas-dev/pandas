from _typeshed import Incomplete, Unused
from collections.abc import Callable, Generator
from typing import Final

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "greedy_color",
    "strategy_connected_sequential",
    "strategy_connected_sequential_bfs",
    "strategy_connected_sequential_dfs",
    "strategy_independent_set",
    "strategy_largest_first",
    "strategy_random_sequential",
    "strategy_saturation_largest_first",
    "strategy_smallest_last",
]

@_dispatchable
def strategy_largest_first(G: Graph[_Node], colors: Unused): ...
@_dispatchable
def strategy_random_sequential(G: Graph[_Node], colors: Unused, seed=None): ...
@_dispatchable
def strategy_smallest_last(G: Graph[_Node], colors: Unused): ...
@_dispatchable
def strategy_independent_set(G: Graph[_Node], colors: Unused) -> Generator[Incomplete, Incomplete]: ...
@_dispatchable
def strategy_connected_sequential_bfs(G: Graph[_Node], colors): ...
@_dispatchable
def strategy_connected_sequential_dfs(G: Graph[_Node], colors): ...
@_dispatchable
def strategy_connected_sequential(G: Graph[_Node], colors: Unused, traversal: str = "bfs") -> Generator[Incomplete]: ...
@_dispatchable
def strategy_saturation_largest_first(G: Graph[_Node], colors) -> Generator[Incomplete, None, Incomplete]: ...

STRATEGIES: Final[dict[str, Callable[..., Incomplete]]]

@_dispatchable
def greedy_color(G: Graph[_Node], strategy="largest_first", interchange: bool = False): ...
