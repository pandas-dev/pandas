from _typeshed import Incomplete, SupportsLenAndGetItem
from collections.abc import Callable, Iterable, Mapping
from typing import Any, Literal, TypeVar

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "traveling_salesman_problem",
    "christofides",
    "asadpour_atsp",
    "greedy_tsp",
    "simulated_annealing_tsp",
    "threshold_accepting_tsp",
]

_SupportsLenAndGetItemT = TypeVar("_SupportsLenAndGetItemT", bound=SupportsLenAndGetItem[Any])

def swap_two_nodes(soln: _SupportsLenAndGetItemT, seed) -> _SupportsLenAndGetItemT: ...
def move_one_node(soln: _SupportsLenAndGetItemT, seed) -> _SupportsLenAndGetItemT: ...
@_dispatchable
def christofides(G: Graph[_Node], weight: str | None = "weight", tree: Graph[_Node] | None = None) -> list[Incomplete]: ...
@_dispatchable
def traveling_salesman_problem(
    G: Graph[_Node],
    weight: str = "weight",
    nodes=None,
    cycle: bool = True,
    method: Callable[..., Incomplete] | None = None,
    **kwargs,
): ...
@_dispatchable
def asadpour_atsp(
    G: DiGraph[_Node], weight: str | None = "weight", seed: int | RandomState | None = None, source: str | None = None
): ...
@_dispatchable
def held_karp_ascent(G: Graph[_Node], weight="weight"): ...
@_dispatchable
def spanning_tree_distribution(G: Graph[_Node], z: Mapping[Incomplete, Incomplete]) -> dict[Incomplete, Incomplete]: ...
@_dispatchable
def greedy_tsp(G: Graph[_Node], weight: str | None = "weight", source=None): ...
@_dispatchable
def simulated_annealing_tsp(
    G: Graph[_Node],
    init_cycle,
    weight: str | None = "weight",
    source=None,
    temp: int | None = 100,
    move="1-1",
    max_iterations: int | None = 10,
    N_inner: int | None = 100,
    alpha=0.01,
    seed: int | RandomState | None = None,
): ...
@_dispatchable
def threshold_accepting_tsp(
    G: Graph[_Node],
    init_cycle: Literal["greedy"] | Iterable[Incomplete],
    weight: str | None = "weight",
    source=None,
    threshold: int | None = 1,
    move="1-1",
    max_iterations: int | None = 10,
    N_inner: int | None = 100,
    alpha=0.1,
    seed: int | RandomState | None = None,
): ...
