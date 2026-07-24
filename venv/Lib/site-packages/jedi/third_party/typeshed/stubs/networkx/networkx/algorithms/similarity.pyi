from _typeshed import Incomplete, SupportsItemAccess
from collections.abc import Callable, Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = [
    "graph_edit_distance",
    "optimal_edit_paths",
    "optimize_graph_edit_distance",
    "optimize_edit_paths",
    "simrank_similarity",
    "panther_similarity",
    "panther_vector_similarity",
    "generate_random_paths",
]

@_dispatchable
def graph_edit_distance(
    G1: Graph[_Node],
    G2: Graph[_Node],
    node_match: Callable[..., Incomplete] | None = None,
    edge_match: Callable[..., Incomplete] | None = None,
    node_subst_cost: Callable[..., Incomplete] | None = None,
    node_del_cost: Callable[..., Incomplete] | None = None,
    node_ins_cost: Callable[..., Incomplete] | None = None,
    edge_subst_cost: Callable[..., Incomplete] | None = None,
    edge_del_cost: Callable[..., Incomplete] | None = None,
    edge_ins_cost: Callable[..., Incomplete] | None = None,
    roots=None,
    upper_bound: float | None = None,
    timeout: float | None = None,
): ...
@_dispatchable
def optimal_edit_paths(
    G1: Graph[_Node],
    G2: Graph[_Node],
    node_match: Callable[..., Incomplete] | None = None,
    edge_match: Callable[..., Incomplete] | None = None,
    node_subst_cost: Callable[..., Incomplete] | None = None,
    node_del_cost: Callable[..., Incomplete] | None = None,
    node_ins_cost: Callable[..., Incomplete] | None = None,
    edge_subst_cost: Callable[..., Incomplete] | None = None,
    edge_del_cost: Callable[..., Incomplete] | None = None,
    edge_ins_cost: Callable[..., Incomplete] | None = None,
    upper_bound: float | None = None,
): ...
@_dispatchable
def optimize_graph_edit_distance(
    G1: Graph[_Node],
    G2: Graph[_Node],
    node_match: Callable[..., Incomplete] | None = None,
    edge_match: Callable[..., Incomplete] | None = None,
    node_subst_cost: Callable[..., Incomplete] | None = None,
    node_del_cost: Callable[..., Incomplete] | None = None,
    node_ins_cost: Callable[..., Incomplete] | None = None,
    edge_subst_cost: Callable[..., Incomplete] | None = None,
    edge_del_cost: Callable[..., Incomplete] | None = None,
    edge_ins_cost: Callable[..., Incomplete] | None = None,
    upper_bound: float | None = None,
) -> Generator[Incomplete]: ...
@_dispatchable
def optimize_edit_paths(
    G1: Graph[_Node],
    G2: Graph[_Node],
    node_match: Callable[..., Incomplete] | None = None,
    edge_match: Callable[..., Incomplete] | None = None,
    node_subst_cost: Callable[..., Incomplete] | None = None,
    node_del_cost: Callable[..., Incomplete] | None = None,
    node_ins_cost: Callable[..., Incomplete] | None = None,
    edge_subst_cost: Callable[..., Incomplete] | None = None,
    edge_del_cost: Callable[..., Incomplete] | None = None,
    edge_ins_cost: Callable[..., Incomplete] | None = None,
    upper_bound: float | None = None,
    strictly_decreasing: bool = True,
    roots=None,
    timeout: float | None = None,
) -> Generator[Incomplete, None, Incomplete]: ...
@_dispatchable
def simrank_similarity(
    G: Graph[_Node],
    source: _Node | None = None,
    target: _Node | None = None,
    importance_factor: float = 0.9,
    max_iterations: int = 1000,
    tolerance: float = 0.0001,
) -> float | dict[Incomplete, Incomplete]: ...
@_dispatchable
def panther_similarity(
    G: Graph[_Node],
    source: _Node,
    k: int = 5,
    path_length: int = 5,
    c: float = 0.5,
    delta: float = 0.1,
    eps: float | None = None,
    weight: str | None = "weight",
    seed: int | RandomState | None = None,
) -> dict[bytes, bytes]: ...
@_dispatchable
def panther_vector_similarity(
    G: Graph[_Node],
    source: _Node,
    *,
    D: int = 10,
    k: int = 5,
    path_length: int = 5,
    c: float = 0.5,
    delta: float = 0.1,
    eps: float | None = None,
    weight: str | None = "weight",
    seed: int | RandomState | None = None,
) -> dict[Incomplete, float]: ...
@_dispatchable
def generate_random_paths(
    G: Graph[_Node],
    sample_size: int,
    path_length: int = 5,
    index_map: SupportsItemAccess[Incomplete, Incomplete] | None = None,
    weight: str | None = "weight",
    seed: int | RandomState | None = None,
    *,
    source: _Node | None = None,
) -> Generator[list[Incomplete]]: ...
