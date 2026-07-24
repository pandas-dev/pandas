from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.algorithms.shortest_paths.weighted import _WeightFunc
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["kernighan_lin_bisection", "spectral_modularity_bipartition", "greedy_node_swap_bipartition"]

@_dispatchable
def kernighan_lin_bisection(
    G: Graph[_Node],
    partition: tuple[Iterable[Incomplete], Iterable[Incomplete]] | None = None,
    max_iter: int = 10,
    weight: str | _WeightFunc[_Node] = "weight",
    seed: int | RandomState | None = None,
) -> tuple[set[Incomplete], set[Incomplete]]: ...
def spectral_modularity_bipartition(G: Graph[_Node]) -> tuple[set[Incomplete], set[Incomplete]]: ...
def greedy_node_swap_bipartition(
    G: Graph[_Node], *, init_split: tuple[set[Incomplete], set[Incomplete]] | None = None, max_iter: int = 10
) -> tuple[set[Incomplete], set[Incomplete]]: ...
