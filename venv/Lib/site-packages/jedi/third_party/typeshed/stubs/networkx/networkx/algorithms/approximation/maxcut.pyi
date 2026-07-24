from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["randomized_partitioning", "one_exchange"]

@_dispatchable
def randomized_partitioning(
    G: Graph[_Node], seed: int | RandomState | None = None, p: float = 0.5, weight: str | None = None
): ...
@_dispatchable
def one_exchange(
    G: Graph[_Node], initial_cut: set[Incomplete] | None = None, seed: int | RandomState | None = None, weight: str | None = None
): ...
