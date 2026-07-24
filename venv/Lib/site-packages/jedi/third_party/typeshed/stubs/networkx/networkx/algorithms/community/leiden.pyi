from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["leiden_communities", "leiden_partitions"]

@_dispatchable
def leiden_communities(
    G: Graph[_Node],
    weight: str | None = "weight",
    resolution: float = 1,
    max_level: int | None = None,
    seed: int | RandomState | None = None,
): ...
@_dispatchable
def leiden_partitions(
    G: Graph[_Node], weight: str | None = "weight", resolution: float = 1, seed: int | RandomState | None = None
): ...
