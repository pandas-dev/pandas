from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["average_clustering"]

@_dispatchable
def average_clustering(G: Graph[_Node], trials: int = 1000, seed: int | RandomState | None = None) -> float: ...
