from _typeshed import Incomplete
from collections import defaultdict

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["k_components"]

@_dispatchable
def k_components(G: Graph[_Node], min_density: float = 0.95) -> defaultdict[Incomplete, list[Incomplete]]: ...
