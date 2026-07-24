from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["spectral_graph_forge"]

@_dispatchable
def spectral_graph_forge(G: Graph[_Node], alpha, transformation: str = "identity", seed=None) -> Graph[Incomplete]: ...
