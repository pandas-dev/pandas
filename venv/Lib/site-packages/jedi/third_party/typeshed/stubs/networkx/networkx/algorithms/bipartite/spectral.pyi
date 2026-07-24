from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["spectral_bipartivity"]

@_dispatchable
def spectral_bipartivity(G: Graph[_Node], nodes=None, weight: str = "weight") -> float | dict[Incomplete, Incomplete]: ...
