from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["communicability", "communicability_exp"]

@_dispatchable
def communicability(G: Graph[_Node]) -> dict[dict[Incomplete, Incomplete], dict[Incomplete, float]]: ...
@_dispatchable
def communicability_exp(G: Graph[_Node]) -> dict[dict[Incomplete, Incomplete], dict[Incomplete, float]]: ...
