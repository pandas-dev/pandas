from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["k_clique_communities"]

@_dispatchable
def k_clique_communities(G: Graph[_Node], k: int, cliques=None) -> Generator[Incomplete]: ...
