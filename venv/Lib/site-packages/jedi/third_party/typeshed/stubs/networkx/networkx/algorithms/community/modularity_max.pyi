from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["greedy_modularity_communities", "naive_greedy_modularity_communities"]

@_dispatchable
def greedy_modularity_communities(
    G: Graph[_Node], weight: str | None = None, resolution: float | None = 1, cutoff: int | None = 1, best_n: int | None = None
) -> list[set[Incomplete]] | list[frozenset[Incomplete]]: ...
@_dispatchable
def naive_greedy_modularity_communities(G: Graph[_Node], resolution: float = 1, weight: str | None = None): ...
