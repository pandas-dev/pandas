from _typeshed import Incomplete

from networkx.utils.backends import _dispatchable

from ..classes.graph import Graph

__all__ = ["partial_duplication_graph", "duplication_divergence_graph"]

@_dispatchable
def partial_duplication_graph(N, n, p, q, seed=None): ...
@_dispatchable
def duplication_divergence_graph(n, p, seed=None) -> Graph[Incomplete]: ...
