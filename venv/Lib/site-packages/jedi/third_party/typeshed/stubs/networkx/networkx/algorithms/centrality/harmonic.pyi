from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["harmonic_centrality"]

@_dispatchable
def harmonic_centrality(
    G: Graph[_Node], nbunch: Iterable[Incomplete] | None = None, distance=None, sources: Iterable[Incomplete] | None = None
) -> dict[Incomplete, int]: ...
