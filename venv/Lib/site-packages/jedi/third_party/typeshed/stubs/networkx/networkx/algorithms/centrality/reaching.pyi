from _typeshed import Incomplete
from collections.abc import Mapping

from networkx.classes.digraph import DiGraph
from networkx.classes.graph import _Node
from networkx.utils.backends import _dispatchable

__all__ = ["global_reaching_centrality", "local_reaching_centrality"]

@_dispatchable
def global_reaching_centrality(G: DiGraph[_Node], weight: str | None = None, normalized: bool | None = True) -> float: ...
@_dispatchable
def local_reaching_centrality(
    G: DiGraph[_Node],
    v: _Node,
    paths: Mapping[Incomplete, Incomplete] | None = None,
    weight: str | None = None,
    normalized: bool | None = True,
) -> float: ...
