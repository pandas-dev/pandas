from _typeshed import ConvertibleToFloat, Incomplete, SupportsItemAccess
from collections.abc import Iterable, Mapping

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["katz_centrality", "katz_centrality_numpy"]

@_dispatchable
def katz_centrality(
    G: Graph[_Node],
    alpha: float | None = 0.1,
    beta: ConvertibleToFloat | Iterable[Incomplete] | None = 1.0,
    max_iter: int | None = 1000,
    tol: float | None = 1e-06,
    nstart: SupportsItemAccess[Incomplete, Incomplete] | None = None,
    normalized: bool | None = True,
    weight: str | None = None,
) -> dict[Incomplete, Incomplete]: ...
@_dispatchable
def katz_centrality_numpy(
    G: Graph[_Node],
    alpha: float = 0.1,
    beta: float | Mapping[Incomplete, Incomplete] | None = 1.0,
    normalized: bool = True,
    weight: str | None = None,
) -> dict[Incomplete, Incomplete]: ...
