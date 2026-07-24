from collections.abc import Iterable

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["birank"]

@_dispatchable
def birank(
    G: Graph[_Node],
    nodes: Iterable[_Node],
    *,
    alpha: float | None = None,
    beta: float | None = None,
    top_personalization: dict[str, int] | None = None,
    bottom_personalization: dict[str, int] | None = None,
    max_iter: int = 100,
    tol: float = 1.0e-6,
    weight: str | None = "weight",
) -> dict[_Node, float]: ...
