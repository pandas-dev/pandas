from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["closeness_vitality"]

@_dispatchable
def closeness_vitality(
    G: Graph[_Node], node=None, weight: str | None = None, wiener_index: float | None = None
) -> float | dict[Incomplete, Incomplete]: ...
