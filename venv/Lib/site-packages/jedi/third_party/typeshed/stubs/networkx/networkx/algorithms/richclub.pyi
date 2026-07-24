from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.random import RandomState

__all__ = ["rich_club_coefficient"]

@_dispatchable
def rich_club_coefficient(
    G: Graph[_Node], normalized: bool = True, Q: float = 100, seed: int | RandomState | None = None
) -> dict[Incomplete, Incomplete]: ...
