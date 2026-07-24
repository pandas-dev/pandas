from _typeshed import Incomplete
from typing import Final

from networkx import DiGraph
from networkx.utils.backends import _dispatchable

__all__ = ["triad_graph"]

TRIAD_EDGES: Final[dict[str, list[str]]]

@_dispatchable
def triad_graph(triad_name) -> DiGraph[Incomplete]: ...
