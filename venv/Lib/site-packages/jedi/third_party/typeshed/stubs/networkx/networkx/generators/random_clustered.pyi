from collections.abc import Iterable
from typing import TypeVar, overload

from networkx import MultiGraph
from networkx.classes.graph import Graph
from networkx.utils.misc import _RandomState

_G = TypeVar("_G", bound=Graph[int])
__all__ = ["random_clustered_graph"]

@overload
def random_clustered_graph(
    joint_degree_sequence: Iterable[tuple[int, int]], create_using: None = None, seed: _RandomState = None
) -> MultiGraph[int]: ...
@overload
def random_clustered_graph(
    joint_degree_sequence: Iterable[tuple[int, int]], create_using: type[_G], seed: _RandomState = None
) -> _G: ...
