from collections.abc import Callable, Hashable
from typing import Literal
from typing_extensions import TypeAlias

from networkx.classes.graph import Graph, _Node

__all__ = ["greedy_source_expansion"]

_Algorithm: TypeAlias = Literal["clauset"]

ALGORITHMS: dict[_Algorithm, Callable[[Graph[Hashable], Hashable, int | None], set[Hashable]]]

def greedy_source_expansion(
    G: Graph[_Node], *, source: _Node, cutoff: int | None = None, method: _Algorithm = "clauset"
) -> set[_Node | None]: ...
