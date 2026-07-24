from _typeshed import Incomplete
from collections.abc import Iterable
from typing import Literal

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["average_degree_connectivity"]

@_dispatchable
def average_degree_connectivity(
    G: Graph[_Node],
    source: Literal["in+out", "out", "in"] = "in+out",
    target: Literal["in+out", "out", "in"] = "in+out",
    nodes: Iterable[Incomplete] | None = None,
    weight: str | None = None,
) -> dict[Incomplete, int | float]: ...
