from _typeshed import Incomplete, SupportsGetItem
from collections.abc import Callable
from typing import Any

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["voronoi_cells"]

@_dispatchable
def voronoi_cells(
    G: Graph[_Node],
    center_nodes: set[Incomplete],
    weight: str | Callable[[Any, Any, SupportsGetItem[str, Any]], float | None] | None = "weight",
) -> dict[Incomplete, set[Incomplete]]: ...
