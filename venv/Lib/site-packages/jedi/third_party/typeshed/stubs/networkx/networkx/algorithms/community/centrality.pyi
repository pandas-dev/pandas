from _typeshed import Incomplete
from collections.abc import Callable, Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["girvan_newman"]

@_dispatchable
def girvan_newman(
    G: Graph[_Node], most_valuable_edge: Callable[..., Incomplete] | None = None
) -> Generator[Incomplete, None, Incomplete]: ...
