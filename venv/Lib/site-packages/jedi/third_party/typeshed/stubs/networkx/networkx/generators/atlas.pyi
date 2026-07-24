import sys
from _typeshed import Incomplete
from typing import Final

from networkx.utils.backends import _dispatchable

from ..classes.graph import Graph

if sys.version_info >= (3, 11):
    from importlib.resources.abc import Traversable
else:
    from importlib.abc import Traversable

__all__ = ["graph_atlas", "graph_atlas_g"]

NUM_GRAPHS: Final = 1253
ATLAS_FILE: Final[Traversable]

@_dispatchable
def graph_atlas(i) -> Graph[Incomplete]: ...
@_dispatchable
def graph_atlas_g() -> list[Graph[Incomplete]]: ...
