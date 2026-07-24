from _typeshed import Incomplete

from networkx.classes.graph import Graph, _Node
from networkx.classes.multidigraph import MultiDiGraph
from networkx.utils.backends import _dispatchable

def write_p2g(G: Graph[_Node], path, encoding: str = "utf-8") -> None: ...
@_dispatchable
def read_p2g(path, encoding: str = "utf-8") -> MultiDiGraph[Incomplete]: ...
@_dispatchable
def parse_p2g(lines) -> MultiDiGraph[Incomplete]: ...
