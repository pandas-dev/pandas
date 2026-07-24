from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = ["generate_multiline_adjlist", "write_multiline_adjlist", "parse_multiline_adjlist", "read_multiline_adjlist"]

def generate_multiline_adjlist(G: Graph[_Node], delimiter: str = " ") -> Generator[str]: ...
def write_multiline_adjlist(G: Graph[_Node], path, delimiter=" ", comments="#", encoding="utf-8") -> None: ...
@_dispatchable
def parse_multiline_adjlist(lines, comments: str = "#", delimiter=None, create_using=None, nodetype=None, edgetype=None): ...
@_dispatchable
def read_multiline_adjlist(
    path, comments: str = "#", delimiter=None, create_using=None, nodetype=None, edgetype=None, encoding: str = "utf-8"
): ...
