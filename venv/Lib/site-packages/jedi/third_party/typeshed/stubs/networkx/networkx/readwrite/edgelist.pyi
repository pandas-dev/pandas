from _typeshed import Incomplete
from collections.abc import Generator

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable

__all__ = [
    "generate_edgelist",
    "write_edgelist",
    "parse_edgelist",
    "read_edgelist",
    "read_weighted_edgelist",
    "write_weighted_edgelist",
]

def generate_edgelist(G: Graph[_Node], delimiter: str = " ", data: bool = True) -> Generator[Incomplete]: ...
def write_edgelist(
    G: Graph[_Node], path, comments: str = "#", delimiter: str = " ", data: bool = True, encoding: str = "utf-8"
) -> None: ...
@_dispatchable
def parse_edgelist(lines, comments: str = "#", delimiter=None, create_using=None, nodetype=None, data: bool = True): ...
@_dispatchable
def read_edgelist(
    path,
    comments: str = "#",
    delimiter=None,
    create_using=None,
    nodetype=None,
    data: bool = True,
    edgetype=None,
    encoding: str = "utf-8",
): ...
def write_weighted_edgelist(
    G: Graph[_Node], path, comments: str = "#", delimiter: str = " ", encoding: str = "utf-8"
) -> None: ...
@_dispatchable
def read_weighted_edgelist(
    path, comments: str = "#", delimiter=None, create_using=None, nodetype=None, encoding: str = "utf-8"
): ...
