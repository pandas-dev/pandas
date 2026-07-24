from _typeshed import OpenBinaryModeUpdating, OpenTextModeReading, OpenTextModeWriting, SupportsWrite
from collections.abc import Callable
from typing import IO, Any, Protocol, TypeVar, type_check_only

from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from pygraphviz.agraph import AGraph  # type: ignore[import-not-found]  # pyright: ignore[reportMissingImports]

__all__ = ["from_agraph", "to_agraph", "write_dot", "read_dot", "graphviz_layout", "pygraphviz_layout", "view_pygraphviz"]

_ModeT_contra = TypeVar("_ModeT_contra", bound=str, contravariant=True)
_FileT_co = TypeVar("_FileT_co", covariant=True)

@type_check_only
class _SupportsOpen(Protocol[_ModeT_contra, _FileT_co]):
    def open(self, *, mode: _ModeT_contra) -> _FileT_co: ...

@_dispatchable
def from_agraph(
    A: AGraph, create_using: Graph[str] | type[Graph[str]] | None = None  # TODO: add overloads on `create_using`
) -> Graph[str]: ...
def to_agraph(N: Graph[_Node]) -> AGraph: ...
def write_dot(
    G: Graph[_Node], path: str | IO[str] | IO[bytes] | _SupportsOpen[OpenTextModeWriting, IO[str] | IO[bytes]]
) -> None: ...
@_dispatchable
def read_dot(path: str | IO[str] | IO[bytes] | _SupportsOpen[OpenTextModeReading, IO[str] | IO[bytes]]) -> Graph[str]: ...
def graphviz_layout(
    G: Graph[_Node], prog: str = "neato", root: str | None = None, args: str = ""
) -> dict[_Node, tuple[float, float]]: ...
def pygraphviz_layout(
    G: Graph[_Node], prog: str = "neato", root: str | None = None, args: str = ""
) -> dict[_Node, tuple[float, float]]: ...
def view_pygraphviz(
    G: Graph[_Node],
    # From implementation looks like Callable could return object since it's always immediately stringified
    # But judging by documentation this seems like an extra runtime safety thing and not intended
    # Leaving as str unless anyone reports a valid use-case
    edgelabel: str | Callable[[dict[str, Any]], str] | None = None,
    prog: str = "dot",
    args: str = "",
    suffix: str = "",
    path: str | SupportsWrite[bytes] | _SupportsOpen[OpenBinaryModeUpdating, SupportsWrite[bytes]] | None = None,
    show: bool = True,
) -> tuple[str, AGraph]: ...
