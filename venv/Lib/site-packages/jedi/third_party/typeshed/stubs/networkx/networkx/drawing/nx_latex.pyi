from _typeshed import StrPath, SupportsWrite
from collections.abc import Collection
from typing_extensions import TypeAlias, TypeVar

from networkx.classes.graph import Graph, _Node

__all__ = ["to_latex_raw", "to_latex", "write_latex"]

# runtime requires a dict but it doesn't mutate it, we use a bounded typevar as
# a values type to make type checkers treat the dict covariantely
_PosT = TypeVar("_PosT", bound=Collection[float] | str)
_Pos: TypeAlias = str | dict[_Node, _PosT]

def to_latex_raw(
    G: Graph[_Node],
    pos: _Pos[_Node, _PosT] = "pos",
    tikz_options: str = "",
    default_node_options: str = "",
    node_options: str | dict[_Node, str] = "node_options",
    node_label: str | dict[_Node, str] = "label",
    default_edge_options: str = "",
    edge_options: str | dict[tuple[_Node, _Node], str] = "edge_options",
    edge_label: str | dict[tuple[_Node, _Node], str] = "label",
    edge_label_options: str | dict[tuple[_Node, _Node], str] = "edge_label_options",
) -> str: ...
def to_latex(
    Gbunch: Graph[_Node] | Collection[Graph[_Node]],
    pos: _Pos[_Node, _PosT] | Collection[_Pos[_Node, _PosT]] = "pos",
    tikz_options: str = "",
    default_node_options: str = "",
    node_options: str | dict[_Node, str] = "node_options",
    node_label: str | dict[_Node, str] = "node_label",
    default_edge_options: str = "",
    edge_options: str | dict[tuple[_Node, _Node], str] = "edge_options",
    edge_label: str | dict[tuple[_Node, _Node], str] = "edge_label",
    edge_label_options: str | dict[tuple[_Node, _Node], str] = "edge_label_options",
    caption: str = "",
    latex_label: str = "",
    sub_captions: Collection[str] | None = None,
    sub_labels: Collection[str] | None = None,
    n_rows: int = 1,
    as_document: bool = True,
    document_wrapper: str = ...,
    figure_wrapper: str = ...,
    subfigure_wrapper: str = ...,
) -> str: ...
def write_latex(
    Gbunch: Graph[_Node] | Collection[Graph[_Node]],
    path: StrPath | SupportsWrite[str],
    *,
    # **options passed to `to_latex`
    pos: _Pos[_Node, _PosT] | Collection[_Pos[_Node, _PosT]] = "pos",
    tikz_options: str = "",
    default_node_options: str = "",
    node_options: str | dict[_Node, str] = "node_options",
    node_label: str | dict[_Node, str] = "node_label",
    default_edge_options: str = "",
    edge_options: str | dict[tuple[_Node, _Node], str] = "edge_options",
    edge_label: str | dict[tuple[_Node, _Node], str] = "edge_label",
    edge_label_options: str | dict[tuple[_Node, _Node], str] = "edge_label_options",
    caption: str = "",
    latex_label: str = "",
    sub_captions: Collection[str] | None = None,
    sub_labels: Collection[str] | None = None,
    n_rows: int = 1,
    as_document: bool = True,
    document_wrapper: str = ...,
    figure_wrapper: str = ...,
    subfigure_wrapper: str = ...,
) -> None: ...
