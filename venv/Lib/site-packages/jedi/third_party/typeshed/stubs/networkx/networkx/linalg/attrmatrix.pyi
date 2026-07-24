from _typeshed import Incomplete
from collections.abc import Collection
from typing import Any, Literal

from networkx._typing import Array2D
from networkx.classes.graph import Graph, _Node
from networkx.utils.backends import _dispatchable
from numpy.typing import DTypeLike
from scipy.sparse import lil_array  # type: ignore[import-untyped]  # pyright: ignore[reportMissingImports]

__all__ = ["attr_matrix", "attr_sparse_matrix"]

@_dispatchable
def attr_matrix(
    G: Graph[_Node],
    edge_attr: str | None = None,
    node_attr: str | None = None,  # runtime also accepts `Callable[[_Node], object]`, but it is not documented
    normalized: bool = False,  # runtime also accepts `Callable[[_Node, _Node], object]`, but it is not documented
    rc_order: Collection[_Node] | None = None,
    dtype: DTypeLike | None = None,
    order: Literal["C", "F"] | None = None,
    # TODO: overload on rc_order and node_attr
    # (rc_order:[node], node_attr:None) -> 2D-array
    # (rc_order:[any], node_attr:str) -> 2D-array
    # (rc_order:None, node_attr:None) -> (2D-array, list[node])
    # (rc_order:None, node_attr:str) -> (2D-array, list[any])
) -> Array2D[Incomplete] | tuple[Array2D[Incomplete], list[_Node] | list[Any]]: ...
@_dispatchable
def attr_sparse_matrix(
    G: Graph[_Node],
    edge_attr: str | None = None,
    node_attr: str | None = None,  # runtime also accepts `Callable[[_Node], object]`, but it is not documented
    normalized: bool = False,  # runtime also accepts `Callable[[_Node, _Node], object]`, but it is not documented
    rc_order: Collection[_Node] | None = None,
    dtype: DTypeLike | None = None,
    # TODO: overload on rc_order and node_attr
    # (rc_order:[node], node_attr:None) -> lil_array
    # (rc_order:[any], node_attr:str) -> lil_array
    # (rc_order:None, node_attr:None) -> (lil_array, list[node])
    # (rc_order:None, node_attr:str) -> (lil_array, list[any])
) -> lil_array | tuple[lil_array, list[_Node] | list[Any]]: ...
