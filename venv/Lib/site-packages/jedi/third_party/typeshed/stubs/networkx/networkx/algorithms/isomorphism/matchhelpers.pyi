from _typeshed import Incomplete
from types import FunctionType

from networkx.utils.backends import _dispatchable

__all__ = [
    "categorical_node_match",
    "categorical_edge_match",
    "categorical_multiedge_match",
    "numerical_node_match",
    "numerical_edge_match",
    "numerical_multiedge_match",
    "generic_node_match",
    "generic_edge_match",
    "generic_multiedge_match",
]

def copyfunc(f, name=None) -> FunctionType: ...
def allclose(x, y, rtol: float = 1.0000000000000001e-05, atol=1e-08) -> bool: ...
@_dispatchable
def categorical_node_match(attr, default): ...

categorical_edge_match: Incomplete

@_dispatchable
def categorical_multiedge_match(attr, default): ...
@_dispatchable
def numerical_node_match(attr, default, rtol: float = 1e-05, atol: float = 1e-08): ...

numerical_edge_match: Incomplete

@_dispatchable
def numerical_multiedge_match(attr, default, rtol: float = 1e-05, atol: float = 1e-08): ...
@_dispatchable
def generic_node_match(attr, default, op): ...

generic_edge_match: Incomplete

@_dispatchable
def generic_multiedge_match(attr, default, op): ...
