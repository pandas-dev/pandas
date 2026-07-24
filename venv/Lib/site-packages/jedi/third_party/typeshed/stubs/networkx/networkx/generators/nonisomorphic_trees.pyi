from _typeshed import Incomplete
from collections.abc import Generator

from networkx.utils.backends import _dispatchable

__all__ = ["nonisomorphic_trees", "number_of_nonisomorphic_trees"]

@_dispatchable
def nonisomorphic_trees(order) -> Generator[list[Incomplete]]: ...
@_dispatchable
def number_of_nonisomorphic_trees(order): ...
