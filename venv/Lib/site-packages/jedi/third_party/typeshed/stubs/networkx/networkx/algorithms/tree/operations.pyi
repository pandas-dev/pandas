from _typeshed import Incomplete
from collections.abc import Iterable

from networkx.utils.backends import _dispatchable

__all__ = ["join_trees"]

@_dispatchable
def join_trees(rooted_trees: Iterable[Incomplete], *, label_attribute: str | None = None, first_label: int | None = 0): ...
