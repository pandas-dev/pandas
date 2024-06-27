from __future__ import annotations

from collections.abc import Iterator
from typing import Callable

from xarray.core.treenode import Tree

"""These iterators are copied from anytree.iterators, with minor modifications."""


class LevelOrderIter(Iterator):
    """Iterate over tree applying level-order strategy starting at `node`.
    This is the iterator used by `DataTree` to traverse nodes.

    Parameters
    ----------
    node : Tree
        Node in a tree to begin iteration at.
    filter_ : Callable, optional
        Function called with every `node` as argument, `node` is returned if `True`.
        Default is to iterate through all ``node`` objects in the tree.
    stop : Callable, optional
        Function that will cause iteration to stop if ``stop`` returns ``True``
        for ``node``.
    maxlevel : int, optional
        Maximum level to descend in the node hierarchy.

    Examples
    --------
    >>> from xarray.core.datatree import DataTree
    >>> from xarray.core.iterators import LevelOrderIter
    >>> f = DataTree(name="f")
    >>> b = DataTree(name="b", parent=f)
    >>> a = DataTree(name="a", parent=b)
    >>> d = DataTree(name="d", parent=b)
    >>> c = DataTree(name="c", parent=d)
    >>> e = DataTree(name="e", parent=d)
    >>> g = DataTree(name="g", parent=f)
    >>> i = DataTree(name="i", parent=g)
    >>> h = DataTree(name="h", parent=i)
    >>> print(f)
    DataTree('f', parent=None)
    ├── DataTree('b')
    │   ├── DataTree('a')
    │   └── DataTree('d')
    │       ├── DataTree('c')
    │       └── DataTree('e')
    └── DataTree('g')
        └── DataTree('i')
            └── DataTree('h')
    >>> [node.name for node in LevelOrderIter(f)]
    ['f', 'b', 'g', 'a', 'd', 'i', 'c', 'e', 'h']
    >>> [node.name for node in LevelOrderIter(f, maxlevel=3)]
    ['f', 'b', 'g', 'a', 'd', 'i']
    >>> [
    ...     node.name
    ...     for node in LevelOrderIter(f, filter_=lambda n: n.name not in ("e", "g"))
    ... ]
    ['f', 'b', 'a', 'd', 'i', 'c', 'h']
    >>> [node.name for node in LevelOrderIter(f, stop=lambda n: n.name == "d")]
    ['f', 'b', 'g', 'a', 'i', 'h']
    """

    def __init__(
        self,
        node: Tree,
        filter_: Callable | None = None,
        stop: Callable | None = None,
        maxlevel: int | None = None,
    ):
        self.node = node
        self.filter_ = filter_
        self.stop = stop
        self.maxlevel = maxlevel
        self.__iter = None

    def __init(self):
        node = self.node
        maxlevel = self.maxlevel
        filter_ = self.filter_ or LevelOrderIter.__default_filter
        stop = self.stop or LevelOrderIter.__default_stop
        children = (
            []
            if LevelOrderIter._abort_at_level(1, maxlevel)
            else LevelOrderIter._get_children([node], stop)
        )
        return self._iter(children, filter_, stop, maxlevel)

    @staticmethod
    def __default_filter(node: Tree) -> bool:
        return True

    @staticmethod
    def __default_stop(node: Tree) -> bool:
        return False

    def __iter__(self) -> Iterator[Tree]:
        return self

    def __next__(self) -> Iterator[Tree]:
        if self.__iter is None:
            self.__iter = self.__init()
        item = next(self.__iter)  # type: ignore[call-overload]
        return item

    @staticmethod
    def _abort_at_level(level: int, maxlevel: int | None) -> bool:
        return maxlevel is not None and level > maxlevel

    @staticmethod
    def _get_children(children: list[Tree], stop: Callable) -> list[Tree]:
        return [child for child in children if not stop(child)]

    @staticmethod
    def _iter(
        children: list[Tree], filter_: Callable, stop: Callable, maxlevel: int | None
    ) -> Iterator[Tree]:
        level = 1
        while children:
            next_children = []
            for child in children:
                if filter_(child):
                    yield child
                next_children += LevelOrderIter._get_children(
                    list(child.children.values()), stop
                )
            children = next_children
            level += 1
            if LevelOrderIter._abort_at_level(level, maxlevel):
                break
