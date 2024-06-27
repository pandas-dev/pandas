"""
String Tree Rendering. Copied from anytree.

Minor changes to `RenderDataTree` include accessing `children.values()`, and
type hints.

"""

from __future__ import annotations

from collections import namedtuple
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from xarray.core.datatree import DataTree

Row = namedtuple("Row", ("pre", "fill", "node"))


class AbstractStyle:
    def __init__(self, vertical: str, cont: str, end: str):
        """
        Tree Render Style.
        Args:
            vertical: Sign for vertical line.
            cont: Chars for a continued branch.
            end: Chars for the last branch.
        """
        super().__init__()
        self.vertical = vertical
        self.cont = cont
        self.end = end
        assert (
            len(cont) == len(vertical) == len(end)
        ), f"'{vertical}', '{cont}' and '{end}' need to have equal length"

    @property
    def empty(self) -> str:
        """Empty string as placeholder."""
        return " " * len(self.end)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"


class ContStyle(AbstractStyle):
    def __init__(self):
        """
        Continued style, without gaps.

        >>> from xarray.core.datatree import DataTree
        >>> from xarray.core.datatree_render import RenderDataTree
        >>> root = DataTree(name="root")
        >>> s0 = DataTree(name="sub0", parent=root)
        >>> s0b = DataTree(name="sub0B", parent=s0)
        >>> s0a = DataTree(name="sub0A", parent=s0)
        >>> s1 = DataTree(name="sub1", parent=root)
        >>> print(RenderDataTree(root))
        DataTree('root', parent=None)
        ├── DataTree('sub0')
        │   ├── DataTree('sub0B')
        │   └── DataTree('sub0A')
        └── DataTree('sub1')
        """
        super().__init__("\u2502   ", "\u251c\u2500\u2500 ", "\u2514\u2500\u2500 ")


class RenderDataTree:
    def __init__(
        self,
        node: DataTree,
        style=ContStyle(),
        childiter: type = list,
        maxlevel: int | None = None,
    ):
        """
        Render tree starting at `node`.
        Keyword Args:
            style (AbstractStyle): Render Style.
            childiter: Child iterator. Note, due to the use of node.children.values(),
                Iterables that change the order of children  cannot be used
                (e.g., `reversed`).
            maxlevel: Limit rendering to this depth.
        :any:`RenderDataTree` is an iterator, returning a tuple with 3 items:
        `pre`
            tree prefix.
        `fill`
            filling for multiline entries.
        `node`
            :any:`NodeMixin` object.
        It is up to the user to assemble these parts to a whole.

        Examples
        --------

        >>> from xarray import Dataset
        >>> from xarray.core.datatree import DataTree
        >>> from xarray.core.datatree_render import RenderDataTree
        >>> root = DataTree(name="root", data=Dataset({"a": 0, "b": 1}))
        >>> s0 = DataTree(name="sub0", parent=root, data=Dataset({"c": 2, "d": 3}))
        >>> s0b = DataTree(name="sub0B", parent=s0, data=Dataset({"e": 4}))
        >>> s0a = DataTree(name="sub0A", parent=s0, data=Dataset({"f": 5, "g": 6}))
        >>> s1 = DataTree(name="sub1", parent=root, data=Dataset({"h": 7}))

        # Simple one line:

        >>> for pre, _, node in RenderDataTree(root):
        ...     print(f"{pre}{node.name}")
        ...
        root
        ├── sub0
        │   ├── sub0B
        │   └── sub0A
        └── sub1

        # Multiline:

        >>> for pre, fill, node in RenderDataTree(root):
        ...     print(f"{pre}{node.name}")
        ...     for variable in node.variables:
        ...         print(f"{fill}{variable}")
        ...
        root
        a
        b
        ├── sub0
        │   c
        │   d
        │   ├── sub0B
        │   │   e
        │   └── sub0A
        │       f
        │       g
        └── sub1
            h

        :any:`by_attr` simplifies attribute rendering and supports multiline:
        >>> print(RenderDataTree(root).by_attr())
        root
        ├── sub0
        │   ├── sub0B
        │   └── sub0A
        └── sub1

        # `maxlevel` limits the depth of the tree:

        >>> print(RenderDataTree(root, maxlevel=2).by_attr("name"))
        root
        ├── sub0
        └── sub1
        """
        if not isinstance(style, AbstractStyle):
            style = style()
        self.node = node
        self.style = style
        self.childiter = childiter
        self.maxlevel = maxlevel

    def __iter__(self) -> Iterator[Row]:
        return self.__next(self.node, tuple())

    def __next(
        self, node: DataTree, continues: tuple[bool, ...], level: int = 0
    ) -> Iterator[Row]:
        yield RenderDataTree.__item(node, continues, self.style)
        children = node.children.values()
        level += 1
        if children and (self.maxlevel is None or level < self.maxlevel):
            children = self.childiter(children)
            for child, is_last in _is_last(children):
                yield from self.__next(child, continues + (not is_last,), level=level)

    @staticmethod
    def __item(
        node: DataTree, continues: tuple[bool, ...], style: AbstractStyle
    ) -> Row:
        if not continues:
            return Row("", "", node)
        else:
            items = [style.vertical if cont else style.empty for cont in continues]
            indent = "".join(items[:-1])
            branch = style.cont if continues[-1] else style.end
            pre = indent + branch
            fill = "".join(items)
            return Row(pre, fill, node)

    def __str__(self) -> str:
        return str(self.node)

    def __repr__(self) -> str:
        classname = self.__class__.__name__
        args = [
            repr(self.node),
            f"style={repr(self.style)}",
            f"childiter={repr(self.childiter)}",
        ]
        return f"{classname}({', '.join(args)})"

    def by_attr(self, attrname: str = "name") -> str:
        """
        Return rendered tree with node attribute `attrname`.

        Examples
        --------

        >>> from xarray import Dataset
        >>> from xarray.core.datatree import DataTree
        >>> from xarray.core.datatree_render import RenderDataTree
        >>> root = DataTree(name="root")
        >>> s0 = DataTree(name="sub0", parent=root)
        >>> s0b = DataTree(
        ...     name="sub0B", parent=s0, data=Dataset({"foo": 4, "bar": 109})
        ... )
        >>> s0a = DataTree(name="sub0A", parent=s0)
        >>> s1 = DataTree(name="sub1", parent=root)
        >>> s1a = DataTree(name="sub1A", parent=s1)
        >>> s1b = DataTree(name="sub1B", parent=s1, data=Dataset({"bar": 8}))
        >>> s1c = DataTree(name="sub1C", parent=s1)
        >>> s1ca = DataTree(name="sub1Ca", parent=s1c)
        >>> print(RenderDataTree(root).by_attr("name"))
        root
        ├── sub0
        │   ├── sub0B
        │   └── sub0A
        └── sub1
            ├── sub1A
            ├── sub1B
            └── sub1C
                └── sub1Ca
        """

        def get() -> Iterator[str]:
            for pre, fill, node in self:
                attr = (
                    attrname(node)
                    if callable(attrname)
                    else getattr(node, attrname, "")
                )
                if isinstance(attr, (list, tuple)):
                    lines = attr
                else:
                    lines = str(attr).split("\n")
                yield f"{pre}{lines[0]}"
                for line in lines[1:]:
                    yield f"{fill}{line}"

        return "\n".join(get())


def _is_last(iterable: Iterable) -> Iterator[tuple[DataTree, bool]]:
    iter_ = iter(iterable)
    try:
        nextitem = next(iter_)
    except StopIteration:
        pass
    else:
        item = nextitem
        while True:
            try:
                nextitem = next(iter_)
                yield item, False
            except StopIteration:
                yield nextitem, True
                break
            item = nextitem
