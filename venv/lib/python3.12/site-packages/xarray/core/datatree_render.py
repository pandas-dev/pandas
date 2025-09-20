"""
String Tree Rendering. Copied from anytree.

Minor changes to `RenderDataTree` include accessing `children.values()`, and
type hints.

"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from math import ceil
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from xarray.core.datatree import DataTree


class Row(NamedTuple):
    pre: str
    fill: str
    node: DataTree | str


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
        assert len(cont) == len(vertical) == len(end), (
            f"'{vertical}', '{cont}' and '{end}' need to have equal length"
        )

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
        >>> root = DataTree.from_dict(
        ...     {
        ...         "/": None,
        ...         "/sub0": None,
        ...         "/sub0/sub0B": None,
        ...         "/sub0/sub0A": None,
        ...         "/sub1": None,
        ...     },
        ...     name="root",
        ... )
        >>> print(RenderDataTree(root))
        <xarray.DataTree 'root'>
        Group: /
        ├── Group: /sub0
        │   ├── Group: /sub0/sub0B
        │   └── Group: /sub0/sub0A
        └── Group: /sub1
        """
        super().__init__("\u2502   ", "\u251c\u2500\u2500 ", "\u2514\u2500\u2500 ")


class RenderDataTree:
    def __init__(
        self,
        node: DataTree,
        style=None,
        childiter: type = list,
        maxlevel: int | None = None,
        maxchildren: int | None = None,
    ):
        """
        Render tree starting at `node`.
        Keyword Args:
            style (AbstractStyle): Render Style.
            childiter: Child iterator. Note, due to the use of node.children.values(),
                Iterables that change the order of children  cannot be used
                (e.g., `reversed`).
            maxlevel: Limit rendering to this depth.
            maxchildren: Limit number of children at each node.
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
        >>> root = DataTree.from_dict(
        ...     {
        ...         "/": Dataset({"a": 0, "b": 1}),
        ...         "/sub0": Dataset({"c": 2, "d": 3}),
        ...         "/sub0/sub0B": Dataset({"e": 4}),
        ...         "/sub0/sub0A": Dataset({"f": 5, "g": 6}),
        ...         "/sub1": Dataset({"h": 7}),
        ...     },
        ...     name="root",
        ... )

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

        # `maxchildren` limits the number of children per node

        >>> print(RenderDataTree(root, maxchildren=1).by_attr("name"))
        root
        ├── sub0
        │   ├── sub0B
        │   ...
        ...

        """
        if style is None:
            style = ContStyle()
        if not isinstance(style, AbstractStyle):
            style = style()
        self.node = node
        self.style = style
        self.childiter = childiter
        self.maxlevel = maxlevel
        self.maxchildren = maxchildren

    def __iter__(self) -> Iterator[Row]:
        return self.__next(self.node, tuple())

    def __next(
        self,
        node: DataTree,
        continues: tuple[bool, ...],
        level: int = 0,
    ) -> Iterator[Row]:
        yield RenderDataTree.__item(node, continues, self.style)
        children = node.children.values()
        level += 1
        if children and (self.maxlevel is None or level < self.maxlevel):
            nchildren = len(children)
            children = self.childiter(children)
            for i, (child, is_last) in enumerate(_is_last(children)):
                if (
                    self.maxchildren is None
                    or i < ceil(self.maxchildren / 2)
                    or i >= ceil(nchildren - self.maxchildren / 2)
                ):
                    yield from self.__next(
                        child,
                        continues + (not is_last,),
                        level=level,
                    )
                if (
                    self.maxchildren is not None
                    and nchildren > self.maxchildren
                    and i == ceil(self.maxchildren / 2)
                ):
                    yield RenderDataTree.__item("...", continues, self.style)

    @staticmethod
    def __item(
        node: DataTree | str, continues: tuple[bool, ...], style: AbstractStyle
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
            f"style={self.style!r}",
            f"childiter={self.childiter!r}",
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
        >>> root = DataTree.from_dict(
        ...     {
        ...         "/sub0/sub0B": Dataset({"foo": 4, "bar": 109}),
        ...         "/sub0/sub0A": None,
        ...         "/sub1/sub1A": None,
        ...         "/sub1/sub1B": Dataset({"bar": 8}),
        ...         "/sub1/sub1C/sub1Ca": None,
        ...     },
        ...     name="root",
        ... )
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
                if isinstance(node, str):
                    yield f"{fill}{node}"
                    continue
                attr = (
                    attrname(node)
                    if callable(attrname)
                    else getattr(node, attrname, "")
                )
                if isinstance(attr, list | tuple):
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
