"""Create links in the HDF5 file.

This module implements containers for soft and external links.  Hard
links doesn't need a container as such as they are the same as regular
nodes (groups or leaves).

Classes:

    SoftLink
    ExternalLink

Functions:

Misc variables:

"""

from __future__ import annotations

from typing import Any, Literal, NoReturn, TYPE_CHECKING
from pathlib import Path

import tables as tb

from . import linkextension
from .node import Node
from .utils import lazyattr
from .attributeset import AttributeSet

if TYPE_CHECKING:
    from .group import Group


def _g_get_link_class(
    parent_id: int, name: str
) -> Literal["ExternalLink", "HardLink", "SoftLink", "UnImplemented"]:
    """Guess the link class."""
    return linkextension._get_link_class(parent_id, name)


class Link(Node):
    """Abstract base class for all PyTables links.

    A link is a node that refers to another node.  The Link class inherits from
    Node class and the links that inherits from Link are SoftLink and
    ExternalLink.  There is not a HardLink subclass because hard links behave
    like a regular Group or Leaf.  Contrarily to other nodes, links cannot have
    HDF5 attributes.  This is an HDF5 library limitation that might be solved
    in future releases.

    See :ref:`LinksTutorial` for a small tutorial on how to work with links.

    .. rubric:: Link attributes

    .. attribute:: target

        The path string to the pointed node.

    """

    # Properties
    @lazyattr
    def _v_attrs(self) -> AttributeSet:
        """Attributes.

        A *NoAttrs* instance replacing the typical *AttributeSet* instance of
        other node objects.  The purpose of *NoAttrs* is to make clear that
        HDF5 attributes are not supported in link nodes.
        """

        class NoAttrs(AttributeSet):
            def __getattr__(self, name: str) -> NoReturn:
                raise KeyError(
                    "you cannot get attributes from this "
                    "`%s` instance" % self.__class__.__name__
                )

            def __setattr__(self, name: str, value: Any) -> NoReturn:
                raise KeyError(
                    "you cannot set attributes to this "
                    "`%s` instance" % self.__class__.__name__
                )

            def _g_close(self) -> None:
                pass

        return NoAttrs(self)

    def __init__(
        self,
        parentnode: Group,
        name: str,
        target: str | None = None,
        _log: bool = False,
    ) -> None:
        self._v_new = target is not None
        self.target = target
        """The path string to the pointed node."""

        super().__init__(parentnode, name, _log)

    # Public and tailored versions for copy, move, rename and remove methods
    def copy(
        self,
        newparent: Group | None = None,
        newname: str | None = None,
        overwrite: bool = False,
        createparents: bool = False,
    ) -> Link:
        """Copy this link and return the new one.

        See :meth:`Node._f_copy` for a complete explanation of the arguments.
        Please note that there is no recursive flag since links do not have
        child nodes.

        """
        newnode = self._f_copy(
            newparent=newparent,
            newname=newname,
            overwrite=overwrite,
            createparents=createparents,
        )
        # Insert references to a `newnode` via `newname`
        newnode._v_parent._g_refnode(newnode, newname, True)
        return newnode

    def move(
        self,
        newparent: Group | None = None,
        newname: str | None = None,
        overwrite: bool = False,
    ) -> None:
        """Move or rename this link.

        See :meth:`Node._f_move` for a complete explanation of the arguments.

        """
        return self._f_move(
            newparent=newparent, newname=newname, overwrite=overwrite
        )

    def remove(self) -> None:
        """Remove this link from the hierarchy."""
        return self._f_remove()

    def rename(
        self, newname: str | None = None, overwrite: bool = False
    ) -> None:
        """Rename this link in place.

        See :meth:`Node._f_rename` for a complete explanation of the arguments.

        """
        return self._f_rename(newname=newname, overwrite=overwrite)

    def __repr__(self):
        return str(self)


class SoftLink(linkextension.SoftLink, Link):
    """Represents a soft link (aka symbolic link).

    A soft link is a reference to another node in the *same* file hierarchy.
    Provided that the target node exists, its attributes and methods can be
    accessed directly from the softlink using the normal `.` syntax.

    Softlinks also have the following public methods/attributes:

        * `target`
        * `dereference()`
        * `copy()`
        * `move()`
        * `remove()`
        * `rename()`
        * `is_dangling()`

    Note that these will override any correspondingly named methods/attributes
    of the target node.

    For backwards compatibility, it is also possible to obtain the target node
    via the `__call__()` special method (this action is called *dereferencing*;
    see below)

    Examples
    --------
    ::

        >>> import numpy as np
        >>> f = tb.open_file('/tmp/test_softlink.h5', 'w')
        >>> a = f.create_array('/', 'A', np.arange(10))
        >>> link_a = f.create_soft_link('/', 'link_A', target='/A')

        # transparent read/write access to a softlinked node
        >>> link_a[0] = -1
        >>> link_a[:], link_a.dtype
        (array([-1,  1,  2,  3,  4,  5,  6,  7,  8,  9]), dtype('int64'))

        # dereferencing a softlink using the __call__() method
        >>> link_a() is a
        True

        # SoftLink.remove() overrides Array.remove()
        >>> link_a.remove()
        >>> print(link_a)
        <closed tables.link.SoftLink at ...>
        >>> a[:], a.dtype
        (array([-1,  1,  2,  3,  4,  5,  6,  7,  8,  9]), dtype('int64'))
        >>> f.close()

    """

    # Class identifier.
    _c_classid = "SOFTLINK"

    # attributes with these names/prefixes are treated as attributes of the
    # SoftLink rather than the target node
    _link_attrnames = (
        "target",
        "dereference",
        "is_dangling",
        "copy",
        "move",
        "remove",
        "rename",
        "__init__",
        "__str__",
        "__repr__",
        "__unicode__",
        "__class__",
        "__dict__",
    )
    _link_attrprefixes = ("_f_", "_c_", "_g_", "_v_")

    def __call__(self) -> Node | None:
        """Dereference `self.target` and return the object.

        Examples
        --------
        ::

            >>> f = tb.open_file('tables/tests/slink.h5')
            >>> f.root.arr2
            /arr2 (SoftLink) -> /arr
            >>> print(f.root.arr2())
            /arr (Array(2,)) ''
            >>> f.close()

        """
        return self.dereference()

    def dereference(self) -> Node | None:
        """Dereference a link."""
        if self._v_isopen:
            target = self.target
            # Check for relative pathnames
            if not self.target.startswith("/"):
                target = self._v_parent._g_join(self.target)
            return self._v_file._get_node(target)
        else:
            return None

    def __getattribute__(self, attrname: str) -> Any:

        # get attribute of the SoftLink itself
        if (
            attrname in SoftLink._link_attrnames
            or attrname[:3] in SoftLink._link_attrprefixes
        ):
            return object.__getattribute__(self, attrname)

        # get attribute of the target node
        elif not self._v_isopen:
            raise tb.ClosedNodeError("the node object is closed")
        elif self.is_dangling():
            return None
        else:
            target_node = self.dereference()
            try:
                # __getattribute__() fails to get children of Groups
                return target_node.__getattribute__(attrname)
            except AttributeError:
                # some node classes (e.g. Array) don't implement __getattr__()
                return target_node.__getattr__(attrname)

    def __setattr__(self, attrname: str, value: Any) -> None:

        # set attribute of the SoftLink itself
        if (
            attrname in SoftLink._link_attrnames
            or attrname[:3] in SoftLink._link_attrprefixes
        ):
            object.__setattr__(self, attrname, value)

        # set attribute of the target node
        elif not self._v_isopen:
            raise tb.ClosedNodeError("the node object is closed")
        elif self.is_dangling():
            raise ValueError("softlink target does not exist")
        else:
            self.dereference().__setattr__(attrname, value)

    def __getitem__(self, key: str) -> Any:
        """Getitem magic method.

        The __getitem__ must be defined in the SoftLink class in order
        for array indexing syntax to work.
        """
        if not self._v_isopen:
            raise tb.ClosedNodeError("the node object is closed")
        elif self.is_dangling():
            raise ValueError("softlink target does not exist")
        else:
            return self.dereference().__getitem__(key)

    def __setitem__(self, key: str, value: Any) -> None:
        """Setitem magic method.

        The __setitem__ method must be defined in the SoftLink class
        in order for array indexing syntax to work.
        """
        if not self._v_isopen:
            raise tb.ClosedNodeError("the node object is closed")
        elif self.is_dangling():
            raise ValueError("softlink target does not exist")
        else:
            self.dereference().__setitem__(key, value)

    def is_dangling(self) -> bool:
        """Return True if the link is dangling."""
        return not (self.dereference() in self._v_file)

    def __str__(self) -> str:
        """Return a short string representation of the link.

        Examples
        --------
        ::

            >>> f = tb.open_file('tables/tests/slink.h5')
            >>> f.root.arr2
            /arr2 (SoftLink) -> /arr
            >>> f.close()

        """
        target = str(self.target)
        # Check for relative pathnames
        if not self.target.startswith("/"):
            target = self._v_parent._g_join(self.target)
        closed = "" if self._v_isopen else "closed "
        dangling = "" if target in self._v_file else " (dangling)"
        return (
            f"{closed}{self._v_pathname} ({self.__class__.__name__}) -> "
            f"{self.target}{dangling}"
        )


class ExternalLink(linkextension.ExternalLink, Link):
    """Represents an external link.

    An external link is a reference to a node in *another* file.
    Getting access to the pointed node (this action is called
    *dereferencing*) is done via the :meth:`__call__` special method
    (see below).

    .. rubric:: ExternalLink attributes

    .. attribute:: extfile

        The external file handler, if the link has been dereferenced.
        In case the link has not been dereferenced yet, its value is
        None.

    """

    # Class identifier.
    _c_classid = "EXTERNALLINK"

    def __init__(
        self,
        parentnode: Group,
        name: str,
        target: str | None = None,
        _log: bool = False,
    ) -> None:
        self.extfile = None
        """The external file handler, if the link has been dereferenced.
        In case the link has not been dereferenced yet, its value is
        None."""
        super().__init__(parentnode, name, target, _log)

    def _get_filename_node(self) -> tuple[str, str]:
        """Return the external filename and nodepath from `self.target`."""
        # This is needed for avoiding the 'C:\\file.h5' filepath notation
        filename, target = self.target.split(":/")
        return filename, "/" + target

    def __call__(self, **kwargs) -> Node:
        """Dereference self.target and return the object.

        You can pass all the arguments supported by the :func:`open_file`
        function (except filename, of course) so as to open the referenced
        external file.

        Examples
        --------
        ::

            >>> f = tb.open_file('tables/tests/elink.h5')
            >>> f.root.pep.pep2
            /pep/pep2 (ExternalLink) -> elink2.h5:/pep
            >>> pep2 = f.root.pep.pep2(mode='r')  # open in 'r'ead mode
            >>> print(pep2)
            /pep (Group) ''
            >>> pep2._v_file.filename       # belongs to referenced file
            'tables/tests/elink2.h5'
            >>> f.close()

        """
        filename, target = self._get_filename_node()

        if not Path(filename).is_absolute():
            # Resolve the external link with respect to this
            # file's directory.  See #306.
            filename = str(Path(self._v_file.filename).parent / filename)

        if self.extfile is None or not self.extfile.isopen:
            self.extfile = tb.open_file(filename, **kwargs)
        else:
            # XXX: implement better consistency checks
            assert self.extfile.filename == filename
            assert self.extfile.mode == kwargs.get("mode", "r")

        return self.extfile._get_node(target)

    def umount(self) -> None:
        """Safely unmount self.extfile, if opened."""
        extfile = self.extfile
        # Close external file, if open
        if extfile is not None and extfile.isopen:
            extfile.close()
            self.extfile = None

    def _f_close(self) -> None:
        """Especific close for external links."""
        self.umount()
        super()._f_close()

    def __str__(self) -> str:
        """Return a short string representation of the link.

        Examples
        --------
        ::

            >>> f = tb.open_file('tables/tests/elink.h5')
            >>> f.root.pep.pep2
            /pep/pep2 (ExternalLink) -> elink2.h5:/pep
            >>> f.close()

        """
        return (
            f"{self._v_pathname} ({self.__class__.__name__}) -> "
            f"{self.target}"
        )
