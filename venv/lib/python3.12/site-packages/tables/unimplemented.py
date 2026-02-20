"""Here is defined the UnImplemented class."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from . import hdf5extension
from .leaf import Leaf
from .node import Node
from .utils import SizeType

if TYPE_CHECKING:
    from .group import Group


class UnImplemented(hdf5extension.UnImplemented, Leaf):
    """Class represents datasets not supported by PyTables in an HDF5 file.

    When reading a generic HDF5 file (i.e. one that has not been created with
    PyTables, but with some other HDF5 library based tool), chances are that
    the specific combination of datatypes or dataspaces in some dataset might
    not be supported by PyTables yet. In such a case, this dataset will be
    mapped into an UnImplemented instance and the user will still be able to
    access the complete object tree of the generic HDF5 file. The user will
    also be able to *read and write the attributes* of the dataset, *access
    some of its metadata*, and perform *certain hierarchy manipulation
    operations* like deleting or moving (but not copying) the node. Of course,
    the user will not be able to read the actual data on it.

    This is an elegant way to allow users to work with generic HDF5 files
    despite the fact that some of its datasets are not supported by
    PyTables. However, if you are really interested in having full access to an
    unimplemented dataset, please get in contact with the developer team.

    This class does not have any public instance variables or methods, except
    those inherited from the Leaf class (see :ref:`LeafClassDescr`).

    """

    # Class identifier.
    _c_classid = "UNIMPLEMENTED"

    def __init__(self, parentnode: Group, name: str) -> None:
        """Create the `UnImplemented` instance."""
        # UnImplemented objects always come from opening an existing node
        # (they can not be created).
        self._v_new = False
        """Is this the first time the node has been created?"""
        self.nrows = SizeType(0)
        """The length of the first dimension of the data."""
        self.shape = (SizeType(0),)
        """The shape of the stored data."""
        self.byteorder: str | None = None
        """The endianness of data in memory ('big', 'little' or
        'irrelevant')."""

        super().__init__(parentnode, name)

    def _g_open(self) -> int:
        (self.shape, self.byteorder, object_id) = self._open_unimplemented()
        try:
            self.nrows = SizeType(self.shape[0])
        except IndexError:
            self.nrows = SizeType(0)
        return object_id

    def _g_copy(
        self,
        newparent: Group,
        newname: str,
        recursive: bool,
        _log: bool = True,
        **kwargs,
    ) -> None:
        """Do nothing.

        This method does nothing, but a ``UserWarning`` is issued.
        Please note that this method *does not return a new node*, but
        ``None``.

        """
        warnings.warn(
            f"UnImplemented node {self._v_pathname!r} does not know how "
            f"to copy itself; skipping"
        )
        return None  # Can you see it?

    def _f_copy(
        self,
        newparent: Group | None = None,
        newname: str | None = None,
        overwrite: bool = False,
        recursive: bool = False,
        createparents: bool = False,
        **kwargs,
    ) -> None:
        """Do nothing.

        This method does nothing, since `UnImplemented` nodes can not
        be copied.  However, a ``UserWarning`` is issued.  Please note
        that this method *does not return a new node*, but ``None``.

        """
        # This also does nothing but warn.
        self._g_copy(newparent, newname, recursive, **kwargs)
        return None  # Can you see it?

    def __repr__(self) -> str:
        return f"""{str(self)}
  NOTE: <The UnImplemented object represents a PyTables unimplemented
         dataset present in the '{self._v_file.filename}' HDF5 file.
         If you want to see this kind of HDF5 dataset implemented in
         PyTables, please contact the developers.>
"""


# Classes reported as H5G_UNKNOWN by HDF5
class Unknown(Node):
    """Class representing nodes reported as *unknown* by the HDF5 library.

    This class does not have any public instance variables or methods, except
    those inherited from the Node class.

    """

    # Class identifier
    _c_classid = "UNKNOWN"

    def __init__(self, parentnode: Group, name: str) -> None:
        """Create the `Unknown` instance."""
        self._v_new = False
        super().__init__(parentnode, name)

    def _g_new(self, parentnode: Group, name: str, init: bool = False) -> None:
        pass

    def _g_open(self) -> int:
        return 0

    def _g_copy(
        self,
        newparent: Group,
        newname: str,
        recursive: bool,
        _log: bool = True,
        **kwargs,
    ) -> None:
        # Silently avoid doing copies of unknown nodes
        return None

    def _g_delete(self, parent: Group) -> None:
        pass

    def __str__(self) -> str:
        pathname = self._v_pathname
        classname = self.__class__.__name__
        return f"{pathname} ({classname})"

    def __repr__(self) -> str:
        return f"""{self!s}
  NOTE: <The Unknown object represents a node which is reported as
         unknown by the underlying HDF5 library, but that might be
         supported in more recent HDF5 versions.>
"""


# These are listed here for backward compatibility with PyTables 0.9.x indexes
class OldIndexArray(UnImplemented):
    """Old IndexArray.

    Provided for compatibility with PyTables 0.9.
    """

    _c_classid = "IndexArray"
