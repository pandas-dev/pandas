"""Here is defined the UnImplemented class."""

import warnings

from . import hdf5extension
from .utils import SizeType
from .node import Node
from .leaf import Leaf


class UnImplemented(hdf5extension.UnImplemented, Leaf):
    """This class represents datasets not supported by PyTables in an HDF5
    file.

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
    _c_classid = 'UNIMPLEMENTED'

    def __init__(self, parentnode, name):
        """Create the `UnImplemented` instance."""

        # UnImplemented objects always come from opening an existing node
        # (they can not be created).
        self._v_new = False
        """Is this the first time the node has been created?"""
        self.nrows = SizeType(0)
        """The length of the first dimension of the data."""
        self.shape = (SizeType(0),)
        """The shape of the stored data."""
        self.byteorder = None
        """The endianness of data in memory ('big', 'little' or
        'irrelevant')."""

        super().__init__(parentnode, name)

    def _g_open(self):
        (self.shape, self.byteorder, object_id) = self._open_unimplemented()
        try:
            self.nrows = SizeType(self.shape[0])
        except IndexError:
            self.nrows = SizeType(0)
        return object_id

    def _g_copy(self, newparent, newname, recursive, _log=True, **kwargs):
        """Do nothing.

        This method does nothing, but a ``UserWarning`` is issued.
        Please note that this method *does not return a new node*, but
        ``None``.

        """

        warnings.warn(
            "UnImplemented node %r does not know how to copy itself; skipping"
            % (self._v_pathname,))
        return None  # Can you see it?

    def _f_copy(self, newparent=None, newname=None,
                overwrite=False, recursive=False, createparents=False,
                **kwargs):
        """Do nothing.

        This method does nothing, since `UnImplemented` nodes can not
        be copied.  However, a ``UserWarning`` is issued.  Please note
        that this method *does not return a new node*, but ``None``.

        """

        # This also does nothing but warn.
        self._g_copy(newparent, newname, recursive, **kwargs)
        return None  # Can you see it?

    def __repr__(self):
        return """{}
  NOTE: <The UnImplemented object represents a PyTables unimplemented
         dataset present in the '{}' HDF5 file.  If you want to see this
         kind of HDF5 dataset implemented in PyTables, please contact the
         developers.>
""".format(str(self), self._v_file.filename)


# Classes reported as H5G_UNKNOWN by HDF5
class Unknown(Node):
    """This class represents nodes reported as *unknown* by the underlying
    HDF5 library.

    This class does not have any public instance variables or methods, except
    those inherited from the Node class.

    """

    # Class identifier
    _c_classid = 'UNKNOWN'

    def __init__(self, parentnode, name):
        """Create the `Unknown` instance."""

        self._v_new = False
        super().__init__(parentnode, name)

    def _g_new(self, parentnode, name, init=False):
        pass

    def _g_open(self):
        return 0

    def _g_copy(self, newparent, newname, recursive, _log=True, **kwargs):
        # Silently avoid doing copies of unknown nodes
        return None

    def _g_delete(self, parent):
        pass

    def __str__(self):
        pathname = self._v_pathname
        classname = self.__class__.__name__
        return f"{pathname} ({classname})"

    def __repr__(self):
        return f"""{self!s}
  NOTE: <The Unknown object represents a node which is reported as
         unknown by the underlying HDF5 library, but that might be
         supported in more recent HDF5 versions.>
"""


# These are listed here for backward compatibility with PyTables 0.9.x indexes
class OldIndexArray(UnImplemented):
    _c_classid = 'IndexArray'
