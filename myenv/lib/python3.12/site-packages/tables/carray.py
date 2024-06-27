"""Here is defined the CArray class."""

import sys

import numpy as np

from .atom import Atom
from .array import Array
from .utils import correct_byteorder, SizeType


# default version for CARRAY objects
# obversion = "1.0"    # Support for time & enumerated datatypes.
obversion = "1.1"    # Numeric and numarray flavors are gone.


class CArray(Array):
    """This class represents homogeneous datasets in an HDF5 file.

    The difference between a CArray and a normal Array (see
    :ref:`ArrayClassDescr`), from which it inherits, is that a CArray
    has a chunked layout and, as a consequence, it supports compression.
    You can use datasets of this class to easily save or load arrays to
    or from disk, with compression support included.

    CArray includes all the instance variables and methods of Array.
    Only those with different behavior are mentioned here.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.

        .. versionchanged:: 3.0
           Renamed from *parentNode* to *parentnode*.

    name : str
        The name of this node in its parent group.
    atom
       An `Atom` instance representing the *type* and *shape* of
       the atomic objects to be saved.

    shape
       The shape of the new array.

    title
       A description for this node (it sets the ``TITLE`` HDF5
       attribute on disk).

    filters
       An instance of the `Filters` class that provides
       information about the desired I/O filters to be applied
       during the life of this object.

    chunkshape
       The shape of the data chunk to be read or written in a
       single HDF5 I/O operation.  Filters are applied to those
       chunks of data.  The dimensionality of `chunkshape` must
       be the same as that of `shape`.  If ``None``, a sensible
       value is calculated (which is recommended).

    byteorder
        The byteorder of the data *on disk*, specified as 'little'
        or 'big'.  If this is not specified, the byteorder is that
        of the platform.

    track_times
        Whether time data associated with the leaf are recorded (object
        access time, raw data modification time, metadata change time, object
        birth time); default True.  Semantics of these times depend on their
        implementation in the HDF5 library: refer to documentation of the
        H5O_info_t data structure.  As of HDF5 1.8.15, only ctime (metadata
        change time) is implemented.

        .. versionadded:: 3.4.3

    Examples
    --------

    See below a small example of the use of the `CArray` class.
    The code is available in ``examples/carray1.py``::

        import numpy as np
        import tables as tb

        fileName = 'carray1.h5'
        shape = (200, 300)
        atom = tb.UInt8Atom()
        filters = tb.Filters(complevel=5, complib='zlib')

        h5f = tb.open_file(fileName, 'w')
        ca = h5f.create_carray(h5f.root, 'carray', atom, shape,
                               filters=filters)

        # Fill a hyperslab in ``ca``.
        ca[10:60, 20:70] = np.ones((50, 50))
        h5f.close()

        # Re-open a read another hyperslab
        h5f = tb.open_file(fileName)
        print(h5f)
        print(h5f.root.carray[8:12, 18:22])
        h5f.close()

    The output for the previous script is something like::

        carray1.h5 (File) ''
        Last modif.: 'Thu Apr 12 10:15:38 2007'
        Object Tree:
        / (RootGroup) ''
        /carray (CArray(200, 300), shuffle, zlib(5)) ''

        [[0 0 0 0]
         [0 0 0 0]
         [0 0 1 1]
         [0 0 1 1]]

    """

    # Class identifier.
    _c_classid = 'CARRAY'

    def __init__(self, parentnode, name,
                 atom=None, shape=None,
                 title="", filters=None,
                 chunkshape=None, byteorder=None,
                 _log=True, track_times=True):

        self.atom = atom
        """An `Atom` instance representing the shape, type of the atomic
        objects to be saved.
        """
        self.shape = None
        """The shape of the stored array."""
        self.extdim = -1  # `CArray` objects are not enlargeable by default
        """The index of the enlargeable dimension."""

        # Other private attributes
        self._v_version = None
        """The object version of this array."""
        self._v_new = new = atom is not None
        """Is this the first time the node has been created?"""
        self._v_new_title = title
        """New title for this node."""
        self._v_convert = True
        """Whether the ``Array`` object must be converted or not."""
        self._v_chunkshape = chunkshape
        """Private storage for the `chunkshape` property of the leaf."""

        # Miscellaneous iteration rubbish.
        self._start = None
        """Starting row for the current iteration."""
        self._stop = None
        """Stopping row for the current iteration."""
        self._step = None
        """Step size for the current iteration."""
        self._nrowsread = None
        """Number of rows read up to the current state of iteration."""
        self._startb = None
        """Starting row for current buffer."""
        self._stopb = None
        """Stopping row for current buffer. """
        self._row = None
        """Current row in iterators (sentinel)."""
        self._init = False
        """Whether we are in the middle of an iteration or not (sentinel)."""
        self.listarr = None
        """Current buffer in iterators."""

        if new:
            if not isinstance(atom, Atom):
                raise ValueError("atom parameter should be an instance of "
                                 "tables.Atom and you passed a %s." %
                                 type(atom))
            if shape is None:
                raise ValueError("you must specify a non-empty shape")
            try:
                shape = tuple(shape)
            except TypeError:
                raise TypeError("`shape` parameter must be a sequence "
                                "and you passed a %s" % type(shape))
            self.shape = tuple(SizeType(s) for s in shape)

            if chunkshape is not None:
                try:
                    chunkshape = tuple(chunkshape)
                except TypeError:
                    raise TypeError(
                        "`chunkshape` parameter must be a sequence "
                        "and you passed a %s" % type(chunkshape))
                if len(shape) != len(chunkshape):
                    raise ValueError(f"the shape ({shape}) and chunkshape "
                                     f"({chunkshape}) ranks must be equal.")
                elif min(chunkshape) < 1:
                    raise ValueError("chunkshape parameter cannot have "
                                     "zero-dimensions.")
                self._v_chunkshape = tuple(SizeType(s) for s in chunkshape)

        # The `Array` class is not abstract enough! :(
        super(Array, self).__init__(parentnode, name, new, filters,
                                    byteorder, _log, track_times)

    def _g_create(self):
        """Create a new array in file (specific part)."""

        if min(self.shape) < 1:
            raise ValueError(
                "shape parameter cannot have zero-dimensions.")
        # Finish the common part of creation process
        return self._g_create_common(self.nrows)

    def _g_create_common(self, expectedrows):
        """Create a new array in file (common part)."""

        self._v_version = obversion

        if self._v_chunkshape is None:
            # Compute the optimal chunk size
            self._v_chunkshape = self._calc_chunkshape(
                expectedrows, self.rowsize, self.atom.size)
        # Compute the optimal nrowsinbuf
        self.nrowsinbuf = self._calc_nrowsinbuf()
        # Correct the byteorder if needed
        if self.byteorder is None:
            self.byteorder = correct_byteorder(self.atom.type, sys.byteorder)

        try:
            # ``self._v_objectid`` needs to be set because would be
            # needed for setting attributes in some descendants later
            # on
            self._v_objectid = self._create_carray(self._v_new_title)
        except Exception:  # XXX
            # Problems creating the Array on disk. Close node and re-raise.
            self.close(flush=0)
            raise

        return self._v_objectid

    def _g_copy_with_stats(self, group, name, start, stop, step,
                           title, filters, chunkshape, _log, **kwargs):
        """Private part of Leaf.copy() for each kind of leaf."""

        (start, stop, step) = self._process_range_read(start, stop, step)
        maindim = self.maindim
        shape = list(self.shape)
        shape[maindim] = len(range(start, stop, step))
        # Now, fill the new carray with values from source
        nrowsinbuf = self.nrowsinbuf
        # The slices parameter for self.__getitem__
        slices = [slice(0, dim, 1) for dim in self.shape]
        # This is a hack to prevent doing unnecessary conversions
        # when copying buffers
        self._v_convert = False
        # Build the new CArray object
        object = CArray(group, name, atom=self.atom, shape=shape,
                        title=title, filters=filters, chunkshape=chunkshape,
                        _log=_log)
        # Start the copy itself
        for start2 in range(start, stop, step * nrowsinbuf):
            # Save the records on disk
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            # Set the proper slice in the main dimension
            slices[maindim] = slice(start2, stop2, step)
            start3 = (start2 - start) // step
            stop3 = start3 + nrowsinbuf
            if stop3 > shape[maindim]:
                stop3 = shape[maindim]
            # The next line should be generalised if, in the future,
            # maindim is designed to be different from 0 in CArrays.
            # See ticket #199.
            object[start3:stop3] = self.__getitem__(tuple(slices))
        # Activate the conversion again (default)
        self._v_convert = True
        nbytes = np.prod(self.shape, dtype=SizeType) * self.atom.size

        return (object, nbytes)
