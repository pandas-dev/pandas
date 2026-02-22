"""Here is defined the EArray class."""

from __future__ import annotations

from typing import TYPE_CHECKING
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from .utils import convert_to_np_atom2, SizeType
from .carray import CArray

if TYPE_CHECKING:
    from .atom import Atom
    from .group import Group
    from .filters import Filters

# default version for EARRAY objects
# obversion = "1.0"    # initial version
# obversion = "1.1"    # support for complex datatypes
# obversion = "1.2"    # This adds support for time datatypes.
# obversion = "1.3"    # This adds support for enumerated datatypes.
obversion = "1.4"  # Numeric and numarray flavors are gone.


class EArray(CArray):
    r"""This class represents extendable, homogeneous datasets in an HDF5 file.

    The main difference between an EArray and a CArray (see
    :ref:`CArrayClassDescr`), from which it inherits, is that the former
    can be enlarged along one of its dimensions, the *enlargeable
    dimension*.  That means that the :attr:`Leaf.extdim` attribute (see
    :class:`Leaf`) of any EArray instance will always be non-negative.
    Multiple enlargeable dimensions might be supported in the future.

    New rows can be added to the end of an enlargeable array by using the
    :meth:`EArray.append` method.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.

        .. versionchanged:: 3.0
           Renamed from *parentNode* to *parentnode*.

    name : str
        The name of this node in its parent group.

    atom
        An `Atom` instance representing the *type* and *shape*
        of the atomic objects to be saved.

    shape
        The shape of the new array.  One (and only one) of
        the shape dimensions *must* be 0.  The dimension being 0
        means that the resulting `EArray` object can be extended
        along it.  Multiple enlargeable dimensions are not supported
        right now.

    title
        A description for this node (it sets the ``TITLE``
        HDF5 attribute on disk).

    filters
        An instance of the `Filters` class that provides information
        about the desired I/O filters to be applied during the life
        of this object.

    expectedrows
        A user estimate about the number of row elements that will
        be added to the growable dimension in the `EArray` node.
        If not provided, the default value is ``EXPECTED_ROWS_EARRAY``
        (see ``tables/parameters.py``).  If you plan to create either
        a much smaller or a much bigger `EArray` try providing a guess;
        this will optimize the HDF5 B-Tree creation and management
        process time and the amount of memory used.

    chunkshape
        The shape of the data chunk to be read or written in a single
        HDF5 I/O operation.  Filters are applied to those chunks of data.
        The dimensionality of `chunkshape` must be the same as that of
        `shape` (beware: no dimension should be 0 this time!).
        If ``None``, a sensible value is calculated based on the
        `expectedrows` parameter (which is recommended).

    byteorder
        The byteorder of the data *on disk*, specified as 'little' or
        'big'. If this is not specified, the byteorder is that of the
        platform.

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
    See below a small example of the use of the `EArray` class.  The
    code is available in ``examples/earray1.py``::

        import numpy as np
        import tables as tb

        fileh = tb.open_file('earray1.h5', mode='w')
        a = tb.StringAtom(itemsize=8)

        # Use ``a`` as the object type for the enlargeable array.
        array_c = fileh.create_earray(fileh.root, 'array_c', a, (0,),
                                      \"Chars\")
        array_c.append(np.array(['a'*2, 'b'*4], dtype='S8'))
        array_c.append(np.array(['a'*6, 'b'*8, 'c'*10], dtype='S8'))

        # Read the string ``EArray`` we have created on disk.
        for s in array_c:
            print('array_c[%s] => %r' % (array_c.nrow, s))
        # Close the file.
        fileh.close()

    The output for the previous script is something like::

        array_c[0] => 'aa'
        array_c[1] => 'bbbb'
        array_c[2] => 'aaaaaa'
        array_c[3] => 'bbbbbbbb'
        array_c[4] => 'cccccccc'

    """

    # Class identifier.
    _c_classid = "EARRAY"

    def __init__(
        self,
        parentnode: Group,
        name: str,
        atom: Atom | None = None,
        shape: Sequence[int] | None = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int | None = None,
        chunkshape: tuple[int, ...] | None = None,
        byteorder: str | None = None,
        _log: bool = True,
        track_times: bool = True,
    ) -> None:

        # Specific of EArray
        if expectedrows is None:
            expectedrows = parentnode._v_file.params["EXPECTED_ROWS_EARRAY"]
        self._v_expectedrows = int(expectedrows)
        """The expected number of rows to be stored in the array."""

        # Call the parent (CArray) init code
        super().__init__(
            parentnode,
            name,
            atom,
            shape,
            title,
            filters,
            chunkshape,
            byteorder,
            _log,
            track_times,
        )

    def _g_create(self) -> int:
        """Create a new array in file (specific part)."""
        # Pre-conditions and extdim computation
        zerodims = np.sum(np.array(self.shape) == 0)
        if zerodims > 0:
            if zerodims == 1:
                self.extdim = list(self.shape).index(0)
            else:
                raise NotImplementedError(
                    "Multiple enlargeable (0-)dimensions are not " "supported."
                )
        else:
            raise ValueError(
                "When creating EArrays, you need to set one of "
                "the dimensions of the Atom instance to zero."
            )

        # Finish the common part of the creation process
        return self._g_create_common(self._v_expectedrows)

    def _check_shape_append(self, nparr: np.ndarray) -> None:
        """Test that nparr shape is consistent with underlying EArray."""
        # Does the array conform to self expandibility?
        myrank = len(self.shape)
        narank = len(nparr.shape) - len(self.atom.shape)
        if myrank != narank:
            raise ValueError(
                (
                    "the ranks of the appended object (%d) and the "
                    "``%s`` EArray (%d) differ"
                )
                % (narank, self._v_pathname, myrank)
            )
        for i in range(myrank):
            if i != self.extdim and self.shape[i] != nparr.shape[i]:
                raise ValueError(
                    (
                        "the shapes of the appended object and the "
                        "``%s`` EArray differ in non-enlargeable "
                        "dimension %d"
                    )
                    % (self._v_pathname, i)
                )

    def append(self, sequence: npt.ArrayLike) -> None:
        """Add a sequence of data to the end of the dataset.

        The sequence must have the same type as the array; otherwise a
        TypeError is raised. In the same way, the dimensions of the
        sequence must conform to the shape of the array, that is, all
        dimensions must match, with the exception of the enlargeable
        dimension, which can be of any length (even 0!).  If the shape
        of the sequence is invalid, a ValueError is raised.

        """
        self._g_check_open()
        self._v_file._check_writable()

        # Convert the sequence into a NumPy object
        nparr = convert_to_np_atom2(sequence, self.atom)
        # Check if it has a consistent shape with underlying EArray
        self._check_shape_append(nparr)
        # If the size of the nparr is zero, don't do anything else
        if nparr.size > 0:
            self._append(nparr)

    def _g_copy_with_stats(
        self,
        group: Group,
        name: str,
        start: int,
        stop: int,
        step: int,
        title: str,
        filters: Filters | None,
        chunkshape: tuple[int, ...] | None,
        _log: bool,
        **kwargs,
    ) -> tuple[EArray, int]:
        """Private part of Leaf.copy() for each kind of leaf."""
        (start, stop, step) = self._process_range_read(start, stop, step)
        # Build the new EArray object
        maindim = self.maindim
        shape = list(self.shape)
        shape[maindim] = 0
        # The number of final rows
        nrows = len(range(start, stop, step))
        # Build the new EArray object
        obj = EArray(
            group,
            name,
            atom=self.atom,
            shape=shape,
            title=title,
            filters=filters,
            expectedrows=nrows,
            chunkshape=chunkshape,
            _log=_log,
        )
        # Now, fill the new earray with values from source
        nrowsinbuf = self.nrowsinbuf
        # The slices parameter for self.__getitem__
        slices = [slice(0, dim, 1) for dim in self.shape]
        # This is a hack to prevent doing unnecessary conversions
        # when copying buffers
        self._v_convert = False
        # Start the copy itself
        for start2 in range(start, stop, step * nrowsinbuf):
            # Save the records on disk
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            # Set the proper slice in the extensible dimension
            slices[maindim] = slice(start2, stop2, step)
            obj._append(self.__getitem__(tuple(slices)))
        # Active the conversion again (default)
        self._v_convert = True
        nbytes = np.prod(self.shape, dtype=SizeType) * self.atom.itemsize

        return (obj, nbytes)
