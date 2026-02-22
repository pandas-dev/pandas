"""Here is defined the Array class."""

from __future__ import annotations

import sys
import operator
from typing import Any, Union, TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from . import hdf5extension
from .leaf import Leaf
from .utils import (
    is_idx,
    convert_to_np_atom2,
    SizeType,
    lazyattr,
    byteorders,
    quantize,
)
from .flavor import flavor_of, array_as_internal, internal_to_flavor
from .filters import Filters

if TYPE_CHECKING:
    from .atom import Atom, EnumAtom
    from .group import Group
    from .misc.enum import Enum

# default version for ARRAY objects
# obversion = "1.0"    # initial version
# obversion = "2.0"    # Added an optional EXTDIM attribute
# obversion = "2.1"    # Added support for complex datatypes
# obversion = "2.2"    # This adds support for time datatypes.
# obversion = "2.3"    # This adds support for enumerated datatypes.
obversion = "2.4"  # Numeric and numarray flavors are gone.

SelectionType = Union[int, slice, list[Union[int, slice]], npt.ArrayLike]


class Array(hdf5extension.Array, Leaf):
    """This class represents homogeneous datasets in an HDF5 file.

    This class provides methods to write or read data to or from array objects
    in the file. This class does not allow you neither to enlarge nor compress
    the datasets on disk; use the EArray class (see :ref:`EArrayClassDescr`) if
    you want enlargeable dataset support or compression features, or CArray
    (see :ref:`CArrayClassDescr`) if you just want compression.

    An interesting property of the Array class is that it remembers the
    *flavor* of the object that has been saved so that if you saved, for
    example, a list, you will get a list during readings afterwards; if you
    saved a NumPy array, you will get a NumPy object, and so forth.

    Note that this class inherits all the public attributes and methods that
    Leaf (see :ref:`LeafClassDescr`) already provides. However, as Array
    instances have no internal I/O buffers, it is not necessary to use the
    flush() method they inherit from Leaf in order to save their internal state
    to disk.  When a writing method call returns, all the data is already on
    disk.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.

        .. versionchanged:: 3.0
           Renamed from *parentNode* to *parentnode*

    name : str
        The name of this node in its parent group.
    obj
        The array or scalar to be saved.  Accepted types are NumPy
        arrays and scalars as well as native Python sequences and
        scalars, provided that values are regular (i.e. they are not
        like ``[[1,2],2]``) and homogeneous (i.e. all the elements are
        of the same type).

        .. versionchanged:: 3.0
           Renamed from *object* into *obj*.
    title
        A description for this node (it sets the ``TITLE`` HDF5 attribute on
        disk).
    byteorder
        The byteorder of the data *on disk*, specified as 'little' or 'big'.
        If this is not specified, the byteorder is that of the given `object`.
    track_times
        Whether time data associated with the leaf are recorded (object
        access time, raw data modification time, metadata change time, object
        birth time); default True.  Semantics of these times depend on their
        implementation in the HDF5 library: refer to documentation of the
        H5O_info_t data structure.  As of HDF5 1.8.15, only ctime (metadata
        change time) is implemented.

        .. versionadded:: 3.4.3

    """

    # Class identifier.
    _c_classid = "ARRAY"

    @lazyattr
    def dtype(self) -> np.dtype:
        """Numpy ``dtype`` most closely matching the one of the array."""
        return self.atom.dtype

    @property
    def nrows(self) -> int:
        """Return the number of rows in the array."""
        if self.shape == ():
            return SizeType(1)  # scalar case
        else:
            return self.shape[self.maindim]

    @property
    def rowsize(self) -> int:
        """Size of the rows in bytes in dimensions orthogonal to *maindim*."""
        maindim = self.maindim
        rowsize = self.atom.size
        for i, dim in enumerate(self.shape):
            if i != maindim:
                rowsize *= dim
        return rowsize

    @property
    def size_in_memory(self) -> int:
        """Size of the array's data in bytes when fully loaded in memory."""
        return self.nrows * self.rowsize

    def __init__(
        self,
        parentnode: Group,
        name: str,
        obj: npt.ArrayLike | None = None,
        title: str = "",
        byteorder: str | None = None,
        _log: bool = True,
        _atom: Atom | EnumAtom | None = None,
        track_times: bool = True,
    ) -> None:

        self._v_version: str | None = None
        """The object version of this array."""
        self._v_new = new = obj is not None
        """Is this the first time the node has been created?"""
        self._v_new_title = title
        """New title for this node."""
        self._obj = obj
        """The object to be stored in the array.  It can be any of numpy,
        list, tuple, string, integer of floating point types, provided
        that they are regular (i.e. they are not like ``[[1, 2], 2]``).

        .. versionchanged:: 3.0
           Renamed form *_object* into *_obj*.

        """

        self._v_convert = True
        """Whether the ``Array`` object must be converted or not."""

        # Miscellaneous iteration rubbish.
        self._start: int | None = None
        """Starting row for the current iteration."""
        self._stop: int | None = None
        """Stopping row for the current iteration."""
        self._step: int | None = None
        """Step size for the current iteration."""
        self._nrowsread: int | None = None
        """Number of rows read up to the current state of iteration."""
        self._startb: int | None = None
        """Starting row for current buffer."""
        self._stopb: int | None = None
        """Stopping row for current buffer. """
        self._row: int | None = None
        """Current row in iterators (sentinel)."""
        self._init = False
        """Whether we are in the middle of an iteration or not (sentinel)."""
        self.listarr: npt.ArrayLike | None = None
        """Current buffer in iterators."""

        # Documented (*public*) attributes.
        self.atom = _atom
        """An Atom (see :ref:`AtomClassDescr`) instance representing the *type*
        and *shape* of the atomic objects to be saved.
        """
        self.shape: list[int] | None = None
        """The shape of the stored array."""
        self.nrow: int | None = None
        """On iterators, this is the index of the current row."""
        self.extdim = -1  # ordinary arrays are not enlargeable
        """The index of the enlargeable dimension."""

        # Ordinary arrays have no filters: leaf is created with default ones.
        super().__init__(
            parentnode, name, new, Filters(), byteorder, _log, track_times
        )

    def _g_create(self) -> int:
        """Save a new array in file."""
        self._v_version = obversion
        try:
            # `Leaf._g_post_init_hook()` should be setting the flavor on disk.
            self._flavor = flavor = flavor_of(self._obj)
            nparr = array_as_internal(self._obj, flavor)
        except Exception:  # XXX
            # Problems converting data. Close the node and re-raise exception.
            self.close(flush=0)
            raise

        # Raise an error in case of unsupported object
        if nparr.dtype.kind in ["V", "U", "O"]:  # in void, unicode, object
            raise TypeError(
                "Array objects cannot currently deal with void, "
                "unicode or object arrays"
            )

        # Decrease the number of references to the object
        self._obj = None

        # Fix the byteorder of data
        nparr = self._g_fix_byteorder_data(nparr, nparr.dtype.byteorder)

        # Create the array on-disk
        try:
            # ``self._v_objectid`` needs to be set because would be
            # needed for setting attributes in some descendants later
            # on
            (self._v_objectid, self.shape, self.atom) = self._create_array(
                nparr, self._v_new_title, self.atom
            )
        except Exception:  # XXX
            # Problems creating the Array on disk. Close node and re-raise.
            self.close(flush=0)
            raise

        # Compute the optimal buffer size
        self.nrowsinbuf = self._calc_nrowsinbuf()
        # Arrays don't have chunkshapes (so, set it to None)
        self._v_chunkshape = None

        return self._v_objectid

    def _g_open(self) -> int:
        """Get the metadata info for an array in file."""
        (oid, self.atom, self.shape, self._v_chunkshape) = self._open_array()

        self.nrowsinbuf = self._calc_nrowsinbuf()

        return oid

    def get_enum(self) -> Enum:
        """Get the enumerated type associated with this array.

        If this array is of an enumerated type, the corresponding Enum instance
        (see :ref:`EnumClassDescr`) is returned. If it is not of an enumerated
        type, a TypeError is raised.

        """
        if self.atom.kind != "enum":
            raise TypeError(
                "array ``%s`` is not of an enumerated type" % self._v_pathname
            )

        return self.atom.enum

    def iterrows(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> tuple | Array:
        """Iterate over the rows of the array.

        This method returns an iterator yielding an object of the current
        flavor for each selected row in the array.  The returned rows are taken
        from the *main dimension*.

        If a range is not supplied, *all the rows* in the array are iterated
        upon - you can also use the :meth:`Array.__iter__` special method for
        that purpose.  If you only want to iterate over a given *range of rows*
        in the array, you may use the start, stop and step parameters.

        Examples
        --------
        ::

            result = [row for row in arrayInstance.iterrows(step=4)]

        .. versionchanged:: 3.0
           If the *start* parameter is provided and *stop* is None then the
           array is iterated from *start* to the last line.
           In PyTables < 3.0 only one element was returned.

        """
        try:
            (self._start, self._stop, self._step) = self._process_range(
                start, stop, step
            )
        except IndexError:
            # If problems with indexes, silently return the null tuple
            return ()
        self._init_loop()
        return self

    def __iter__(self) -> Array:
        """Iterate over the rows of the array.

        This is equivalent to calling :meth:`Array.iterrows` with default
        arguments, i.e. it iterates over *all the rows* in the array.

        Examples
        --------
        ::

            result = [row[2] for row in array]

        Which is equivalent to::

            result = [row[2] for row in array.iterrows()]

        """
        if not self._init:
            # If the iterator is called directly, assign default variables
            self._start = 0
            self._stop = self.nrows
            self._step = 1
            # and initialize the loop
            self._init_loop()
        return self

    def _init_loop(self) -> None:
        """Initialize the __iter__ iterator."""
        self._nrowsread = self._start
        self._startb = self._start
        self._row = -1  # Sentinel
        self._init = True  # Sentinel
        self.nrow = SizeType(self._start - self._step)  # row number

    def __next__(self) -> Any:
        """Get the next element of the array during an iteration.

        The element is returned as an object of the current flavor.

        """
        # this could probably be sped up for long iterations by reusing the
        # listarr buffer
        if self._nrowsread >= self._stop:
            self._init = False
            self.listarr = None  # fixes issue #308
            raise StopIteration  # end of iteration
        else:
            # Read a chunk of rows
            if self._row + 1 >= self.nrowsinbuf or self._row < 0:
                self._stopb = self._startb + self._step * self.nrowsinbuf
                # Protection for reading more elements than needed
                if self._stopb > self._stop:
                    self._stopb = self._stop
                listarr = self._read(self._startb, self._stopb, self._step)
                # Swap the axes to easy the return of elements
                if self.extdim > 0:
                    listarr = listarr.swapaxes(self.extdim, 0)
                self.listarr = internal_to_flavor(listarr, self.flavor)
                self._row = -1
                self._startb = self._stopb
            self._row += 1
            self.nrow += self._step
            self._nrowsread += self._step
            # Fixes bug #968132
            # if self.listarr.shape:
            if self.shape:
                return self.listarr[self._row]
            else:
                return self.listarr  # Scalar case

    def _interpret_indexing(
        self,
        keys: SelectionType,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
        """Implement the common part of `__getitem__` and `__setitem__`."""
        maxlen = len(self.shape)
        shape = (maxlen,)
        startl = np.empty(shape=shape, dtype=SizeType)
        stopl = np.empty(shape=shape, dtype=SizeType)
        stepl = np.empty(shape=shape, dtype=SizeType)
        stop_none = np.zeros(shape=shape, dtype=SizeType)
        if not isinstance(keys, tuple):
            keys = (keys,)
        nkeys = len(keys)
        dim = 0
        # Here is some problem when dealing with [...,...] params
        # but this is a bit weird way to pass parameters anyway
        for key in keys:
            ellipsis = 0  # Sentinel
            if isinstance(key, type(Ellipsis)):
                ellipsis = 1
                for diml in range(dim, len(self.shape) - (nkeys - dim) + 1):
                    startl[dim] = 0
                    stopl[dim] = self.shape[diml]
                    stepl[dim] = 1
                    dim += 1
            elif dim >= maxlen:
                raise IndexError(
                    "Too many indices for object '%s'" % self._v_pathname
                )
            elif is_idx(key):
                key = operator.index(key)

                # Protection for index out of range
                if key >= self.shape[dim]:
                    raise IndexError("Index out of range")
                if key < 0:
                    # To support negative values (Fixes bug #968149)
                    key += self.shape[dim]
                start, stop, step = self._process_range(
                    key, key + 1, 1, dim=dim
                )
                stop_none[dim] = 1
            elif isinstance(key, slice):
                start, stop, step = self._process_range(
                    key.start, key.stop, key.step, dim=dim
                )
            else:
                raise TypeError("Non-valid index or slice: %s" % key)
            if not ellipsis:
                startl[dim] = start
                stopl[dim] = stop
                stepl[dim] = step
                dim += 1

        # Complete the other dimensions, if needed
        if dim < len(self.shape):
            for diml in range(dim, len(self.shape)):
                startl[dim] = 0
                stopl[dim] = self.shape[diml]
                stepl[dim] = 1
                dim += 1

        # Compute the shape for the container properly. Fixes #1288792
        shape = []
        for dim in range(len(self.shape)):
            new_dim = len(range(startl[dim], stopl[dim], stepl[dim]))
            if not (new_dim == 1 and stop_none[dim]):
                shape.append(new_dim)

        return startl, stopl, stepl, shape

    def _fancy_selection(self, args: list[int | list[int]]) -> tuple[
        list[tuple[int, int, int, int, str]],
        tuple[int, np.ndarray] | None,
        tuple[int, ...],
    ]:
        """Perform a NumPy-style fancy selection in `self`.

        Implements advanced NumPy-style selection operations in
        addition to the standard slice-and-int behavior.

        Indexing arguments may be ints, slices or lists of indices.

        Note: This is a backport from the h5py project.

        """
        # Internal functions

        def validate_number(num: int, length: int) -> None:
            """Validate a list member for the given axis length."""
            try:
                num = int(num)
            except TypeError:
                raise TypeError("Illegal index: %r" % num)
            if num > length - 1:
                raise IndexError("Index out of bounds: %d" % num)

        def expand_ellipsis(
            args: tuple[int | list[int], ...], rank: int
        ) -> list:
            """Expand ellipsis objects and fill in missing axes."""
            n_el = sum(1 for arg in args if arg is Ellipsis)
            if n_el > 1:
                raise IndexError("Only one ellipsis may be used.")
            elif n_el == 0 and len(args) != rank:
                args = args + (Ellipsis,)

            final_args = []
            n_args = len(args)
            for idx, arg in enumerate(args):
                if arg is Ellipsis:
                    final_args.extend((slice(None),) * (rank - n_args + 1))
                else:
                    final_args.append(arg)

            if len(final_args) > rank:
                raise IndexError("Too many indices.")

            return final_args

        def translate_slice(exp: slice, length: int) -> tuple[int, int, int]:
            """Given a slice object, return a 3-tuple (start, count, step).

            This is for use with the hyperslab selection routines.

            """
            start, stop, step = exp.start, exp.stop, exp.step
            if start is None:
                start = 0
            else:
                start = int(start)
            if stop is None:
                stop = length
            else:
                stop = int(stop)
            if step is None:
                step = 1
            else:
                step = int(step)

            if step < 1:
                raise IndexError("Step must be >= 1 (got %d)" % step)
            if stop == start:
                raise IndexError("Zero-length selections are not allowed")
            if stop < start:
                raise IndexError("Reverse-order selections are not allowed")
            if start < 0:
                start = length + start
            if stop < 0:
                stop = length + stop

            if not 0 <= start <= (length - 1):
                raise IndexError(
                    "Start index %s out of range (0-%d)" % (start, length - 1)
                )
            if not 1 <= stop <= length:
                raise IndexError(
                    "Stop index %s out of range (1-%d)" % (stop, length)
                )

            count = (stop - start) // step
            if (stop - start) % step != 0:
                count += 1

            if start + count > length:
                raise IndexError(
                    f"Selection out of bounds ({start + count}; "
                    f"axis has {length})"
                )

            return start, count, step

        # Main code for _fancy_selection
        mshape = []
        selection = []

        if not isinstance(args, tuple):
            args = (args,)

        args = expand_ellipsis(args, len(self.shape))

        list_seen = False
        reorder = None
        for idx, (exp, length) in enumerate(zip(args, self.shape)):
            if isinstance(exp, slice):
                start, count, step = translate_slice(exp, length)
                selection.append((start, count, step, idx, "AND"))
                mshape.append(count)
            else:
                try:
                    exp = list(exp)
                except TypeError:
                    exp = [exp]  # Handle scalar index as a list of length 1
                    mshape.append(0)  # Keep track of scalar index for NumPy
                else:
                    mshape.append(len(exp))
                if len(exp) == 0:
                    raise IndexError(
                        f"Empty selections are not allowed (axis {idx})"
                    )
                elif len(exp) > 1:
                    if list_seen:
                        raise IndexError("Only one selection list is allowed")
                    else:
                        list_seen = True
                else:
                    if not isinstance(exp[0], (int, np.integer)) or (
                        isinstance(exp[0], np.ndarray)
                        and not np.issubdtype(exp[0].dtype, np.integer)
                    ):
                        raise TypeError("Only integer coordinates allowed.")

                nexp = np.asarray(exp, dtype="i8")
                # Convert negative values
                nexp = np.where(nexp < 0, length + nexp, nexp)
                # Check whether the list is ordered or not
                # (only one unordered list is allowed)
                if len(nexp) != len(np.unique(nexp)):
                    raise IndexError(
                        "Selection lists cannot have repeated values. "
                        "To see how to handle this, please see "
                        "https://github.com/PyTables/PyTables/issues/1149"
                    )
                neworder = nexp.argsort()
                if (
                    neworder.shape != (len(exp),)
                    or np.sum(np.abs(neworder - np.arange(len(exp)))) != 0
                ):
                    if reorder is not None:
                        raise IndexError(
                            "Only one selection list can be unordered"
                        )
                    corrected_idx = sum(1 for x in mshape if x != 0) - 1
                    reorder = (corrected_idx, neworder)
                    nexp = nexp[neworder]
                for select_idx in range(len(nexp) + 1):
                    # This crazy piece of code performs a list selection
                    # using HDF5 hyperslabs.
                    # For each index, perform a "NOTB" selection on every
                    # portion of *this axis* which falls *outside* the list
                    # selection.  For this to work, the input array MUST be
                    # monotonically increasing.
                    if select_idx < len(nexp):
                        validate_number(nexp[select_idx], length)
                    if select_idx == 0:
                        start = 0
                        count = nexp[0]
                    elif select_idx == len(nexp):
                        start = nexp[-1] + 1
                        count = length - start
                    else:
                        start = nexp[select_idx - 1] + 1
                        count = nexp[select_idx] - start
                    if count > 0:
                        selection.append((start, count, 1, idx, "NOTB"))

        mshape = tuple(x for x in mshape if x != 0)
        return selection, reorder, mshape

    def __getitem__(self, key: SelectionType) -> list | np.ndarray:
        """Get a row, a range of rows or a slice from the array.

        The set of tokens allowed for the key is the same as that for extended
        slicing in Python (including the Ellipsis or ... token).  The result is
        an object of the current flavor; its shape depends on the kind of slice
        used as key and the shape of the array itself.

        Furthermore, NumPy-style fancy indexing, where a list of indices in a
        certain axis is specified, is also supported.  Note that only one list
        per selection is supported right now.  Finally, NumPy-style point and
        boolean selections are supported as well.

        Examples
        --------
        ::

            array1 = array[4]                       # simple selection
            array2 = array[4:1000:2]                # slice selection
            array3 = array[1, ..., ::2, 1:4, 4:]    # general slice selection
            array4 = array[1, [1,5,10], ..., -1]    # fancy selection
            array5 = array[np.where(array[:] > 4)]  # point selection
            array6 = array[array[:] > 4]            # boolean selection

        """
        self._g_check_open()

        try:
            # First, try with a regular selection
            startl, stopl, stepl, shape = self._interpret_indexing(key)
            arr = self._read_slice(startl, stopl, stepl, shape)
        except TypeError:
            # Then, try with a point-wise selection
            try:
                coords = self._point_selection(key)
                arr = self._read_coords(coords)
            except TypeError:
                # Finally, try with a fancy selection
                selection, reorder, shape = self._fancy_selection(key)
                arr = self._read_selection(selection, reorder, shape)

        if self.flavor == "numpy" or not self._v_convert:
            return arr

        return internal_to_flavor(arr, self.flavor)

    def __setitem__(self, key: SelectionType, value: Any) -> None:
        """Set a row, a range of rows or a slice in the array.

        It takes different actions depending on the type of the key parameter:
        if it is an integer, the corresponding array row is set to value (the
        value is broadcast when needed).  If key is a slice, the row slice
        determined by it is set to value (as usual, if the slice to be updated
        exceeds the actual shape of the array, only the values in the existing
        range are updated).

        If value is a multidimensional object, then its shape must be
        compatible with the shape determined by key, otherwise, a ValueError
        will be raised.

        Furthermore, NumPy-style fancy indexing, where a list of indices in a
        certain axis is specified, is also supported.  Note that only one list
        per selection is supported right now.  Finally, NumPy-style point and
        boolean selections are supported as well.

        Examples
        --------
        ::

            a1[0] = 333        # assign an integer to an Integer Array row
            a2[0] = 'b'        # assign a string to a string Array row
            a3[1:4] = 5        # broadcast 5 to slice 1:4
            a4[1:4:2] = 'xXx'  # broadcast 'xXx' to slice 1:4:2

            # General slice update (a5.shape = (4,3,2,8,5,10).
            a5[1, ..., ::2, 1:4, 4:] = numpy.arange(1728, shape=(4,3,2,4,3,6))
            a6[1, [1,5,10], ..., -1] = arr    # fancy selection
            a7[np.where(a6[:] > 4)] = 4       # point selection + broadcast
            a8[arr > 4] = arr2                # boolean selection

        """
        self._g_check_open()

        # Create an array compliant with the specified slice
        nparr = convert_to_np_atom2(value, self.atom)
        if nparr.size == 0:
            return

        # truncate data if least_significant_digit filter is set
        # TODO: add the least_significant_digit attribute to the array on disk
        if (
            self.filters.least_significant_digit is not None
            and not np.issubdtype(nparr.dtype, np.signedinteger)
        ):
            nparr = quantize(nparr, self.filters.least_significant_digit)

        try:
            startl, stopl, stepl, shape = self._interpret_indexing(key)
            self._write_slice(startl, stopl, stepl, shape, nparr)
        except TypeError:
            # Then, try with a point-wise selection
            try:
                coords = self._point_selection(key)
                self._write_coords(coords, nparr)
            except TypeError:
                selection, reorder, shape = self._fancy_selection(key)
                self._write_selection(selection, reorder, shape, nparr)

    def _check_shape(
        self, nparr: np.ndarray, slice_shape: tuple[int, ...]
    ) -> np.ndarray:
        """Test that nparr shape is consistent with underlying object.

        If not, try creating a new nparr object, using broadcasting if
        necessary.

        """
        if nparr.shape != (slice_shape + self.atom.dtype.shape):
            # Create an array compliant with the specified shape
            narr = np.empty(shape=slice_shape, dtype=self.atom.dtype)

            # Assign the value to it. It will raise a ValueError exception
            # if the objects cannot be broadcast to a single shape.
            narr[...] = nparr
            return narr
        else:
            return nparr

    def _read_slice(
        self,
        startl: np.ndarray,
        stopl: np.ndarray,
        stepl: np.ndarray,
        shape: list[int],
    ) -> np.ndarray:
        """Read a slice based on `startl`, `stopl` and `stepl`."""
        nparr = np.empty(dtype=self.atom.dtype, shape=shape)
        # Protection against reading empty arrays
        if 0 not in shape:
            # Arrays that have non-zero dimensionality
            self._g_read_slice(startl, stopl, stepl, nparr)
        # For zero-shaped arrays, return the scalar
        if nparr.shape == ():
            nparr = nparr[()]
        return nparr

    def _read_coords(self, coords: np.ndarray) -> np.ndarray:
        """Read a set of points defined by `coords`."""
        nparr = np.empty(dtype=self.atom.dtype, shape=len(coords))
        if len(coords) > 0:
            self._g_read_coords(coords, nparr)
        # For zero-shaped arrays, return the scalar
        if nparr.shape == ():
            nparr = nparr[()]
        return nparr

    def _read_selection(
        self,
        selection: list[tuple[int, int, int, int, str]],
        reorder: tuple[int, npt.ArrayLike] | None,
        shape: tuple[int, ...],
    ) -> np.ndarray:
        """Read a `selection`.

        Reorder if necessary.

        """
        # Create the container for the slice
        nparr = np.empty(dtype=self.atom.dtype, shape=shape)
        # Arrays that have non-zero dimensionality
        self._g_read_selection(selection, nparr)
        # For zero-shaped arrays, return the scalar
        if nparr.shape == ():
            nparr = nparr[()]
        elif reorder is not None:
            # We need to reorder the array
            idx, neworder = reorder
            k = [slice(None)] * len(shape)
            k[idx] = neworder.argsort()
            # Apparently, a copy is not needed here, but doing it
            # for symmetry with the `_write_selection()` method.
            nparr = nparr[tuple(k)].copy()
        return nparr

    def _write_slice(
        self,
        startl: np.ndarray,
        stopl: np.ndarray,
        stepl: np.ndarray,
        shape: list[int],
        nparr: np.ndarray,
    ) -> None:
        """Write `nparr` in a slice based on `startl`, `stopl` and `stepl`."""
        nparr = self._check_shape(nparr, tuple(shape))
        countl = ((stopl - startl - 1) // stepl) + 1
        self._g_write_slice(startl, stepl, countl, nparr)

    def _write_coords(self, coords: np.ndarray, nparr: np.ndarray) -> None:
        """Write `nparr` values in points defined by `coords` coordinates."""
        if len(coords) > 0:
            nparr = self._check_shape(nparr, (len(coords),))
            self._g_write_coords(coords, nparr)

    def _write_selection(
        self,
        selection: list[tuple[int, int, int, int, str]],
        reorder: tuple[int, npt.ArrayLike] | None,
        shape: tuple[int, ...],
        nparr: np.ndarray,
    ) -> None:
        """Write `nparr` in `selection`.

        Reorder if necessary.

        """
        nparr = self._check_shape(nparr, tuple(shape))
        # Check whether we should reorder the array
        if reorder is not None:
            idx, neworder = reorder
            k = [slice(None)] * len(shape)
            k[idx] = neworder
            # For a reason a don't understand well, we need a copy of
            # the reordered array
            nparr = nparr[tuple(k)].copy()
        self._g_write_selection(selection, nparr)

    def _read(
        self, start: int, stop: int, step: int, out: np.ndarray | None = None
    ) -> np.ndarray:
        """Read the array from disk without slice or flavor processing."""
        nrowstoread = len(range(start, stop, step))
        shape = list(self.shape)
        if shape:
            shape[self.maindim] = nrowstoread
        if out is None:
            arr = np.empty(dtype=self.atom.dtype, shape=shape)
        else:
            bytes_required = self.rowsize * nrowstoread
            # if buffer is too small, it will segfault
            if bytes_required != out.nbytes:
                raise ValueError(
                    f"output array size invalid, got {out.nbytes}"
                    f" bytes, need {bytes_required} bytes"
                )
            if not out.flags["C_CONTIGUOUS"]:
                raise ValueError("output array not C contiguous")
            arr = out
        # Protection against reading empty arrays
        if 0 not in shape:
            # Arrays that have non-zero dimensionality
            self._read_array(start, stop, step, arr)
        # data is always read in the system byteorder
        # if the out array's byteorder is different, do a byteswap
        if (
            out is not None
            and byteorders[arr.dtype.byteorder] != sys.byteorder
        ):
            arr.byteswap(True)
        return arr

    def read(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        out: np.ndarray | None = None,
    ) -> np.ndarray:
        """Get data in the array as an object of the current flavor.

        The start, stop and step parameters can be used to select only a
        *range of rows* in the array.  Their meanings are the same as in
        the built-in range() Python function, except that negative values
        of step are not allowed yet. Moreover, if only start is specified,
        then stop will be set to start + 1. If you do not specify neither
        start nor stop, then *all the rows* in the array are selected.

        The out parameter may be used to specify a NumPy array to receive
        the output data.  Note that the array must have the same size as
        the data selected with the other parameters.  Note that the array's
        datatype is not checked and no type casting is performed, so if it
        does not match the datatype on disk, the output will not be correct.
        Also, this parameter is only valid when the array's flavor is set
        to 'numpy'.  Otherwise, a TypeError will be raised.

        When data is read from disk in NumPy format, the output will be
        in the current system's byteorder, regardless of how it is stored
        on disk.
        The exception is when an output buffer is supplied, in which case
        the output will be in the byteorder of that output buffer.

        .. versionchanged:: 3.0
           Added the *out* parameter.

        """
        self._g_check_open()
        if out is not None and self.flavor != "numpy":
            msg = (
                f"Optional 'out' argument may only be supplied if array "
                f"flavor is 'numpy', currently is {self.flavor}"
            )
            raise TypeError(msg)
        (start, stop, step) = self._process_range_read(start, stop, step)
        arr = self._read(start, stop, step, out)
        return internal_to_flavor(arr, self.flavor)

    def _g_copy_with_stats(
        self,
        group: Group,
        name: str,
        start: int,
        stop: int,
        step: int,
        title: str,
        filters: Filters,
        chunkshape: tuple[int, ...],
        _log: bool,
        **kwargs,
    ) -> tuple[Array, int]:
        """Private part of Leaf.copy() for each kind of leaf."""
        # Compute the correct indices.
        (start, stop, step) = self._process_range_read(start, stop, step)
        # Get the slice of the array
        # (non-buffered version)
        if self.shape:
            arr = self[start:stop:step]
        else:
            arr = self[()]
        # Build the new Array object.  Use the _atom reserved keyword
        # just in case the array is being copied from a native HDF5
        # with atomic types different from scalars.
        # For details, see #275 of trac.
        object_ = Array(
            group, name, arr, title=title, _log=_log, _atom=self.atom
        )
        nbytes = np.prod(self.shape, dtype=SizeType) * self.atom.size

        return (object_, nbytes)

    def __repr__(self) -> str:
        """Provide more metainfo in addition to standard __str__."""
        return f"""{self}
  atom := {self.atom!r}
  maindim := {self.maindim!r}
  flavor := {self.flavor!r}
  byteorder := {self.byteorder!r}
  chunkshape := {self.chunkshape!r}"""


class ImageArray(Array):
    """Array containing an image.

    This class has no additional behaviour or functionality compared to
    that of an ordinary array.  It simply enables the user to open an
    ``IMAGE`` HDF5 node as a normal `Array` node in PyTables.

    """

    # Class identifier.
    _c_classid = "IMAGE"
