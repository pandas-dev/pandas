"""Here is defined the VLArray class."""

from __future__ import annotations

import sys
import operator
from typing import Any, NoReturn, TYPE_CHECKING
from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

from . import hdf5extension
from .atom import ObjectAtom, VLStringAtom, VLUnicodeAtom
from .leaf import Leaf, calc_chunksize
from .utils import (
    convert_to_np_atom,
    convert_to_np_atom2,
    idx2long,
    correct_byteorder,
    SizeType,
    is_idx,
    lazyattr,
)
from .flavor import internal_to_flavor

if TYPE_CHECKING:
    from .atom import Atom, Enum
    from .group import Group
    from .filters import Filters

# default version for VLARRAY objects
# obversion = "1.0"    # initial version
# obversion = "1.0"    # add support for complex datatypes
# obversion = "1.1"    # This adds support for time datatypes.
# obversion = "1.2"    # This adds support for enumerated datatypes.
# obversion = "1.3"     # Introduced 'PSEUDOATOM' attribute.
obversion = "1.4"  # Numeric and numarray flavors are gone.


class VLArray(hdf5extension.VLArray, Leaf):
    """This class represents variable length (ragged) arrays in an HDF5 file.

    Instances of this class represent array objects in the object tree
    with the property that their rows can have a *variable* number of
    homogeneous elements, called *atoms*. Like Table datasets (see
    :ref:`TableClassDescr`), variable length arrays can have only one
    dimension, and the elements (atoms) of their rows can be fully
    multidimensional.

    When reading a range of rows from a VLArray, you will *always* get
    a Python list of objects of the current flavor (each of them for a
    row), which may have different lengths.

    This class provides methods to write or read data to or from
    variable length array objects in the file. Note that it also
    inherits all the public attributes and methods that Leaf (see
    :ref:`LeafClassDescr`) already provides.

    .. note::

          VLArray objects also support compression although compression
          is only performed on the data structures used internally by
          the HDF5 to take references of the location of the variable
          length data. Data itself (the raw data) are not compressed
          or filtered.

          Please refer to the `VLTypes Technical Note
          <https://support.hdfgroup.org/HDF5/doc/TechNotes/VLTypes.html>`_
          for more details on the topic.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.
    name : str
        The name of this node in its parent group.
    atom
        An `Atom` instance representing the *type* and *shape* of the atomic
        objects to be saved.
    title
        A description for this node (it sets the ``TITLE`` HDF5 attribute on
        disk).
    filters
        An instance of the `Filters` class that provides information about the
        desired I/O filters to be applied during the life of this object.
    expectedrows
        A user estimate about the number of row elements that will
        be added to the growable dimension in the `VLArray` node.
        If not provided, the default value is ``EXPECTED_ROWS_VLARRAY``
        (see ``tables/parameters.py``).  If you plan to create either
        a much smaller or a much bigger `VLArray` try providing a guess;
        this will optimize the HDF5 B-Tree creation and management
        process time and the amount of memory used.

        .. versionadded:: 3.0

    chunkshape
        The shape of the data chunk to be read or written in a single HDF5 I/O
        operation.  Filters are applied to those chunks of data.  The
        dimensionality of `chunkshape` must be 1.  If ``None``, a sensible
        value is calculated (which is recommended).
    byteorder
        The byteorder of the data *on disk*, specified as 'little' or 'big'.
        If this is not specified, the byteorder is that of the platform.

    track_times
        Whether time data associated with the leaf are recorded (object
        access time, raw data modification time, metadata change time, object
        birth time); default True.  Semantics of these times depend on their
        implementation in the HDF5 library: refer to documentation of the
        H5O_info_t data structure.  As of HDF5 1.8.15, only ctime (metadata
        change time) is implemented.

        .. versionadded:: 3.4.3


    .. versionchanged:: 3.0
       *parentNode* renamed into *parentnode*.

    .. versionchanged:: 3.0
       The *expectedsizeinMB* parameter has been replaced by *expectedrows*.

    Examples
    --------
    See below a small example of the use of the VLArray class.  The code is
    available in :file:`examples/vlarray1.py`::

        import numpy as np
        import tables as tb

        # Create a VLArray:
        fileh = tb.open_file('vlarray1.h5', mode='w')
        vlarray = fileh.create_vlarray(
            fileh.root,
            'vlarray1',
            tb.Int32Atom(shape=()),
            "ragged array of ints",
            filters=tb.Filters(1))

        # Append some (variable length) rows:
        vlarray.append(np.array([5, 6]))
        vlarray.append(np.array([5, 6, 7]))
        vlarray.append([5, 6, 9, 8])

        # Now, read it through an iterator:
        print('-->', vlarray.title)
        for x in vlarray:
            print('%s[%d]--> %s' % (vlarray.name, vlarray.nrow, x))

        # Now, do the same with native Python strings.
        vlarray2 = fileh.create_vlarray(
            fileh.root,
            'vlarray2',
            tb.StringAtom(itemsize=2),
            "ragged array of strings",
            filters=tb.Filters(1))
        vlarray2.flavor = 'python'

        # Append some (variable length) rows:
        print('-->', vlarray2.title)
        vlarray2.append(['5', '66'])
        vlarray2.append(['5', '6', '77'])
        vlarray2.append(['5', '6', '9', '88'])

        # Now, read it through an iterator:
        for x in vlarray2:
            print('%s[%d]--> %s' % (vlarray2.name, vlarray2.nrow, x))

        # Close the file.
        fileh.close()

    The output for the previous script is something like::

        --> ragged array of ints
        vlarray1[0]--> [5 6]
        vlarray1[1]--> [5 6 7]
        vlarray1[2]--> [5 6 9 8]
        --> ragged array of strings
        vlarray2[0]--> ['5', '66']
        vlarray2[1]--> ['5', '6', '77']
        vlarray2[2]--> ['5', '6', '9', '88']


    .. rubric:: VLArray attributes

    The instance variables below are provided in addition to those in
    Leaf (see :ref:`LeafClassDescr`).

    .. attribute:: atom

        An Atom (see :ref:`AtomClassDescr`)
        instance representing the *type* and
        *shape* of the atomic objects to be
        saved. You may use a *pseudo-atom* for
        storing a serialized object or variable length string per row.

    .. attribute:: flavor

        The type of data object read from this leaf.

        Please note that when reading several rows of VLArray data,
        the flavor only applies to the *components* of the returned
        Python list, not to the list itself.

    .. attribute:: nrow

        On iterators, this is the index of the current row.

    .. attribute:: nrows

        The current number of rows in the array.

    .. attribute:: extdim

       The index of the enlargeable dimension (always 0 for vlarrays).

    """

    # Class identifier.
    _c_classid = "VLARRAY"

    @lazyattr
    def dtype(self) -> np.dtype:
        """Return the NumPy ``dtype`` that most closely matches this array."""
        return self.atom.dtype

    @property
    def shape(self) -> tuple[int]:
        """Return the shape of the stored array."""
        return (self.nrows,)

    @property
    def size_on_disk(self) -> NoReturn:
        """Return the size on disk of the `VLArray` object.

        The HDF5 library does not include a function to determine size_on_disk
        for variable-length arrays.  Accessing this attribute will raise a
        NotImplementedError.
        """
        raise NotImplementedError("size_on_disk not implemented for VLArrays")

    @property
    def size_in_memory(self) -> int:
        """Size of the array's data in bytes.

        .. note::

            When data is stored in a VLArray using the ObjectAtom type,
            it is first serialized using pickle, and then converted to
            a NumPy array suitable for storage in an HDF5 file.
            This attribute will return the size of that NumPy
            representation.  If you wish to know the size of the Python
            objects after they are loaded from disk, you can use this
            `ActiveState recipe
            <http://code.activestate.com/recipes/577504/>`_.
        """
        return self._get_memory_size()

    def __init__(
        self,
        parentnode: Group,
        name: str,
        atom: Atom | None = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int | None = None,
        chunkshape: tuple[int, ...] | None = None,
        byteorder: str | None = None,
        _log: bool = True,
        track_times: bool = True,
    ) -> None:

        self._v_version: str | None = None
        """The object version of this array."""

        self._v_new = new = atom is not None
        """Is this the first time the node has been created?"""

        self._v_new_title = title
        """New title for this node."""

        self._v_new_filters = filters
        """New filter properties for this array."""

        if expectedrows is None:
            expectedrows = parentnode._v_file.params["EXPECTED_ROWS_VLARRAY"]
        self._v_expectedrows = int(expectedrows)
        """The expected number of rows to be stored in the array.

        .. versionadded:: 3.0

        """

        self._v_chunkshape: tuple[int, ...] | None = None
        """Private storage for the `chunkshape` property of Leaf."""

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
        self.atom = atom
        """
        An Atom (see :ref:`AtomClassDescr`) instance representing the
        *type* and *shape* of the atomic objects to be saved. You may
        use a *pseudo-atom* for storing a serialized object or
        variable length string per row.
        """
        self.nrow: int | None = None
        """On iterators, this is the index of the current row."""

        self.nrows: int | None = None
        """The current number of rows in the array."""

        self.extdim = 0  # VLArray only have one dimension currently
        """The index of the enlargeable dimension (always 0 for vlarrays)."""

        # Check the chunkshape parameter
        if new and chunkshape is not None:
            if isinstance(chunkshape, (int, np.integer)):
                chunkshape = (chunkshape,)
            try:
                chunkshape = tuple(chunkshape)
            except TypeError:
                raise TypeError(
                    "`chunkshape` parameter must be an integer or sequence "
                    "and you passed a %s" % type(chunkshape)
                )
            if len(chunkshape) != 1:
                raise ValueError(
                    f"`chunkshape` rank (length) must be 1: {chunkshape!r}"
                )
            self._v_chunkshape = tuple(SizeType(s) for s in chunkshape)

        super().__init__(
            parentnode, name, new, filters, byteorder, _log, track_times
        )

    def _g_post_init_hook(self) -> None:
        super()._g_post_init_hook()
        self.nrowsinbuf = 100  # maybe enough for most applications

    # This is too specific for moving it into Leaf
    def _calc_chunkshape(self, expectedrows: int) -> tuple[int]:
        """Calculate the size for the HDF5 chunk."""
        # For computing the chunkshape for HDF5 VL types, we have to
        # choose the itemsize of the *each* element of the atom and
        # not the size of the entire atom.  I don't know why this
        # should be like this, perhaps I should report this to the
        # HDF5 list.
        # F. Alted 2006-11-23
        # elemsize = self.atom.atomsize()
        elemsize = self._basesize

        # AV 2013-05-03
        # This is just a quick workaround tha allows to change the API for
        # PyTables 3.0 release and remove the expected_mb parameter.
        # The algorithm for computing the chunkshape should be rewritten as
        # requested by gh-35.
        expected_mb = expectedrows * elemsize / 1024**2

        chunksize = calc_chunksize(expected_mb)

        # Set the chunkshape
        chunkshape = chunksize // elemsize
        # Safeguard against itemsizes being extremely large
        if chunkshape == 0:
            chunkshape = 1
        return (SizeType(chunkshape),)

    def _g_create(self) -> int:
        """Create a variable length array (ragged array)."""
        atom = self.atom
        self._v_version = obversion
        # Check for zero dims in atom shape (not allowed in VLArrays)
        zerodims = np.sum(np.array(atom.shape) == 0)
        if zerodims > 0:
            raise ValueError(
                "When creating VLArrays, none of the dimensions "
                "of the Atom instance can be zero."
            )

        if not hasattr(atom, "size"):  # it is a pseudo-atom
            self._atomicdtype = atom.base.dtype
            self._atomicsize = atom.base.size
            self._basesize = atom.base.itemsize
        else:
            self._atomicdtype = atom.dtype
            self._atomicsize = atom.size
            self._basesize = atom.itemsize
        self._atomictype = atom.type
        self._atomicshape = atom.shape

        # Compute the optimal chunkshape, if needed
        if self._v_chunkshape is None:
            self._v_chunkshape = self._calc_chunkshape(self._v_expectedrows)

        self.nrows = SizeType(0)  # No rows at creation time

        # Correct the byteorder if needed
        if self.byteorder is None:
            self.byteorder = correct_byteorder(atom.type, sys.byteorder)

        # After creating the vlarray, ``self._v_objectid`` needs to be
        # set because it is needed for setting attributes afterwards.
        self._v_objectid = self._create_array(self._v_new_title)

        # Add an attribute in case we have a pseudo-atom so that we
        # can retrieve the proper class after a re-opening operation.
        if not hasattr(atom, "size"):  # it is a pseudo-atom
            self.attrs.PSEUDOATOM = atom.kind

        return self._v_objectid

    def _g_open(self) -> int:
        """Get the metadata info for an array in file."""
        self._v_objectid, self.nrows, self._v_chunkshape, atom = (
            self._open_array()
        )

        # Check if the atom can be a PseudoAtom
        if "PSEUDOATOM" in self.attrs:
            kind = self.attrs.PSEUDOATOM
            if kind == "vlstring":
                atom = VLStringAtom()
            elif kind == "vlunicode":
                atom = VLUnicodeAtom()
            elif kind == "object":
                atom = ObjectAtom()
            else:
                raise ValueError("pseudo-atom name ``%s`` not known." % kind)
        elif self._v_file.format_version[:1] == "1":
            flavor1x = self.attrs.FLAVOR
            if flavor1x == "VLString":
                atom = VLStringAtom()
            elif flavor1x == "Object":
                atom = ObjectAtom()

        self.atom = atom
        return self._v_objectid

    def _getnobjects(self, nparr: np.ndarray) -> int:
        """Return the number of objects in a NumPy array."""
        # Check for zero dimensionality array
        zerodims = np.sum(np.array(nparr.shape) == 0)
        if zerodims > 0:
            # No objects to be added
            return 0
        shape = nparr.shape
        atom_shape = self.atom.shape
        shapelen = len(nparr.shape)
        if isinstance(atom_shape, tuple):
            atomshapelen = len(self.atom.shape)
        else:
            atom_shape = (self.atom.shape,)
            atomshapelen = 1
        diflen = shapelen - atomshapelen
        if shape == atom_shape:
            nobjects = 1
        elif diflen == 1 and shape[diflen:] == atom_shape:
            # Check if the leading dimensions are all ones
            # if shape[:diflen-1] == (1,)*(diflen-1):
            #    nobjects = shape[diflen-1]
            #    shape = shape[diflen:]
            # It's better to accept only inputs with the exact dimensionality
            # i.e. a dimensionality only 1 element larger than atom
            nobjects = shape[0]
            shape = shape[1:]
        elif atom_shape == (1,) and shapelen == 1:
            # Case where shape = (N,) and shape_atom = 1 or (1,)
            nobjects = shape[0]
        else:
            raise ValueError(
                "The object '%s' is composed of elements with "
                "shape '%s', which is not compatible with the "
                "atom shape ('%s')." % (nparr, shape, atom_shape)
            )
        return nobjects

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

    def append(self, sequence: npt.ArrayLike) -> None:
        """Add a sequence of data to the end of the dataset.

        This method appends the objects in the sequence to a *single row* in
        this array. The type and shape of individual objects must be compliant
        with the atoms in the array. In the case of serialized objects and
        variable length strings, the object or string to append is itself the
        sequence.

        """
        self._g_check_open()
        self._v_file._check_writable()

        # Prepare the sequence to convert it into a NumPy object
        atom = self.atom
        if not hasattr(atom, "size"):  # it is a pseudo-atom
            sequence = atom.toarray(sequence)
            statom = atom.base
        else:
            try:  # fastest check in most cases
                len(sequence)
            except TypeError:
                raise TypeError("argument is not a sequence")
            statom = atom

        if len(sequence) > 0:
            # The sequence needs to be copied to make the operation safe
            # to in-place conversion.
            nparr = convert_to_np_atom2(sequence, statom)
            nobjects = self._getnobjects(nparr)
        else:
            nobjects = 0
            nparr = None

        self._append(nparr, nobjects)
        self.nrows += 1

    def iterrows(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> VLArray:
        """Iterate over the rows of the array.

        This method returns an iterator yielding an object of the current
        flavor for each selected row in the array.

        If a range is not supplied, *all the rows* in the array are iterated
        upon. You can also use the :meth:`VLArray.__iter__` special method for
        that purpose.  If you only want to iterate over a given *range of rows*
        in the array, you may use the start, stop and step parameters.

        Examples
        --------
        ::

            for row in vlarray.iterrows(step=4):
                print('%s[%d]--> %s' % (vlarray.name, vlarray.nrow, row))

        .. versionchanged:: 3.0
           If the *start* parameter is provided and *stop* is None then the
           array is iterated from *start* to the last line.
           In PyTables < 3.0 only one element was returned.

        """
        (self._start, self._stop, self._step) = self._process_range(
            start, stop, step
        )
        self._init_loop()
        return self

    def __iter__(self) -> VLArray:
        """Iterate over the rows of the array.

        This is equivalent to calling :meth:`VLArray.iterrows` with default
        arguments, i.e. it iterates over *all the rows* in the array.

        Examples
        --------
        ::

            result = [row for row in vlarray]

        Which is equivalent to::

            result = [row for row in vlarray.iterrows()]

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

    def __next__(self) -> list | np.ndarray:
        """Get the next element of the array during an iteration.

        The element is returned as a list of objects of the current
        flavor.

        """
        if self._nrowsread >= self._stop:
            self._init = False
            raise StopIteration  # end of iteration
        else:
            # Read a chunk of rows
            if self._row + 1 >= self.nrowsinbuf or self._row < 0:
                self._stopb = self._startb + self._step * self.nrowsinbuf
                self.listarr = self.read(self._startb, self._stopb, self._step)
                self._row = -1
                self._startb = self._stopb
            self._row += 1
            self.nrow += self._step
            self._nrowsread += self._step
            return self.listarr[self._row]

    def __getitem__(
        self, key: int | slice | Sequence[int] | np.ndarray
    ) -> list:
        """Get a row or a range of rows from the array.

        If key argument is an integer, the corresponding array row is returned
        as an object of the current flavor.  If key is a slice, the range of
        rows determined by it is returned as a list of objects of the current
        flavor.

        In addition, NumPy-style point selections are supported.  In
        particular, if key is a list of row coordinates, the set of rows
        determined by it is returned.  Furthermore, if key is an array of
        boolean values, only the coordinates where key is True are returned.
        Note that for the latter to work it is necessary that key list would
        contain exactly as many rows as the array has.

        Examples
        --------
        ::

            a_row = vlarray[4]
            a_list = vlarray[4:1000:2]
            a_list2 = vlarray[[0,2]]   # get list of coords
            a_list3 = vlarray[[0,-2]]  # negative values accepted
            a_list4 = vlarray[np.array([True,...,False])]  # array of bools

        """
        self._g_check_open()
        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= self.nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += self.nrows
            (start, stop, step) = self._process_range(key, key + 1, 1)
            return self.read(start, stop, step)[0]
        elif isinstance(key, slice):
            start, stop, step = self._process_range(
                key.start, key.stop, key.step
            )
            return self.read(start, stop, step)
        # Try with a boolean or point selection
        elif type(key) in (list, tuple) or isinstance(key, np.ndarray):
            coords = self._point_selection(key)
            return self._read_coordinates(coords)
        else:
            raise IndexError(f"Invalid index or slice: {key!r}")

    def _assign_values(self, coords: Sequence[int], values: Sequence) -> None:
        """Assign the `values` to the positions stated in `coords`."""
        for nrow, value in zip(coords, values):
            if nrow >= self.nrows:
                raise IndexError("First index out of range")
            if nrow < 0:
                # To support negative values
                nrow += self.nrows
            object_ = value
            # Prepare the object to convert it into a NumPy object
            atom = self.atom
            if not hasattr(atom, "size"):  # it is a pseudo-atom
                object_ = atom.toarray(object_)
                statom = atom.base
            else:
                statom = atom
            value = convert_to_np_atom(object_, statom)
            nobjects = self._getnobjects(value)

            # Get the previous value
            nrow = idx2long(nrow)  # To convert any possible numpy scalar value
            nparr = self._read_array(nrow, nrow + 1, 1)[0]
            nobjects = len(nparr)
            if len(value) > nobjects:
                raise ValueError(
                    "Length of value (%s) is larger than number "
                    "of elements in row (%s)" % (len(value), nobjects)
                )
            try:
                nparr[:] = value
            except Exception as exc:  # XXX
                raise ValueError(
                    "Value parameter:\n'%r'\n"
                    "cannot be converted into an array object "
                    "compliant vlarray[%s] row: \n'%r'\n"
                    "The error was: <%s>" % (value, nrow, nparr[:], exc)
                )

            if nparr.size > 0:
                self._modify(nrow, nparr, nobjects)

    def __setitem__(
        self,
        key: int | slice | Sequence[int] | np.ndarray,
        value: Any,
    ) -> None:
        """Set a row, or set of rows, in the array.

        It takes different actions depending on the type of the *key*
        parameter: if it is an integer, the corresponding table row is
        set to *value* (a record or sequence capable of being converted
        to the table structure).  If *key* is a slice, the row slice
        determined by it is set to *value* (a record array or sequence
        of rows capable of being converted to the table structure).

        In addition, NumPy-style point selections are supported.  In
        particular, if key is a list of row coordinates, the set of rows
        determined by it is set to value.  Furthermore, if key is an array of
        boolean values, only the coordinates where key is True are set to
        values from value.  Note that for the latter to work it is necessary
        that key list would contain exactly as many rows as the table has.

        .. note::

            When updating the rows of a VLArray object which uses a
            pseudo-atom, there is a problem: you can only update values
            with *exactly* the same size in bytes than the original row.
            This is very difficult to meet with object pseudo-atoms,
            because :mod:`pickle` applied on a Python object does not
            guarantee to return the same number of bytes than over another
            object, even if they are of the same class.
            This effectively limits the kinds of objects than can be
            updated in variable-length arrays.

        Examples
        --------
        ::

            vlarray[0] = vlarray[0] * 2 + 3
            vlarray[99] = arange(96) * 2 + 3

            # Negative values for the index are supported.
            vlarray[-99] = vlarray[5] * 2 + 3
            vlarray[1:30:2] = list_of_rows
            vlarray[[1,3]] = new_1_and_3_rows

        """
        self._g_check_open()
        self._v_file._check_writable()

        if is_idx(key):
            # If key is not a sequence, convert to it
            coords = [key]
            value = [value]
        elif isinstance(key, slice):
            start, stop, step = self._process_range(
                key.start, key.stop, key.step
            )
            coords = range(start, stop, step)
        # Try with a boolean or point selection
        elif type(key) in (list, tuple) or isinstance(key, np.ndarray):
            coords = self._point_selection(key)
        else:
            raise IndexError(f"Invalid index or slice: {key!r}")

        # Do the assignment row by row
        self._assign_values(coords, value)

    # Accessor for the _read_array method in superclass
    def read(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int = 1,
    ) -> list:
        """Get data in the array as a list of objects of the current flavor.

        Please note that, as the lengths of the different rows are variable,
        the returned value is a *Python list* (not an array of the current
        flavor), with as many entries as specified rows in the range
        parameters.

        The start, stop and step parameters can be used to select only a
        *range of rows* in the array.  Their meanings are the same as in
        the built-in range() Python function, except that negative values
        of step are not allowed yet. Moreover, if only start is specified,
        then stop will be set to start + 1. If you do not specify neither
        start nor stop, then *all the rows* in the array are selected.

        """
        self._g_check_open()
        start, stop, step = self._process_range_read(start, stop, step)
        if start == stop:
            listarr = []
        else:
            listarr = self._read_array(start, stop, step)

        atom = self.atom
        if not hasattr(atom, "size"):  # it is a pseudo-atom
            outlistarr = [atom.fromarray(arr) for arr in listarr]
        else:
            # Convert the list to the right flavor
            flavor = self.flavor
            outlistarr = [internal_to_flavor(arr, flavor) for arr in listarr]
        return outlistarr

    def _read_coordinates(self, coords: Sequence[int]) -> list[list]:
        """Read rows specified in `coords`."""
        rows = []
        for coord in coords:
            rows.append(self.read(idx2long(coord), idx2long(coord) + 1, 1)[0])
        return rows

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
    ) -> tuple[VLArray, int]:
        """Private part of Leaf.copy() for each kind of leaf."""
        # Build the new VLArray object
        obj = VLArray(
            group,
            name,
            self.atom,
            title=title,
            filters=filters,
            expectedrows=self._v_expectedrows,
            chunkshape=chunkshape,
            _log=_log,
        )

        # Now, fill the new vlarray with values from the old one
        # This is not buffered because we cannot forsee the length
        # of each record. So, the safest would be a copy row by row.
        # In the future, some analysis can be done in order to buffer
        # the copy process.
        nrowsinbuf = 1
        (start, stop, step) = self._process_range_read(start, stop, step)
        # Optimized version (no conversions, no type and shape checks, etc...)
        nrowscopied = SizeType(0)
        nbytes = 0
        if not hasattr(self.atom, "size"):  # it is a pseudo-atom
            atomsize = self.atom.base.size
        else:
            atomsize = self.atom.size
        for start2 in range(start, stop, step * nrowsinbuf):
            # Save the records on disk
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            nparr = self._read_array(start=start2, stop=stop2, step=step)[0]
            nobjects = nparr.shape[0]
            obj._append(nparr, nobjects)
            nbytes += nobjects * atomsize
            nrowscopied += 1
        obj.nrows = nrowscopied
        return (obj, nbytes)

    def __repr__(self) -> str:
        """`VLArray` string representation.

        Provides more metainfo w.r.t standard __str__.
        """
        return f"""{self}
  atom = {self.atom!r}
  byteorder = {self.byteorder!r}
  nrows = {self.nrows}
  flavor = {self.flavor!r}"""
