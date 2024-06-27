"""Here is defined the Leaf class."""

import warnings
import math

import numpy as np
try:
    import cpuinfo
    cpu_info = cpuinfo.get_cpu_info()
except ImportError:
    cpu_info = {}

from .flavor import (check_flavor, internal_flavor, toarray,
                     alias_map as flavor_alias_map)
from .node import Node
from .filters import Filters
from .utils import byteorders, lazyattr, SizeType
from .exceptions import PerformanceWarning


def csformula(expected_mb):
    """Return the fitted chunksize for expected_mb."""

    # For a basesize of 8 KB, this will return:
    # 8 KB for datasets <= 1 MB
    # 1 MB for datasets >= 10 TB
    basesize = 8 * 1024  # 8 KB is a good minimum
    return basesize * int(2**math.log10(expected_mb))


def limit_es(expected_mb):
    """Protection against creating too small or too large chunks."""

    if expected_mb < 1:        # < 1 MB
        expected_mb = 1
    elif expected_mb > 10**7:  # > 10 TB
        expected_mb = 10**7
    return expected_mb


def calc_chunksize(expected_mb):
    """Compute the optimum HDF5 chunksize for I/O purposes.

    Rational: HDF5 takes the data in bunches of chunksize length to
    write the on disk. A BTree in memory is used to map structures on
    disk. The more chunks that are allocated for a dataset the larger
    the B-tree. Large B-trees take memory and causes file storage
    overhead as well as more disk I/O and higher contention for the meta
    data cache.  You have to balance between memory and I/O overhead
    (small B-trees) and time to access to data (big B-trees).

    The tuning of the chunksize parameter affects the performance and
    the memory consumed. This is based on my own experiments and, as
    always, your mileage may vary.

    """

    expected_mb = limit_es(expected_mb)
    zone = int(math.log10(expected_mb))
    expected_mb = 10**zone
    chunksize = csformula(expected_mb)
    # XXX: Multiply by 8 seems optimal for sequential access
    return chunksize * 8


class Leaf(Node):
    """Abstract base class for all PyTables leaves.

    A leaf is a node (see the Node class in :class:`Node`) which hangs from a
    group (see the Group class in :class:`Group`) but, unlike a group, it can
    not have any further children below it (i.e. it is an end node).

    This definition includes all nodes which contain actual data (datasets
    handled by the Table - see :ref:`TableClassDescr`, Array -
    see :ref:`ArrayClassDescr`, CArray - see :ref:`CArrayClassDescr`, EArray -
    see :ref:`EArrayClassDescr`, and VLArray - see :ref:`VLArrayClassDescr`
    classes) and unsupported nodes (the UnImplemented
    class - :ref:`UnImplementedClassDescr`) these classes do in fact inherit
    from Leaf.


    .. rubric:: Leaf attributes

    These instance variables are provided in addition to those in Node
    (see :ref:`NodeClassDescr`):

    .. attribute:: byteorder

        The byte ordering of the leaf data *on disk*.  It will be either
        ``little`` or ``big``.

    .. attribute:: dtype

        The NumPy dtype that most closely matches this leaf type.

    .. attribute:: extdim

        The index of the enlargeable dimension (-1 if none).

    .. attribute:: nrows

        The length of the main dimension of the leaf data.

    .. attribute:: nrowsinbuf

        The number of rows that fit in internal input buffers.

        You can change this to fine-tune the speed or memory
        requirements of your application.

    .. attribute:: shape

        The shape of data in the leaf.

    """

    # These are a little hard to override, but so are properties.
    attrs = Node._v_attrs
    """The associated AttributeSet instance - see :ref:`AttributeSetClassDescr`
    (This is an easier-to-write alias of :attr:`Node._v_attrs`."""
    title = Node._v_title
    """A description for this node
    (This is an easier-to-write alias of :attr:`Node._v_title`)."""

    @property
    def name(self):
        """The name of this node in its parent group (This is an
        easier-to-write alias of :attr:`Node._v_name`)."""
        return self._v_name

    @property
    def chunkshape(self):
        """The HDF5 chunk size for chunked leaves (a tuple).

        This is read-only because you cannot change the chunk size of a
        leaf once it has been created.
        """
        return getattr(self, '_v_chunkshape', None)

    @property
    def object_id(self):
        """A node identifier, which may change from run to run.
        (This is an easier-to-write alias of :attr:`Node._v_objectid`).

        .. versionchanged:: 3.0
           The *objectID* property has been renamed into *object_id*.

        """
        return self._v_objectid

    @property
    def ndim(self):
        """The number of dimensions of the leaf data.

        .. versionadded: 2.4"""
        return len(self.shape)

    @lazyattr
    def filters(self):
        """Filter properties for this leaf.

        See Also
        --------
        Filters

        """

        return Filters._from_leaf(self)

    @property
    def track_times(self):
        """Whether timestamps for the leaf are recorded

        If the leaf is not a dataset, this will fail with HDF5ExtError.

        The track times dataset creation property does not seem to
        survive closing and reopening as of HDF5 1.8.17.  Currently,
        it may be more accurate to test whether the ctime for the
        dataset is 0:
        track_times = (leaf._get_obj_timestamps().ctime == 0)
        """
        return self._get_obj_track_times()

    @property
    def maindim(self):
        """The dimension along which iterators work.

        Its value is 0 (i.e. the first dimension) when the dataset is not
        extendable, and self.extdim (where available) for extendable ones.
        """

        if self.extdim < 0:
            return 0  # choose the first dimension
        return self.extdim

    @property
    def flavor(self):
        """The type of data object read from this leaf.

        It can be any of 'numpy' or 'python'.

        You can (and are encouraged to) use this property to get, set
        and delete the FLAVOR HDF5 attribute of the leaf. When the leaf
        has no such attribute, the default flavor is used..
        """

        return self._flavor

    @flavor.setter
    def flavor(self, flavor):
        self._v_file._check_writable()
        check_flavor(flavor)
        self._v_attrs.FLAVOR = self._flavor = flavor  # logs the change

    @flavor.deleter
    def flavor(self):
        del self._v_attrs.FLAVOR
        self._flavor = internal_flavor

    @property
    def size_on_disk(self):
        """
        The size of this leaf's data in bytes as it is stored on disk.  If the
        data is compressed, this shows the compressed size.  In the case of
        uncompressed, chunked data, this may be slightly larger than the amount
        of data, due to partially filled chunks.
        """
        return self._get_storage_size()

    def __init__(self, parentnode, name,
                 new=False, filters=None,
                 byteorder=None, _log=True,
                 track_times=True):
        self._v_new = new
        """Is this the first time the node has been created?"""
        self.nrowsinbuf = None
        """
        The number of rows that fits in internal input buffers.

        You can change this to fine-tune the speed or memory
        requirements of your application.
        """
        self._flavor = None
        """Private storage for the `flavor` property."""

        if new:
            # Get filter properties from parent group if not given.
            if filters is None:
                filters = parentnode._v_filters
            self.__dict__['filters'] = filters  # bypass the property

            if byteorder not in (None, 'little', 'big'):
                raise ValueError(
                    "the byteorder can only take 'little' or 'big' values "
                    "and you passed: %s" % byteorder)
            self.byteorder = byteorder
            """The byte ordering of the leaf data *on disk*."""

        self._want_track_times = track_times

        # Existing filters need not be read since `filters`
        # is a lazy property that automatically handles their loading.

        super().__init__(parentnode, name, _log)

    def __len__(self):
        """Return the length of the main dimension of the leaf data.

        Please note that this may raise an OverflowError on 32-bit platforms
        for datasets having more than 2**31-1 rows.  This is a limitation of
        Python that you can work around by using the nrows or shape attributes.

        """

        return self.nrows

    def __str__(self):
        """The string representation for this object is its pathname in the
        HDF5 object tree plus some additional metainfo."""

        filters = []
        if self.filters.fletcher32:
            filters.append("fletcher32")
        if self.filters.complevel:
            if self.filters.shuffle:
                filters.append("shuffle")
            if self.filters.bitshuffle:
                filters.append("bitshuffle")
            filters.append(f"{self.filters.complib}({self.filters.complevel})")
        return (f"{self._v_pathname} ({self.__class__.__name__}"
                f"{self.shape}{', '.join(filters)}) {self._v_title!r}")

    def _g_post_init_hook(self):
        """Code to be run after node creation and before creation logging.

        This method gets or sets the flavor of the leaf.

        """

        super()._g_post_init_hook()
        if self._v_new:  # set flavor of new node
            if self._flavor is None:
                self._flavor = internal_flavor
            else:  # flavor set at creation time, do not log
                if self._v_file.params['PYTABLES_SYS_ATTRS']:
                    self._v_attrs._g__setattr('FLAVOR', self._flavor)
        else:  # get flavor of existing node (if any)
            if self._v_file.params['PYTABLES_SYS_ATTRS']:
                flavor = getattr(self._v_attrs, 'FLAVOR', internal_flavor)
                self._flavor = flavor_alias_map.get(flavor, flavor)
            else:
                self._flavor = internal_flavor

    def _calc_chunkshape(self, expectedrows, rowsize, itemsize):
        """Calculate the shape for the HDF5 chunk."""

        # In case of a scalar shape, return the unit chunksize
        if self.shape == ():
            return (SizeType(1),)

        # Compute the chunksize
        MB = 1024 * 1024
        expected_mb = (expectedrows * rowsize) // MB
        chunksize = calc_chunksize(expected_mb)
        complib = self.filters.complib
        if (complib is not None and
            complib.startswith("blosc2") and
            self._c_classid in ('TABLE', 'CARRAY', 'EARRAY')):
            # Blosc2 can introspect into blocks, so we can increase the
            # chunksize for improving HDF5 perf for its internal btree.
            # For the time being, this has been implemented efficiently
            # just for tables, but in the future *Array objects could also
            # be included.
            # Use a decent default value for chunksize
            chunksize *= 16
            # Now, go explore the L3 size and try to find a smarter chunksize
            if 'l3_cache_size' in cpu_info:
                # In general, is a good idea to set the chunksize equal to L3
                l3_cache_size = cpu_info['l3_cache_size']
                # cpuinfo sometimes returns cache sizes as strings (like,
                # "4096 KB"), so refuse the temptation to guess and use the
                # value only when it is an actual int.
                # Also, sometimes cpuinfo does not return a correct L3 size;
                # so in general, enforcing L3 > L2 is a good sanity check.
                l2_cache_size = cpu_info.get('l2_cache_size', "Not found")
                if (type(l3_cache_size) is int and
                    type(l2_cache_size) is int and
                    l3_cache_size > l2_cache_size):
                    chunksize = l3_cache_size
            # In Blosc2, the chunksize cannot be larger than 2 GB - BLOSC2_MAX_BUFFERSIZE
            if chunksize > 2**31 - 32:
                chunksize = 2**31 - 32

        maindim = self.maindim
        # Compute the chunknitems
        chunknitems = chunksize // itemsize
        # Safeguard against itemsizes being extremely large
        if chunknitems == 0:
            chunknitems = 1
        chunkshape = list(self.shape)
        # Check whether trimming the main dimension is enough
        chunkshape[maindim] = 1
        newchunknitems = np.prod(chunkshape, dtype=SizeType)
        if newchunknitems <= chunknitems:
            chunkshape[maindim] = chunknitems // newchunknitems
        else:
            # No, so start trimming other dimensions as well
            for j in range(len(chunkshape)):
                # Check whether trimming this dimension is enough
                chunkshape[j] = 1
                newchunknitems = np.prod(chunkshape, dtype=SizeType)
                if newchunknitems <= chunknitems:
                    chunkshape[j] = chunknitems // newchunknitems
                    break
            else:
                # Ops, we ran out of the loop without a break
                # Set the last dimension to chunknitems
                chunkshape[-1] = chunknitems

        return tuple(SizeType(s) for s in chunkshape)

    def _calc_nrowsinbuf(self):
        """Calculate the number of rows that fits on a PyTables buffer."""

        params = self._v_file.params
        # Compute the nrowsinbuf
        rowsize = self.rowsize
        buffersize = params['IO_BUFFER_SIZE']
        if rowsize != 0:
            nrowsinbuf = buffersize // rowsize
        else:
            nrowsinbuf = 1

        # Safeguard against row sizes being extremely large
        if nrowsinbuf == 0:
            nrowsinbuf = 1
            # If rowsize is too large, issue a Performance warning
            maxrowsize = params['BUFFER_TIMES'] * buffersize
            if rowsize > maxrowsize:
                warnings.warn("""\
The Leaf ``%s`` is exceeding the maximum recommended rowsize (%d bytes);
be ready to see PyTables asking for *lots* of memory and possibly slow
I/O.  You may want to reduce the rowsize by trimming the value of
dimensions that are orthogonal (and preferably close) to the *main*
dimension of this leave.  Alternatively, in case you have specified a
very small/large chunksize, you may want to increase/decrease it."""
                              % (self._v_pathname, maxrowsize),
                              PerformanceWarning)
        return nrowsinbuf

    # This method is appropriate for calls to __getitem__ methods
    def _process_range(self, start, stop, step, dim=None, warn_negstep=True):
        if dim is None:
            nrows = self.nrows  # self.shape[self.maindim]
        else:
            nrows = self.shape[dim]

        if warn_negstep and step and step < 0:
            raise ValueError("slice step cannot be negative")

        # if start is not None: start = long(start)
        # if stop is not None: stop = long(stop)
        # if step is not None: step = long(step)

        return slice(start, stop, step).indices(int(nrows))

    # This method is appropriate for calls to read() methods
    def _process_range_read(self, start, stop, step, warn_negstep=True):
        nrows = self.nrows
        if start is not None and stop is None and step is None:
            # Protection against start greater than available records
            # nrows == 0 is a special case for empty objects
            if 0 < nrows <= start:
                raise IndexError("start of range (%s) is greater than "
                                 "number of rows (%s)" % (start, nrows))
            step = 1
            if start == -1:  # corner case
                stop = nrows
            else:
                stop = start + 1
        # Finally, get the correct values (over the main dimension)
        start, stop, step = self._process_range(start, stop, step,
                                                warn_negstep=warn_negstep)
        return (start, stop, step)

    def _g_copy(self, newparent, newname, recursive, _log=True, **kwargs):
        # Compute default arguments.
        start = kwargs.pop('start', None)
        stop = kwargs.pop('stop', None)
        step = kwargs.pop('step', None)
        title = kwargs.pop('title', self._v_title)
        filters = kwargs.pop('filters', self.filters)
        chunkshape = kwargs.pop('chunkshape', self.chunkshape)
        copyuserattrs = kwargs.pop('copyuserattrs', True)
        stats = kwargs.pop('stats', None)
        if chunkshape == 'keep':
            chunkshape = self.chunkshape  # Keep the original chunkshape
        elif chunkshape == 'auto':
            chunkshape = None             # Will recompute chunkshape

        # Fix arguments with explicit None values for backwards compatibility.
        if title is None:
            title = self._v_title
        if filters is None:
            filters = self.filters

        # Create a copy of the object.
        (new_node, bytes) = self._g_copy_with_stats(
            newparent, newname, start, stop, step,
            title, filters, chunkshape, _log, **kwargs)

        # Copy user attributes if requested (or the flavor at least).
        if copyuserattrs:
            self._v_attrs._g_copy(new_node._v_attrs, copyclass=True)
        elif 'FLAVOR' in self._v_attrs:
            if self._v_file.params['PYTABLES_SYS_ATTRS']:
                new_node._v_attrs._g__setattr('FLAVOR', self._flavor)
        new_node._flavor = self._flavor  # update cached value

        # Update statistics if needed.
        if stats is not None:
            stats['leaves'] += 1
            stats['bytes'] += bytes

        return new_node

    def _g_fix_byteorder_data(self, data, dbyteorder):
        """Fix the byteorder of data passed in constructors."""
        dbyteorder = byteorders[dbyteorder]
        # If self.byteorder has not been passed as an argument of
        # the constructor, then set it to the same value of data.
        if self.byteorder is None:
            self.byteorder = dbyteorder
        # Do an additional in-place byteswap of data if the in-memory
        # byteorder doesn't match that of the on-disk.  This is the only
        # place that we have to do the conversion manually. In all the
        # other cases, it will be HDF5 the responsible for doing the
        # byteswap properly.
        if dbyteorder in ['little', 'big']:
            if dbyteorder != self.byteorder:
                # if data is not writeable, do a copy first
                if not data.flags.writeable:
                    data = data.copy()
                data.byteswap(True)
        else:
            # Fix the byteorder again, no matter which byteorder have
            # specified the user in the constructor.
            self.byteorder = "irrelevant"
        return data

    def _point_selection(self, key):
        """Perform a point-wise selection.

        `key` can be any of the following items:

        * A boolean array with the same shape than self. Those positions
          with True values will signal the coordinates to be returned.

        * A numpy array (or list or tuple) with the point coordinates.
          This has to be a two-dimensional array of size len(self.shape)
          by num_elements containing a list of zero-based values
          specifying the coordinates in the dataset of the selected
          elements. The order of the element coordinates in the array
          specifies the order in which the array elements are iterated
          through when I/O is performed. Duplicate coordinate locations
          are not checked for.

        Return the coordinates array.  If this is not possible, raise a
        `TypeError` so that the next selection method can be tried out.

        This is useful for whatever `Leaf` instance implementing a
        point-wise selection.

        """
        input_key = key
        if type(key) in (list, tuple):
            if isinstance(key, tuple) and len(key) > len(self.shape):
                raise IndexError(f"Invalid index or slice: {key!r}")
            # Try to convert key to a numpy array.  If not possible,
            # a TypeError will be issued (to be catched later on).
            try:
                key = toarray(key)
            except ValueError:
                raise TypeError(f"Invalid index or slice: {key!r}")
        elif not isinstance(key, np.ndarray):
            raise TypeError(f"Invalid index or slice: {key!r}")

        # Protection against empty keys
        if len(key) == 0:
            return np.array([], dtype="i8")

        if key.dtype.kind == 'b':
            if not key.shape == self.shape:
                raise IndexError(
                    "Boolean indexing array has incompatible shape")
            # Get the True coordinates (64-bit indices!)
            coords = np.asarray(key.nonzero(), dtype='i8')
            coords = np.transpose(coords)
        elif key.dtype.kind == 'i' or key.dtype.kind == 'u':
            if len(key.shape) > 2:
                raise IndexError(
                    "Coordinate indexing array has incompatible shape")
            elif len(key.shape) == 2:
                if key.shape[0] != len(self.shape):
                    raise IndexError(
                        "Coordinate indexing array has incompatible shape")
                coords = np.asarray(key, dtype="i8")
                coords = np.transpose(coords)
            else:
                # For 1-dimensional datasets
                coords = np.asarray(key, dtype="i8")

            # handle negative indices
            base = coords if coords.base is None else coords.base
            if base is input_key:
                # never modify the original "key" data
                coords = coords.copy()

            idx = coords < 0
            coords[idx] = (coords + self.shape)[idx]

            # bounds check
            if np.any(coords < 0) or np.any(coords >= self.shape):
                raise IndexError("Index out of bounds")
        else:
            raise TypeError("Only integer coordinates allowed.")
        # We absolutely need a contiguous array
        if not coords.flags.contiguous:
            coords = coords.copy()
        return coords

    # Tree manipulation
    def remove(self):
        """Remove this node from the hierarchy.

        This method has the behavior described
        in :meth:`Node._f_remove`. Please note that there is no recursive flag
        since leaves do not have child nodes.

        """

        self._f_remove(False)

    def rename(self, newname):
        """Rename this node in place.

        This method has the behavior described in :meth:`Node._f_rename()`.

        """

        self._f_rename(newname)

    def move(self, newparent=None, newname=None,
             overwrite=False, createparents=False):
        """Move or rename this node.

        This method has the behavior described in :meth:`Node._f_move`

        """

        self._f_move(newparent, newname, overwrite, createparents)

    def copy(self, newparent=None, newname=None,
             overwrite=False, createparents=False, **kwargs):
        """Copy this node and return the new one.

        This method has the behavior described in :meth:`Node._f_copy`. Please
        note that there is no recursive flag since leaves do not have child
        nodes.

        .. warning::

            Note that unknown parameters passed to this method will be
            ignored, so may want to double check the spelling of these
            (i.e. if you write them incorrectly, they will most probably
            be ignored).

        Parameters
        ----------
        title
            The new title for the destination. If omitted or None, the original
            title is used.
        filters : Filters
            Specifying this parameter overrides the original filter properties
            in the source node. If specified, it must be an instance of the
            Filters class (see :ref:`FiltersClassDescr`). The default is to
            copy the filter properties from the source node.
        copyuserattrs
            You can prevent the user attributes from being copied by setting
            this parameter to False. The default is to copy them.
        start, stop, step : int
            Specify the range of rows to be copied; the default is to copy all
            the rows.
        stats
            This argument may be used to collect statistics on the copy
            process. When used, it should be a dictionary with keys 'groups',
            'leaves' and 'bytes' having a numeric value. Their values will be
            incremented to reflect the number of groups, leaves and bytes,
            respectively, that have been copied during the operation.
        chunkshape
            The chunkshape of the new leaf.  It supports a couple of special
            values.  A value of keep means that the chunkshape will be the same
            than original leaf (this is the default).  A value of auto means
            that a new shape will be computed automatically in order to ensure
            best performance when accessing the dataset through the main
            dimension.  Any other value should be an integer or a tuple
            matching the dimensions of the leaf.

        """

        return self._f_copy(
            newparent, newname, overwrite, createparents, **kwargs)

    def truncate(self, size):
        """Truncate the main dimension to be size rows.

        If the main dimension previously was larger than this size, the extra
        data is lost.  If the main dimension previously was shorter, it is
        extended, and the extended part is filled with the default values.

        The truncation operation can only be applied to *enlargeable* datasets,
        else a TypeError will be raised.

        """

        # A non-enlargeable arrays (Array, CArray) cannot be truncated
        if self.extdim < 0:
            raise TypeError("non-enlargeable datasets cannot be truncated")
        self._g_truncate(size)

    def isvisible(self):
        """Is this node visible?

        This method has the behavior described in :meth:`Node._f_isvisible()`.

        """

        return self._f_isvisible()

    # Attribute handling
    def get_attr(self, name):
        """Get a PyTables attribute from this node.

        This method has the behavior described in :meth:`Node._f_getattr`.

        """

        return self._f_getattr(name)

    def set_attr(self, name, value):
        """Set a PyTables attribute for this node.

        This method has the behavior described in :meth:`Node._f_setattr()`.

        """

        self._f_setattr(name, value)

    def del_attr(self, name):
        """Delete a PyTables attribute from this node.

        This method has the behavior described in :meth:`Node_f_delAttr`.

        """

        self._f_delattr(name)

    # Data handling
    def flush(self):
        """Flush pending data to disk.

        Saves whatever remaining buffered data to disk. It also releases
        I/O buffers, so if you are filling many datasets in the same
        PyTables session, please call flush() extensively so as to help
        PyTables to keep memory requirements low.

        """

        self._g_flush()

    def _f_close(self, flush=True):
        """Close this node in the tree.

        This method has the behavior described in :meth:`Node._f_close`.
        Besides that, the optional argument flush tells whether to flush
        pending data to disk or not before closing.

        """

        if not self._v_isopen:
            return  # the node is already closed or not initialized

        # Only do a flush in case the leaf has an IO buffer.  The
        # internal buffers of HDF5 will be flushed afterwards during the
        # self._g_close() call.  Avoiding an unnecessary flush()
        # operation accelerates the closing for the unbuffered leaves.
        if flush and hasattr(self, "_v_iobuf"):
            self.flush()

        # Close the dataset and release resources
        self._g_close()

        # Close myself as a node.
        super()._f_close()

    def close(self, flush=True):
        """Close this node in the tree.

        This method is completely equivalent to :meth:`Leaf._f_close`.

        """

        self._f_close(flush)
