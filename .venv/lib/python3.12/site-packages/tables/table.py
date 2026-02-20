"""Here is defined the Table class."""

from __future__ import annotations

import sys
import math
import weakref
import operator
import warnings
import functools
from time import perf_counter as clock
from typing import Any, Literal, TYPE_CHECKING
from pathlib import Path
from collections.abc import Callable, Generator, Iterator, Sequence

import numpy as np
import numexpr as ne
import numpy.typing as npt

from . import tableextension
from .atom import Atom
from .leaf import Leaf
from .path import join_path, split_path
from .index import (
    OldIndex,
    default_index_filters,
    default_auto_index,
    Index,
    IndexesDescG,
    IndexesTableG,
)
from .utils import is_idx, lazyattr, SizeType
from .utils import NailedDict as CacheDict
from .flavor import flavor_of, array_as_internal, internal_to_flavor
from .conditions import compile_condition
from .exceptions import (
    NodeError,
    HDF5ExtError,
    PerformanceWarning,
    OldIndexWarning,
    NoSuchNodeError,
)
from .description import IsDescription, Description, Col, descr_from_dtype
from .utilsextension import get_nested_field
from .lrucacheextension import ObjectCache, NumCache

profile = False
# profile = True  # Uncomment for profiling
if profile:
    from .utils import show_stats

if TYPE_CHECKING:
    from .node import Node
    from .group import Group
    from .filters import Filters
    from .misc.enum import Enum
    from .conditions import CompiledCondition

# 2.2: Added support for complex types. Introduced in version 0.9.
# 2.2.1: Added support for time types.
# 2.3: Changed the indexes naming schema.
# 2.4: Changed indexes naming schema (again).
# 2.5: Added the FIELD_%d_FILL attributes.
# 2.6: Added the FLAVOR attribute (optional).
# 2.7: Numeric and numarray flavors are gone.
obversion = "2.7"  # The Table VERSION number


# Maps NumPy types to the types used by Numexpr.
_nxtype_from_nptype = {
    np.bool_: bool,
    np.int8: ne.necompiler.int_,
    np.int16: ne.necompiler.int_,
    np.int32: ne.necompiler.int_,
    np.int64: ne.necompiler.long_,
    np.uint8: ne.necompiler.int_,
    np.uint16: ne.necompiler.int_,
    np.uint32: ne.necompiler.long_,
    np.uint64: ne.necompiler.long_,
    np.float32: float,
    np.float64: ne.necompiler.double,
    np.complex64: complex,
    np.complex128: complex,
    np.bytes_: bytes,
}

_nxtype_from_nptype[np.str_] = str

if hasattr(np, "float16"):
    _nxtype_from_nptype[np.float16] = float  # XXX: check
if hasattr(np, "float96"):
    _nxtype_from_nptype[np.float96] = ne.necompiler.double  # XXX: check
if hasattr(np, "float128"):
    _nxtype_from_nptype[np.float128] = ne.necompiler.double  # XXX: check
if hasattr(np, "complex192"):
    _nxtype_from_nptype[np.complex192] = complex  # XXX: check
if hasattr(np, "complex256"):
    _nxtype_from_nptype[np.complex256] = complex  # XXX: check


# The NumPy scalar type corresponding to `SizeType`.
_npsizetype = np.array(SizeType(0)).dtype.type


def _index_name_of(node: Node) -> str:
    return "_i_%s" % node._v_name


def _index_pathname_of(node: Node) -> str:
    node_parent_path = split_path(node._v_pathname)[0]
    return join_path(node_parent_path, _index_name_of(node))


def _index_pathname_of_column(table: Table, colpathname: str) -> str:
    return join_path(_index_pathname_of(table), colpathname)


# The next are versions that work with just paths (i.e. we don't need
# a node instance for using them, which can be critical in certain
# situations)


def _index_name_of_(nodeName: str) -> str:  # noqa: N803
    return "_i_%s" % nodeName


def _index_pathname_of_(nodePath: str) -> str:  # noqa: N803
    node_parent_path, node_name = split_path(nodePath)
    return join_path(node_parent_path, _index_name_of_(node_name))


def _index_pathname_of_column_(table_path: str, colpathname: str) -> str:
    return join_path(_index_pathname_of_(table_path), colpathname)


def restorecache(self: Table) -> None:
    """Restore the chunk cache."""
    # Define a cache for sparse table reads
    params = self._v_file.params
    chunksize = self._v_chunkshape[0]
    nslots = params["TABLE_MAX_SIZE"] / (chunksize * self._v_dtype.itemsize)
    self._chunkcache = NumCache(
        (nslots, chunksize), self._v_dtype, "table chunk cache"
    )
    self._seqcache = ObjectCache(
        params["ITERSEQ_MAX_SLOTS"],
        params["ITERSEQ_MAX_SIZE"],
        "Iter sequence cache",
    )
    self._dirtycache = False


def _table__where_indexed(
    self: Table,
    compiled: CompiledCondition,
    condition: str,
    condvars: dict[str, Column | np.ndarray],
    start: int,
    stop: int,
    step: int,
) -> Iterator[tableextension.Row] | np.ndarray:
    if profile:
        tref = clock()
    if profile:
        show_stats("Entering table_whereIndexed", tref)
    self._use_index = True
    # Clean the table caches for indexed queries if needed
    if self._dirtycache:
        restorecache(self)

    # Get the values in expression that are not columns
    values = []
    for key, value in condvars.items():
        if isinstance(value, np.ndarray):
            values.append((key, value.item()))
    # Build a key for the sequence cache
    seqkey = (condition, tuple(values), (start, stop, step))
    # Do a lookup in sequential cache for this query
    nslot = self._seqcache.getslot(seqkey)
    if nslot >= 0:
        # Get the row sequence from the cache
        seq = self._seqcache.getitem(nslot)
        if len(seq) == 0:
            return iter([])
        # seq is a list.
        seq = np.array(seq, dtype="int64")
        # Correct the ranges in cached sequence
        if (start, stop, step) != (0, self.nrows, 1):
            seq = seq[
                (seq >= start) & (seq < stop) & ((seq - start) % step == 0)
            ]
        return self.itersequence(seq)
    else:
        # No luck.  self._seqcache will be populated
        # in the iterator if possible. (Row._finish_riterator)
        self._seqcache_key = seqkey

    # Compute the chunkmap for every index in indexed expression
    idxexprs = compiled.index_expressions
    strexpr = compiled.string_expression
    cmvars = {}
    tcoords = 0
    for i, idxexpr in enumerate(idxexprs):
        var, ops, lims = idxexpr
        col = condvars[var]
        index = col.index
        assert index is not None, "the chosen column is not indexed"
        assert not index.dirty, "the chosen column has a dirty index"

        # Get the number of rows that the indexed condition yields.
        range_ = index.get_lookup_range(ops, lims)
        ncoords = index.search(range_)
        tcoords += ncoords
        if index.reduction == 1 and ncoords == 0:
            # No values from index condition, thus the chunkmap should be empty
            nrowsinchunk = self.chunkshape[0]
            nchunks = math.ceil(self.nrows / nrowsinchunk)
            chunkmap = np.zeros(shape=nchunks, dtype="bool")
        else:
            # Get the chunkmap from the index
            chunkmap = index.get_chunkmap()
        # Assign the chunkmap to the cmvars dictionary
        cmvars["e%d" % i] = chunkmap

    if index.reduction == 1 and tcoords == 0:
        # No candidates found in any indexed expression component, so leave now
        self._seqcache.setitem(seqkey, [], 1)
        return iter([])

    # Compute the final chunkmap
    chunkmap = ne.evaluate(strexpr, cmvars)
    if not chunkmap.any():
        # The chunkmap is all False, so the result is empty
        self._seqcache.setitem(seqkey, [], 1)
        return iter([])

    if profile:
        show_stats("Exiting table_whereIndexed", tref)
    return chunkmap


def create_indexes_table(table: Table) -> IndexesTableG:
    """Create indexes for a table."""
    itgroup = IndexesTableG(
        table._v_parent,
        _index_name_of(table),
        "Indexes container for table " + table._v_pathname,
        new=True,
    )
    return itgroup


def create_indexes_descr(
    igroup: Group, dname: str, iname: str, filters: Filters | None
) -> IndexesDescG:
    """Create indexes descriptor."""
    idgroup = IndexesDescG(
        igroup,
        iname,
        "Indexes container for sub-description " + dname,
        filters=filters,
        new=True,
    )
    return idgroup


def _column__create_index(
    self: Column,
    optlevel: int,
    kind: str,
    filters: Filters | None,
    tmp_dir: str,
    blocksizes: tuple[int, int, int, int],
    verbose: bool,
) -> int:
    name = self.name
    table = self.table
    dtype = self.dtype
    descr = self.descr
    index = self.index
    get_node = table._v_file._get_node

    # Warn if the index already exists
    if index:
        raise ValueError(
            "%s for column '%s' already exists. If you want to "
            "re-create it, please, try with reindex() method "
            "better" % (str(index), str(self.pathname))
        )

    # Check that the datatype is indexable.
    if dtype.str[1:] == "u8":
        raise NotImplementedError(
            "indexing 64-bit unsigned integer columns "
            "is not supported yet, sorry"
        )
    if dtype.kind == "c":
        raise TypeError("complex columns can not be indexed")
    if dtype.shape != ():
        raise TypeError("multidimensional columns can not be indexed")

    # Get the indexes group for table, and if not exists, create it
    try:
        itgroup = get_node(_index_pathname_of(table))
    except NoSuchNodeError:
        itgroup = create_indexes_table(table)

    # Create the necessary intermediate groups for descriptors
    idgroup = itgroup
    dname = ""
    pathname = descr._v_pathname
    if pathname != "":
        inames = pathname.split("/")
        for iname in inames:
            if dname == "":
                dname = iname
            else:
                dname += "/" + iname
            try:
                idgroup = get_node(f"{itgroup._v_pathname}/{dname}")
            except NoSuchNodeError:
                idgroup = create_indexes_descr(idgroup, dname, iname, filters)

    # Create the atom
    assert dtype.shape == ()
    atom = Atom.from_dtype(np.dtype((dtype, (0,))))

    # Protection on tables larger than the expected rows (perhaps the
    # user forgot to pass this parameter to the Table constructor?)
    expectedrows = table._v_expectedrows
    if table.nrows > expectedrows:
        expectedrows = table.nrows

    # Create the index itself
    index = Index(
        idgroup,
        name,
        atom=atom,
        title="Index for %s column" % name,
        kind=kind,
        optlevel=optlevel,
        filters=filters,
        tmp_dir=tmp_dir,
        expectedrows=expectedrows,
        byteorder=table.byteorder,
        blocksizes=blocksizes,
    )

    table._set_column_indexing(self.pathname, True)

    # Feed the index with values

    # Add rows to the index if necessary
    if table.nrows > 0:
        indexedrows = table._add_rows_to_index(
            self.pathname, 0, table.nrows, lastrow=True, update=False
        )
    else:
        indexedrows = 0
    index.dirty = False
    table._indexedrows = indexedrows
    table._unsaved_indexedrows = table.nrows - indexedrows

    # Optimize the index that has been already filled-up
    index.optimize(verbose=verbose)

    # We cannot do a flush here because when reindexing during a
    # flush, the indexes are created anew, and that creates a nested
    # call to flush().
    # table.flush()

    return indexedrows


class _ColIndexes(dict):
    """Provide a nice representation of column indexes."""

    def __repr__(self) -> str:
        """Return a detailed Description column representation."""
        rep = [f'  "{k}": {v}' for k, v in self.items()]
        return "{\n  %s}" % (",\n  ".join(rep))


class Table(tableextension.Table, Leaf):
    """This class represents heterogeneous datasets in an HDF5 file.

    Tables are leaves (see the Leaf class in :ref:`LeafClassDescr`) whose data
    consists of a unidimensional sequence of *rows*, where each row contains
    one or more *fields*.  Fields have an associated unique *name* and
    *position*, with the first field having position 0.  All rows have the same
    fields, which are arranged in *columns*.

    Fields can have any type supported by the Col class (see
    :ref:`ColClassDescr`) and its descendants, which support multidimensional
    data.  Moreover, a field can be *nested* (to an arbitrary depth), meaning
    that it includes further fields inside.  A field named x inside a nested
    field a in a table can be accessed as the field a/x (its *path name*) from
    the table.

    The structure of a table is declared by its description, which is made
    available in the Table.description attribute (see :class:`Table`).

    This class provides new methods to read, write and search table data
    efficiently.  It also provides special Python methods to allow accessing
    the table as a normal sequence or array (with extended slicing supported).

    PyTables supports *in-kernel* searches working simultaneously on several
    columns using complex conditions.  These are faster than selections using
    Python expressions.  See the :meth:`Table.where` method for more
    information on in-kernel searches.

    Non-nested columns can be *indexed*.  Searching an indexed column can be
    several times faster than searching a non-nested one.  Search methods
    automatically take advantage of indexing where available.

    When iterating a table, an object from the Row (see :ref:`RowClassDescr`)
    class is used.  This object allows to read and write data one row at a
    time, as well as to perform queries which are not supported by in-kernel
    syntax (at a much lower speed, of course).

    Objects of this class support access to individual columns via *natural
    naming* through the :attr:`Table.cols` accessor.  Nested columns are
    mapped to Cols instances, and non-nested ones to Column instances.
    See the Column class in :ref:`ColumnClassDescr` for examples of this
    feature.

    Parameters
    ----------
    parentnode
        The parent :class:`Group` object.

        .. versionchanged:: 3.0
           Renamed from *parentNode* to *parentnode*.

    name : str
        The name of this node in its parent group.
    description
        An IsDescription subclass or a dictionary where the keys are the field
        names, and the values the type definitions. In addition, a pure NumPy
        dtype is accepted.  If None, the table metadata is read from disk,
        else, it's taken from previous parameters.
    title
        Sets a TITLE attribute on the HDF5 table entity.
    filters : Filters
        An instance of the Filters class that provides information about the
        desired I/O filters to be applied during the life of this object.
    expectedrows
        A user estimate about the number of rows that will be on table. If not
        provided, the default value is ``EXPECTED_ROWS_TABLE`` (see
        ``tables/parameters.py``).  If you plan to save bigger tables, try
        providing a guess; this will optimize the HDF5 B-Tree creation and
        management process time and memory used.
    chunkshape
        The shape of the data chunk to be read or written as a single HDF5 I/O
        operation. The filters are applied to those chunks of data. Its rank
        for tables has to be 1.  If ``None``, a sensible value is calculated
        based on the `expectedrows` parameter (which is recommended).
    byteorder
        The byteorder of the data *on-disk*, specified as 'little' or 'big'. If
        this is not specified, the byteorder is that of the platform, unless
        you passed a recarray as the `description`, in which case the recarray
        byteorder will be chosen.
    track_times
        Whether time data associated with the leaf are recorded (object
        access time, raw data modification time, metadata change time, object
        birth time); default True.  Semantics of these times depend on their
        implementation in the HDF5 library: refer to documentation of the
        H5O_info_t data structure.  As of HDF5 1.8.15, only ctime (metadata
        change time) is implemented.

        .. versionadded:: 3.4.3

    Notes
    -----
    The instance variables below are provided in addition to those in
    Leaf (see :ref:`LeafClassDescr`).  Please note that there are several
    col* dictionaries to ease retrieving information about a column
    directly by its path name, avoiding the need to walk through
    Table.description or Table.cols.


    .. rubric:: Table attributes

    .. attribute:: coldescrs

        Maps the name of a column to its Col description (see
        :ref:`ColClassDescr`).

    .. attribute:: coldflts

        Maps the name of a column to its default value.

    .. attribute:: coldtypes

        Maps the name of a column to its NumPy data type.

    .. attribute:: colindexed

        Is the column which name is used as a key indexed?

    .. attribute:: colinstances

        Maps the name of a column to its Column (see
        :ref:`ColumnClassDescr`) or Cols (see :ref:`ColsClassDescr`)
        instance.

    .. attribute:: colnames

        A list containing the names of *top-level* columns in the table.

    .. attribute:: colpathnames

        A list containing the pathnames of *bottom-level* columns in
        the table.

        These are the leaf columns obtained when walking the table
        description left-to-right, bottom-first. Columns inside a
        nested column have slashes (/) separating name components in
        their pathname.

    .. attribute:: cols

        A Cols instance that provides *natural naming* access to
        non-nested (Column, see :ref:`ColumnClassDescr`) and nested
        (Cols, see :ref:`ColsClassDescr`) columns.

    .. attribute:: coltypes

        Maps the name of a column to its PyTables data type.

    .. attribute:: description

        A Description instance (see :ref:`DescriptionClassDescr`)
        reflecting the structure of the table.

    .. attribute:: extdim

        The index of the enlargeable dimension (always 0 for tables).

    .. attribute:: indexed

        Does this table have any indexed columns?

    .. attribute:: nrows

        The current number of rows in the table.

    """

    # Class identifier.
    _c_classid = "TABLE"

    @lazyattr
    def row(self) -> tableextension.Row:
        """Row instance (see :ref:`RowClassDescr`) associated to the Table."""
        return tableextension.Row(self)

    @lazyattr
    def dtype(self) -> np.dtype:
        """Numpy ``dtype`` that most closely matches this table."""
        return self.description._v_dtype

    @property
    def shape(self) -> tuple[int]:
        """Shape of this table."""
        return (self.nrows,)

    @property
    def rowsize(self) -> int:
        """Size in bytes of each row in the table."""
        return self.description._v_dtype.itemsize

    @property
    def size_in_memory(self) -> int:
        """Size of the table's data in bytes when fully loaded in memory.

        This may be used in combination with size_on_disk to calculate
        the compression ratio of the data.
        """
        return self.nrows * self.rowsize

    @lazyattr
    def _v_iobuf(self) -> np.ndarray:
        """Buffer for doing I/O."""
        return self._get_container(self.nrowsinbuf)

    @lazyattr
    def _v_wdflts(self) -> np.ndarray | None:
        """Contain the default values for writing in recarray format."""
        # First, do a check to see whether we need to set default values
        # different from 0 or not.
        for coldflt in self.coldflts.values():
            if isinstance(coldflt, np.ndarray) or coldflt:
                break
        else:
            # No default different from 0 found.  Returning None.
            return None
        wdflts = self._get_container(1)
        for colname, coldflt in self.coldflts.items():
            ra = get_nested_field(wdflts, colname)
            ra[:] = coldflt
        return wdflts

    @lazyattr
    def _colunaligned(self) -> frozenset:
        """Pathnames of unaligned, *unidimensional* columns."""
        colunaligned, rarr = [], self._get_container(0)
        for colpathname in self.colpathnames:
            carr = get_nested_field(rarr, colpathname)
            if not carr.flags.aligned and carr.ndim == 1:
                colunaligned.append(colpathname)
        return frozenset(colunaligned)

    # **************** WARNING! ***********************
    # This function can be called during the destruction time of a table
    # so measures have been taken so that it doesn't have to revive
    # another node (which can fool the LRU cache). The solution devised
    # has been to add a cache for autoindex (Table._autoindex), populate
    # it in creation time of the cache (which is a safe period) and then
    # update the cache whenever it changes.
    # This solves the error when running test_indexes.py ManyNodesTestCase.
    # F. Alted 2007-04-20
    # **************************************************

    @property
    def autoindex(self) -> bool:
        """Is True if the `Table` automatically keep column indexes up to date.

        Setting this value states whether existing indexes should be
        automatically updated after an append operation or recomputed
        after an index-invalidating operation (i.e. removal and
        modification of rows).  The default is true.

        This value gets into effect whenever a column is altered.  If you
        don't have automatic indexing activated and you want to do an
        immediate update use `Table.flush_rows_to_index()`; for an immediate
        reindexing of invalidated indexes, use `Table.reindex_dirty()`.

        This value is persistent.

        .. versionchanged:: 3.0
           The *autoIndex* property has been renamed into *autoindex*.
        """
        if self._autoindex is None:
            try:
                indexgroup = self._v_file._get_node(_index_pathname_of(self))
            except NoSuchNodeError:
                self._autoindex = default_auto_index  # update cache
                return self._autoindex
            else:
                self._autoindex = indexgroup.auto  # update cache
                return self._autoindex
        else:
            # The value is in cache, return it
            return self._autoindex

    @autoindex.setter
    def autoindex(self, auto: bool) -> None:
        auto = bool(auto)
        try:
            indexgroup = self._v_file._get_node(_index_pathname_of(self))
        except NoSuchNodeError:
            indexgroup = create_indexes_table(self)
        indexgroup.auto = auto
        # Update the cache in table instance as well
        self._autoindex = auto

    @property
    def indexedcolpathnames(self) -> list[str]:
        """List of pathnames of indexed columns in the table."""
        return [
            _colpname
            for _colpname in self.colpathnames
            if self.colindexed[_colpname]
        ]

    @property
    def colindexes(self) -> _ColIndexes:
        """Return a dictionary with the indexes of the indexed columns."""
        return _ColIndexes(
            (_colpname, self.cols._f_col(_colpname).index)
            for _colpname in self.colpathnames
            if self.colindexed[_colpname]
        )

    @property
    def _dirtyindexes(self) -> bool:
        """Whether some index in table is dirty."""
        return self._condition_cache._nailcount > 0

    def __init__(
        self,
        parentnode: Group,
        name: str,
        description: (
            dict | type[IsDescription] | Description | npt.DTypeLike | None
        ) = None,
        title: str = "",
        filters: Filters | None = None,
        expectedrows: int | None = None,
        chunkshape: int | tuple[int] | None = None,
        byteorder: Literal["little", "big", None] = None,
        _log: bool = True,
        track_times: bool = True,
    ) -> None:

        self._v_new = new = description is not None
        """Is this the first time the node has been created?"""
        self._v_new_title = title
        """New title for this node."""
        self._v_new_filters = filters
        """New filter properties for this node."""
        self.extdim = 0  # Tables only have one dimension currently
        """The index of the enlargeable dimension (always 0 for tables)."""
        self._v_recarray = None
        """A structured array to be stored in the table."""
        self._rabyteorder: str | None = None
        """The computed byteorder of the self._v_recarray."""
        if expectedrows is None:
            expectedrows = parentnode._v_file.params["EXPECTED_ROWS_TABLE"]
        self._v_expectedrows: int = expectedrows
        """The expected number of rows to be stored in the table."""
        self.nrows = SizeType(0)
        """The current number of rows in the table."""
        self.description: Description | None = None
        """A Description instance (see :ref:`DescriptionClassDescr`)
        reflecting the structure of the table."""
        self._time64colnames: list[str] = []
        """The names of ``Time64`` columns."""
        self._strcolnames: list[str] = []
        """The names of ``String`` columns."""
        self._colenums: dict[str, Enum] = {}
        """Maps the name of an enumerated column to its ``Enum`` instance."""
        self._v_chunkshape: tuple[int] | None = None
        """Private storage for the `chunkshape` property of the leaf."""

        self.indexed = False
        """Does this table have any indexed columns?"""
        self._indexedrows = 0
        """Number of rows indexed in disk."""
        self._unsaved_indexedrows = 0
        """Number of rows indexed in memory but still not in disk."""
        self._listoldindexes: list[str] = []
        """The list of columns with old indexes."""
        self._autoindex: bool | None = None
        """Private variable that caches the value for autoindex."""

        self.colnames: list[str] = []
        """A list containing the names of *top-level* columns in the table."""
        self.colpathnames: list[str] = []
        """A list containing the pathnames of *bottom-level* columns in the
        table.

        These are the leaf columns obtained when walking the
        table description left-to-right, bottom-first.  Columns inside a
        nested column have slashes (/) separating name components in
        their pathname.
        """
        self.colinstances: dict[str, Column | Cols] = {}
        """Maps the name of a column to its Column (see
        :ref:`ColumnClassDescr`) or Cols (see :ref:`ColsClassDescr`)
        instance."""
        self.coldescrs: dict[str, Col] = {}
        """Maps the name of a column to its Col description (see
        :ref:`ColClassDescr`)."""
        self.coltypes: dict[str, str] = {}
        """Maps the name of a column to its PyTables data type."""
        self.coldtypes: dict[str, np.dtype] = {}
        """Maps the name of a column to its NumPy data type."""
        self.coldflts: dict[str, Any] = {}
        """Maps the name of a column to its default value."""
        self.colindexed: dict[str, bool] = {}
        """Is the column which name is used as a key indexed?"""

        self._use_index = False
        """Whether an index can be used or not in a search.  Boolean."""
        self._where_condition: tuple[Callable, Any, Any] | None = None
        """Condition function and argument list for selection of values."""
        self._seqcache_key = None
        """The key under which to save a query's results (list of row indexes)
        or None to not save."""
        max_slots = parentnode._v_file.params["COND_CACHE_SLOTS"]
        self._condition_cache = CacheDict(max_slots)
        """Cache of already compiled conditions."""
        self._exprvars_cache: dict[str, list[str]] = {}
        """Cache of variables participating in numexpr expressions."""
        self._enabled_indexing_in_queries = True
        """Is indexing enabled in queries?  *Use only for testing.*"""
        self._empty_array_cache: dict[np.dtype, np.ndarray] = {}
        """Cache of empty arrays."""

        self._v_dtype: np.dtype | None = None
        """The NumPy datatype fopr this table."""
        self.cols: Cols | None = None
        """
        A Cols instance that provides *natural naming* access to non-nested
        (Column, see :ref:`ColumnClassDescr`) and nested (Cols, see
        :ref:`ColsClassDescr`) columns.
        """
        self._dirtycache = True
        """Whether the data caches are dirty or not. Initially set to yes."""
        self._descflavor: str | None = None
        """Temporarily keeps the flavor of a description with data."""

        # Initialize this object in case is a new Table

        # Try purely descriptive description objects.
        if new and isinstance(description, dict):
            # Dictionary case
            self.description = Description(
                description, ptparams=parentnode._v_file.params
            )
        elif new and (
            type(description) is type(IsDescription)
            and issubclass(description, IsDescription)
        ):
            # IsDescription subclass case
            descr = description()
            self.description = Description(
                descr.columns, ptparams=parentnode._v_file.params
            )
        elif new and isinstance(description, Description):
            # It is a Description instance already
            self.description = description

        # No description yet?
        if new and self.description is None:
            # Try NumPy dtype instances
            if isinstance(description, np.dtype):
                tup = descr_from_dtype(
                    description, ptparams=parentnode._v_file.params
                )
                self.description, self._rabyteorder = tup

        # No description yet?
        if new and self.description is None:
            # Try structured array description objects.
            try:
                self._descflavor = flavor = flavor_of(description)
            except TypeError:  # probably not an array
                pass
            else:
                if flavor == "python":
                    nparray = np.rec.array(description)
                else:
                    nparray = array_as_internal(description, flavor)
                self.nrows = nrows = SizeType(nparray.size)
                # If `self._v_recarray` is set, it will be used as the
                # initial buffer.
                if nrows > 0:
                    self._v_recarray = nparray
                tup = descr_from_dtype(
                    nparray.dtype, ptparams=parentnode._v_file.params
                )
                self.description, self._rabyteorder = tup

        # No description yet?
        if new and self.description is None:
            raise TypeError(
                "the ``description`` argument is not of a supported type: "
                "``IsDescription`` subclass, ``Description`` instance, "
                "dictionary, or structured array"
            )

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
        # We are putting here the index-related issues
        # as well as filling general info for table
        # This is needed because we need first the index objects created

        # First, get back the flavor of input data (if any) for
        # `Leaf._g_post_init_hook()`.
        self._flavor, self._descflavor = self._descflavor, None
        super()._g_post_init_hook()

        # Create a cols accessor.
        self.cols = Cols(self, self.description)

        # Place the `Cols` and `Column` objects into `self.colinstances`.
        colinstances, cols = self.colinstances, self.cols
        for colpathname in self.description._v_pathnames:
            colinstances[colpathname] = cols._g_col(colpathname)

        if self._v_new:
            # Columns are never indexed on creation.
            self.colindexed = {cpn: False for cpn in self.colpathnames}
            return

        # The following code is only for opened tables.

        # Do the indexes group exist?
        indexesgrouppath = _index_pathname_of(self)
        igroup = indexesgrouppath in self._v_file
        oldindexes = False
        for colobj in self.description._f_walk(type="Col"):
            colname = colobj._v_pathname
            # Is this column indexed?
            if igroup:
                indexname = _index_pathname_of_column(self, colname)
                indexed = indexname in self._v_file
                self.colindexed[colname] = indexed
                if indexed:
                    column = self.cols._g_col(colname)
                    indexobj = column.index
                    if isinstance(indexobj, OldIndex):
                        indexed = False  # Not a vaild index
                        oldindexes = True
                        self._listoldindexes.append(colname)
                    else:
                        # Tell the condition cache about columns with dirty
                        # indexes.
                        if indexobj.dirty:
                            self._condition_cache.nail()
            else:
                indexed = False
                self.colindexed[colname] = False
            if indexed:
                self.indexed = True

        if oldindexes:  # this should only appear under 2.x Pro
            warnings.warn(
                "table ``%s`` has column indexes with PyTables 1.x format. "
                "Unfortunately, this format is not supported in "
                "PyTables 2.x series. Note that you can use the "
                "``ptrepack`` utility in order to recreate the indexes. "
                "The 1.x indexed columns found are: %s"
                % (self._v_pathname, self._listoldindexes),
                OldIndexWarning,
            )

        # It does not matter to which column 'indexobj' belongs,
        # since their respective index objects share
        # the same number of elements.
        if self.indexed:
            self._indexedrows = indexobj.nelements
            self._unsaved_indexedrows = self.nrows - self._indexedrows
            # Put the autoindex value in a cache variable
            self._autoindex = self.autoindex

    def _calc_nrowsinbuf(self) -> int:
        """Calculate the number of rows that fits on a PyTables buffer."""
        params = self._v_file.params
        # Compute the nrowsinbuf
        rowsize = self.rowsize
        buffersize = params["IO_BUFFER_SIZE"]
        if rowsize != 0:
            nrowsinbuf = buffersize // rowsize
            # The number of rows in buffer needs to be an exact multiple of
            # chunkshape[0] for queries using indexed columns.
            # Fixes #319 and probably #409 too.
            nrowsinbuf -= nrowsinbuf % self.chunkshape[0]
        else:
            nrowsinbuf = 1

        # tableextension.pyx performs an assertion
        # to make sure nrowsinbuf is greater than or
        # equal to the chunksize.
        # See gh-206 and gh-238
        if self.chunkshape is not None:
            if nrowsinbuf < self.chunkshape[0]:
                nrowsinbuf = self.chunkshape[0]

        # Safeguard against row sizes being extremely large
        if nrowsinbuf == 0:
            nrowsinbuf = 1
            # If rowsize is too large, issue a Performance warning
            maxrowsize = params["BUFFER_TIMES"] * buffersize
            if rowsize > maxrowsize:
                warnings.warn(
                    f"""\
The Table ``{self._v_pathname}`` is exceeding the maximum recommended rowsize
({maxrowsize} bytes);
be ready to see PyTables asking for *lots* of memory and possibly slow
I/O.  You may want to reduce the rowsize by trimming the value of
dimensions that are orthogonal (and preferably close) to the *main*
dimension of this leave.  Alternatively, in case you have specified a
very small/large chunksize, you may want to increase/decrease it.""",
                    PerformanceWarning,
                )
        return nrowsinbuf

    def _getemptyarray(self, dtype: np.dtype) -> np.ndarray:
        # Acts as a cache for empty arrays
        key = dtype
        if key in self._empty_array_cache:
            return self._empty_array_cache[key]
        else:
            self._empty_array_cache[key] = arr = np.empty(shape=0, dtype=key)
            return arr

    def _get_container(self, shape: int) -> np.ndarray:
        """Get the appropriate buffer for data depending on table nestedness."""
        # This is *much* faster than the numpy.rec.array counterpart
        return np.empty(shape=shape, dtype=self._v_dtype)

    def _get_type_col_names(self, type_: str) -> list[str]:
        """Return a list containing 'type_' column names."""
        return [
            colobj._v_pathname
            for colobj in self.description._f_walk("Col")
            if colobj.type == type_
        ]

    def _get_enum_map(self) -> dict[str, Enum]:
        """Return mapping from enumerated column names to `Enum` instances."""
        enum_map = {}
        for colobj in self.description._f_walk("Col"):
            if colobj.kind == "enum":
                enum_map[colobj._v_pathname] = colobj.enum
        return enum_map

    def _g_create(self) -> int:
        """Create a new table on disk."""
        # Warning against assigning too much columns...
        # F. Alted 2005-06-05
        max_columns = self._v_file.params["MAX_COLUMNS"]
        if len(self.description._v_names) > max_columns:
            warnings.warn(
                "table ``%s`` is exceeding the recommended "
                "maximum number of columns (%d); "
                "be ready to see PyTables asking for *lots* of memory "
                "and possibly slow I/O" % (self._v_pathname, max_columns),
                PerformanceWarning,
            )

        # 1. Create the HDF5 table (some parameters need to be computed).

        # Fix the byteorder of the recarray and update the number of
        # expected rows if necessary
        if self._v_recarray is not None:
            self._v_recarray = self._g_fix_byteorder_data(
                self._v_recarray, self._rabyteorder
            )
            if len(self._v_recarray) > self._v_expectedrows:
                self._v_expectedrows = len(self._v_recarray)
        # Compute a sensible chunkshape
        if self._v_chunkshape is None:
            self._v_chunkshape = self._calc_chunkshape(
                self._v_expectedrows, self.rowsize, self.rowsize
            )
        # Correct the byteorder, if still needed
        if self.byteorder is None:
            self.byteorder = sys.byteorder

        # Cache some data which is already in the description.
        # This is necessary to happen before creation time in order
        # to be able to populate the self._v_wdflts
        self._cache_description_data()

        # After creating the table, ``self._v_objectid`` needs to be
        # set because it is needed for setting attributes afterwards.
        self._v_objectid = self._create_table(
            self._v_new_title, self.filters.complib or "", obversion
        )
        self._v_recarray = None  # not useful anymore
        self._rabyteorder = None  # not useful anymore

        # 2. Compute or get chunk shape and buffer size parameters.
        self.nrowsinbuf = self._calc_nrowsinbuf()

        # 3. Get field fill attributes from the table description and
        #    set them on disk.
        if self._v_file.params["PYTABLES_SYS_ATTRS"]:
            set_attr = self._v_attrs._g__setattr
            for i, colobj in enumerate(self.description._f_walk(type="Col")):
                fieldname = "FIELD_%d_FILL" % i
                set_attr(fieldname, colobj.dflt)

        return self._v_objectid

    def _g_open(self) -> int:
        """Open a table from disk and read the metadata on it.

        Creates an user description on the flight to easy the access to
        the actual data.

        """
        # 1. Open the HDF5 table and get some data from it.
        self._v_objectid, description, chunksize = self._get_info()
        self._v_expectedrows = self.nrows  # the actual number of rows

        # 2. Create an instance description to host the record fields.
        validate = not self._v_file._isPTFile  # only for non-PyTables files
        self.description = Description(
            description, validate=validate, ptparams=self._v_file.params
        )

        # 3. Compute or get chunk shape and buffer size parameters.
        if chunksize == 0:
            self._v_chunkshape = self._calc_chunkshape(
                self._v_expectedrows, self.rowsize, self.rowsize
            )
        else:
            self._v_chunkshape = (chunksize,)
        self.nrowsinbuf = self._calc_nrowsinbuf()

        # 4. If there are field fill attributes, get them from disk and
        #    set them in the table description.
        if self._v_file.params["PYTABLES_SYS_ATTRS"]:
            if "FIELD_0_FILL" in self._v_attrs._f_list("sys"):
                i = 0
                get_attr = self._v_attrs.__getattr__
                for objcol in self.description._f_walk(type="Col"):
                    colname = objcol._v_pathname
                    # Get the default values for each column
                    fieldname = "FIELD_%s_FILL" % i
                    defval = get_attr(fieldname)
                    if defval is not None:
                        objcol.dflt = defval
                    else:
                        warnings.warn(
                            "could not load default value "
                            "for the ``%s`` column of table ``%s``; "
                            "using ``%r`` instead"
                            % (colname, self._v_pathname, objcol.dflt)
                        )
                        defval = objcol.dflt
                    i += 1

                # Set also the correct value in the desc._v_dflts dictionary
                for descr in self.description._f_walk(type="Description"):
                    for name in descr._v_names:
                        objcol = descr._v_colobjects[name]
                        if isinstance(objcol, Col):
                            descr._v_dflts[objcol._v_name] = objcol.dflt

        # 5. Cache some data which is already in the description.
        self._cache_description_data()

        return self._v_objectid

    def _cache_description_data(self) -> None:
        """Cache some data which is already in the description.

        Some information is extracted from `self.description` to build
        some useful (but redundant) structures:

        * `self.colnames`
        * `self.colpathnames`
        * `self.coldescrs`
        * `self.coltypes`
        * `self.coldtypes`
        * `self.coldflts`
        * `self._v_dtype`
        * `self._time64colnames`
        * `self._strcolnames`
        * `self._colenums`

        """
        self.colnames = list(self.description._v_names)
        self.colpathnames = [
            col._v_pathname
            for col in self.description._f_walk()
            if not hasattr(col, "_v_names")
        ]  # bottom-level

        # Find ``time64`` column names.
        self._time64colnames = self._get_type_col_names("time64")
        # Find ``string`` column names.
        self._strcolnames = self._get_type_col_names("string")
        # Get a mapping of enumerated columns to their `Enum` instances.
        self._colenums = self._get_enum_map()

        # Get info about columns
        for colobj in self.description._f_walk(type="Col"):
            colname = colobj._v_pathname
            # Get the column types, types and defaults
            self.coldescrs[colname] = colobj
            self.coltypes[colname] = colobj.type
            self.coldtypes[colname] = colobj.dtype
            self.coldflts[colname] = colobj.dflt

        # Assign _v_dtype for this table
        self._v_dtype = self.description._v_dtype

    def _get_column_instance(self, colpathname: str):
        """Get the instance of the column with the given `colpathname`.

        If the column does not exist in the table, a `KeyError` is
        raised.

        """
        try:
            return functools.reduce(
                getattr, colpathname.split("/"), self.description
            )
        except AttributeError:
            raise KeyError(
                "table ``%s`` does not have a column named ``%s``"
                % (self._v_pathname, colpathname)
            )

    _check_column = _get_column_instance

    def _disable_indexing_in_queries(self) -> None:
        """Force queries not to use indexing.

        *Use only for testing.*

        """
        if not self._enabled_indexing_in_queries:
            return  # already disabled
        # The nail avoids setting/getting compiled conditions in/from
        # the cache where indexing is used.
        self._condition_cache.nail()
        self._enabled_indexing_in_queries = False

    def _enable_indexing_in_queries(self) -> None:
        """Allow queries to use indexing.

        *Use only for testing.*

        """
        if self._enabled_indexing_in_queries:
            return  # already enabled
        self._condition_cache.unnail()
        self._enabled_indexing_in_queries = True

    def _required_expr_vars(
        self,
        expression: str,
        uservars: dict[str, Column | np.ndarray] | None,
        depth: int = 1,
    ) -> dict[str, Column | np.ndarray]:
        """Get the variables required by the `expression`.

        A new dictionary defining the variables used in the `expression`
        is returned.  Required variables are first looked up in the
        `uservars` mapping, then in the set of top-level columns of the
        table.  Unknown variables cause a `NameError` to be raised.

        When `uservars` is `None`, the local and global namespace where
        the API callable which uses this method is called is sought
        instead.  This mechanism will not work as expected if this
        method is not used *directly* from an API callable.  To disable
        this mechanism, just specify a mapping as `uservars`.

        Nested columns and columns from other tables are not allowed
        (`TypeError` and `ValueError` are raised, respectively).  Also,
        non-column variable values are converted to NumPy arrays.

        `depth` specifies the depth of the frame in order to reach local
        or global variables.

        """
        # Get the names of variables used in the expression.
        exprvarscache = self._exprvars_cache
        if expression not in exprvarscache:
            # Protection against growing the cache too much
            if len(exprvarscache) > 256:
                # Remove 10 (arbitrary) elements from the cache
                for k in list(exprvarscache)[:10]:
                    del exprvarscache[k]
            cexpr = compile(expression, "<string>", "eval")
            exprvars = [
                var
                for var in cexpr.co_names
                if var not in ["None", "False", "True"]
                and var not in ne.expressions.functions
            ]
            exprvarscache[expression] = exprvars
        else:
            exprvars = exprvarscache[expression]

        # Get the local and global variable mappings of the user frame
        # if no mapping has been explicitly given for user variables.
        user_locals, user_globals = {}, {}
        if uservars is None:
            # We use specified depth to get the frame where the API
            # callable using this method is called.  For instance:
            #
            # * ``table._required_expr_vars()`` (depth 0) is called by
            # * ``table._where()`` (depth 1) is called by
            # * ``table.where()`` (depth 2) is called by
            # * user-space functions (depth 3)
            user_frame = sys._getframe(depth)
            user_locals = user_frame.f_locals
            user_globals = user_frame.f_globals

        colinstances = self.colinstances
        tblfile, tblpath = self._v_file, self._v_pathname
        # Look for the required variables first among the ones
        # explicitly provided by the user, then among implicit columns,
        # then among external variables (only if no explicit variables).
        reqvars: dict[str, np.ndarray] = {}
        for var in exprvars:
            # Get the value.
            if uservars is not None and var in uservars:
                val = uservars[var]
            elif var in colinstances:
                val = colinstances[var]
            elif uservars is None and var in user_locals:
                val = user_locals[var]
            elif uservars is None and var in user_globals:
                val = user_globals[var]
            else:
                raise NameError("name ``%s`` is not defined" % var)

            # Check the value.
            if hasattr(val, "pathname"):  # non-nested column
                if val.shape[1:] != ():
                    raise NotImplementedError(
                        f"variable ``{var}`` refers to a multidimensional "
                        f"column, not yet supported in conditions, sorry"
                    )
                if (
                    val._table_file is not tblfile
                    or val._table_path != tblpath
                ):
                    raise ValueError(
                        f"variable ``{var}`` refers to a column "
                        f"which is not part of table ``{tblpath}``"
                    )
                if val.dtype.str[1:] == "u8":
                    raise NotImplementedError(
                        f"variable ``{var}`` refers to a 64-bit unsigned "
                        f"integer column, not yet supported in conditions, "
                        f"sorry; please use regular Python selections"
                    )
            elif hasattr(val, "_v_colpathnames"):  # nested column
                raise TypeError(
                    f"variable ``{var}`` refers to a nested column, "
                    f"not allowed in conditions"
                )
            else:  # only non-column values are converted to arrays
                # XXX: not 100% sure about this
                if isinstance(val, str):
                    val = np.asarray(val.encode("ascii"))
                else:
                    val = np.asarray(val)
            reqvars[var] = val
        return reqvars

    def _get_condition_key(
        self,
        condition: str,
        condvars: dict[str, Column],
    ) -> tuple[
        str,
        tuple[str, ...],
        tuple[str, ...],
        tuple[str, ...],
        tuple[Any, ...],
    ]:
        """Get the condition cache key for `condition` with `condvars`.

        Currently, the key is a tuple of `condition`, column variables
        names, normal variables names, column paths and variable paths
        (all are tuples).

        """
        # Variable names for column and normal variables.
        colnames, varnames = [], []
        # Column paths and types for each of the previous variable.
        colpaths, vartypes = [], []
        for var, val in condvars.items():
            if hasattr(val, "pathname"):  # column
                colnames.append(var)
                colpaths.append(val.pathname)
            else:  # array
                try:
                    varnames.append(var)
                    vartypes.append(ne.necompiler.getType(val))  # expensive
                except ValueError:
                    # This is more clear than the error given by Numexpr.
                    raise TypeError(
                        "variable ``%s`` has data type ``%s``, "
                        "not allowed in conditions" % (var, val.dtype.name)
                    )
        colnames, varnames = tuple(colnames), tuple(varnames)
        colpaths, vartypes = tuple(colpaths), tuple(vartypes)
        condkey = (condition, colnames, varnames, colpaths, vartypes)
        return condkey

    def _compile_condition(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray],
    ) -> CompiledCondition:
        """Compile the `condition` and extract usable index conditions.

        This method returns an instance of ``CompiledCondition``.  See
        the ``compile_condition()`` function in the ``conditions``
        module for more information about the compilation process.

        This method makes use of the condition cache when possible.

        """
        # Look up the condition in the condition cache.
        condcache = self._condition_cache
        condkey = self._get_condition_key(condition, condvars)
        compiled = condcache.get(condkey)
        if compiled:
            return compiled.with_replaced_vars(condvars)  # bingo!

        # Bad luck, the condition must be parsed and compiled.
        # Fortunately, the key provides some valuable information. ;)
        (condition, colnames, varnames, colpaths, vartypes) = condkey

        # Extract more information from referenced columns.

        # start with normal variables
        typemap = dict(list(zip(varnames, vartypes)))
        indexedcols = []
        for colname in colnames:
            col = condvars[colname]

            # Extract types from *all* the given variables.
            coltype = col.dtype.type
            typemap[colname] = _nxtype_from_nptype[coltype]

            # Get the set of columns with usable indexes.
            if (
                self._enabled_indexing_in_queries  # no in-kernel searches
                and self.colindexed[col.pathname]
                and not col.index.dirty
            ):
                indexedcols.append(colname)

        indexedcols = frozenset(indexedcols)
        # Now let ``compile_condition()`` do the Numexpr-related job.
        compiled = compile_condition(condition, typemap, indexedcols)

        # Check that there actually are columns in the condition.
        if not set(compiled.parameters).intersection(set(colnames)):
            raise ValueError(
                f"there are no columns taking part in "
                f"condition ``{condition}``"
            )

        # Store the compiled condition in the cache and return it.
        condcache[condkey] = compiled
        return compiled.with_replaced_vars(condvars)

    def will_query_use_indexing(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray] | None = None,
    ) -> frozenset:
        """Return True if the query for the condition will use indexing.

        The meaning of the condition and *condvars* arguments is the same as in
        the :meth:`Table.where` method. If condition can use indexing, this
        method returns a frozenset with the path names of the columns whose
        index is usable. Otherwise, it returns an empty list.

        This method is mainly intended for testing. Keep in mind that changing
        the set of indexed columns or their dirtiness may make this method
        return different values for the same arguments at different times.

        """
        # Compile the condition and extract usable index conditions.
        condvars = self._required_expr_vars(condition, condvars, depth=2)
        compiled = self._compile_condition(condition, condvars)
        # Return the columns in indexed expressions
        idxcols = [condvars[var].pathname for var in compiled.index_variables]
        return frozenset(idxcols)

    def where(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray] | None = None,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> Iterator[tableextension.Row]:
        r"""Iterate over values fulfilling a condition.

        This method returns a Row iterator (see :ref:`RowClassDescr`) which
        only selects rows in the table that satisfy the given condition (an
        expression-like string).

        The condvars mapping may be used to define the variable names appearing
        in the condition. condvars should consist of identifier-like strings
        pointing to Column (see :ref:`ColumnClassDescr`) instances *of this
        table*, or to other values (which will be converted to arrays). A
        default set of condition variables is provided where each top-level,
        non-nested column with an identifier-like name appears. Variables in
        condvars override the default ones.

        When condvars is not provided or None, the current local and global
        namespace is sought instead of condvars. The previous mechanism is
        mostly intended for interactive usage. To disable it, just specify a
        (maybe empty) mapping as condvars.

        If a range is supplied (by setting some of the start, stop or step
        parameters), only the rows in that range and fulfilling the condition
        are used. The meaning of the start, stop and step parameters is the
        same as for Python slices.

        When possible, indexed columns participating in the condition will be
        used to speed up the search. It is recommended that you place the
        indexed columns as left and out in the condition as possible. Anyway,
        this method has always better performance than regular Python
        selections on the table.

        You can mix this method with regular Python selections in order to
        support even more complex queries. It is strongly recommended that you
        pass the most restrictive condition as the parameter to this method if
        you want to achieve maximum performance.

        .. warning::

            When in the middle of a table row iterator, you should not
            use methods that can change the number of rows in the table
            (like :meth:`Table.append` or :meth:`Table.remove_rows`) or
            unexpected errors will happen.

        Examples
        --------
        ::

            passvalues = [ row['col3'] for row in
                           table.where('(col1 > 0) & (col2 <= 20)', step=5)
                           if your_function(row['col2']) ]
            print("Values that pass the cuts:", passvalues)

        .. note::

            A special care should be taken when the query condition includes
            string literals.

            Let's assume that the table ``table`` has the following
            structure::

                class Record(IsDescription):
                    col1 = StringCol(4)  # 4-character String of bytes
                    col2 = IntCol()
                    col3 = FloatCol()

            The type of "col1" corresponds to strings of bytes.

            Any condition involving "col1" should be written using the
            appropriate type for string literals in order to avoid
            :exc:`TypeError`\ s.

            The code below will fail with a :exc:`TypeError`::

                condition = 'col1 == "AAAA"'
                for record in table.where(condition):  # TypeError in Python3
                    # do something with "record"

            The reason is that in Python 3 "condition" implies a comparison
            between a string of bytes ("col1" contents) and a unicode literal
            ("AAAA").

            The correct way to write the condition is::

                condition = 'col1 == b"AAAA"'

        .. versionchanged:: 3.0
           The start, stop and step parameters now behave like in slice.

        """
        return self._where(condition, condvars, start, stop, step)

    def _where(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray] | None,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> Iterator[tableextension.Row]:
        """Low-level counterpart of `self.where()`."""
        if profile:
            tref = clock()
        if profile:
            show_stats("Entering table._where", tref)
        # Adjust the slice to be used.
        (start, stop, step) = self._process_range_read(start, stop, step)
        if start >= stop:  # empty range, reset conditions
            self._use_index = False
            self._where_condition = None
            return iter([])

        # Compile the condition and extract usable index conditions.
        condvars = self._required_expr_vars(condition, condvars, depth=3)
        compiled = self._compile_condition(condition, condvars)

        # Can we use indexes?
        if compiled.index_expressions:
            chunkmap = _table__where_indexed(
                self, compiled, condition, condvars, start, stop, step
            )
            if not isinstance(chunkmap, np.ndarray):
                # If it is not a NumPy array it should be an iterator
                # Reset conditions
                self._use_index = False
                self._where_condition = None
                # ...and return the iterator
                return chunkmap
        else:
            chunkmap = None  # default to an in-kernel query

        args = [condvars[param] for param in compiled.parameters]
        self._where_condition = (compiled.function, args, compiled.kwargs)
        row = tableextension.Row(self)
        if profile:
            show_stats("Exiting table._where", tref)
        return row._iter(start, stop, step, chunkmap=chunkmap)

    def read_where(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray] | None = None,
        field=None,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> np.ndarray:
        """Read table data fulfilling the given *condition*.

        This method is similar to :meth:`Table.read`, having their common
        arguments and return values the same meanings. However, only the rows
        fulfilling the *condition* are included in the result.

        The meaning of the other arguments is the same as in the
        :meth:`Table.where` method.

        """
        self._g_check_open()
        coords = [
            p.nrow for p in self._where(condition, condvars, start, stop, step)
        ]
        self._where_condition = None  # reset the conditions
        if len(coords) > 1:
            cstart, cstop = coords[0], coords[-1] + 1
            if cstop - cstart == len(coords):
                # Chances for monotonically increasing row values. Refine.
                inc_seq = np.all(np.arange(cstart, cstop) == np.array(coords))
                if inc_seq:
                    return self.read(cstart, cstop, field=field)
        return self.read_coordinates(coords, field)

    def append_where(
        self,
        dstTable: Table,  # noqa: N803
        condition: str | None = None,
        condvars: dict[str, Column | np.ndarray] | None = None,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> int:
        """Append rows fulfilling the condition to the dstTable table.

        dstTable must be capable of taking the rows resulting from the query,
        i.e. it must have columns with the expected names and compatible
        types. The meaning of the other arguments is the same as in the
        :meth:`Table.where` method.

        The number of rows appended to dstTable is returned as a result.

        .. versionchanged:: 3.0
           The *whereAppend* method has been renamed into *append_where*.

        """
        self._g_check_open()

        # Check that the destination file is not in read-only mode.
        dstTable._v_file._check_writable()

        # Row objects do not support nested columns, so we must iterate
        # over the flat column paths.  When rows support nesting,
        # ``self.colnames`` can be directly iterated upon.
        col_names = tuple(col_name for col_name in self.colpathnames)
        dst_row = dstTable.row
        nrows = 0
        if condition is not None:
            src_rows = self._where(condition, condvars, start, stop, step)
        else:
            src_rows = self.iterrows(start, stop, step)
        for src_row in src_rows:
            for col_name in col_names:
                dst_row[col_name] = src_row[col_name]
            dst_row.append()
            nrows += 1
        dstTable.flush()
        return nrows

    def get_where_list(
        self,
        condition: str,
        condvars: dict[str, Column | np.ndarray] | None = None,
        sort: bool = False,
        start: str | None = None,
        stop: str | None = None,
        step: str | None = None,
    ) -> np.ndarray:
        """Get the row coordinates fulfilling the given condition.

        The coordinates are returned as a list of the current flavor.  sort
        means that you want to retrieve the coordinates ordered. The default is
        to not sort them.

        The meaning of the other arguments is the same as in the
        :meth:`Table.where` method.

        """
        self._g_check_open()

        coords = [
            p.nrow for p in self._where(condition, condvars, start, stop, step)
        ]
        coords = np.array(coords, dtype=SizeType)
        # Reset the conditions
        self._where_condition = None
        if sort:
            coords = np.sort(coords)
        return internal_to_flavor(coords, self.flavor)

    def itersequence(self, sequence: Sequence) -> Iterator[tableextension.Row]:
        """Iterate over a sequence of row coordinates."""
        if not hasattr(sequence, "__getitem__"):
            raise TypeError(
                "Wrong 'sequence' parameter type. Only sequences "
                "are suported."
            )
        # start, stop and step are necessary for the new iterator for
        # coordinates, and perhaps it would be useful to add them as
        # parameters in the future (not now, because I've just removed
        # the `sort` argument for 2.1).
        #
        # *Important note*: Negative values for step are not supported
        # for the general case, but only for the itersorted() and
        # read_sorted() purposes!  The self._process_range_read will raise
        # an appropriate error.
        # F. Alted 2008-09-18
        # A.V. 20130513: _process_range_read --> _process_range
        (start, stop, step) = self._process_range(None, None, None)
        if (start > stop) or (len(sequence) == 0):
            return iter([])
        row = tableextension.Row(self)
        return row._iter(start, stop, step, coords=sequence)

    def _check_sortby_csi(
        self, sortby: Column | str, check_csi: bool
    ) -> Index:
        if isinstance(sortby, Column):
            icol = sortby
        elif isinstance(sortby, str):
            icol = self.cols._f_col(sortby)
        else:
            raise TypeError(
                f"`sortby` can only be a `Column` or string object, "
                f"but you passed an object of type: {type(sortby)}"
            )
        if icol.is_indexed and icol.index.kind == "full":
            if check_csi and not icol.index.is_csi:
                # The index exists, but it is not a CSI one.
                raise ValueError(
                    "Field `{sortby}` must have associated a CSI index "
                    "in table `{self}`, but the existing one is not. "
                )
            return icol.index
        else:
            raise ValueError(
                f"Field `{sortby}` must have associated a 'full' index "
                f"in table `{self}`."
            )

    def itersorted(
        self,
        sortby: Column | str,
        checkCSI: bool = False,  # noqa: N803
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> Iterator[tableextension.Row]:
        """Iterate table data in the order of the index of sortby column.

        The sortby column must have associated a full index.  If you want to
        ensure a fully sorted order, the index must be a CSI one.  You may want
        to use the checkCSI argument in order to explicitly check for the
        existence of a CSI index.

        The meaning of the start, stop and step arguments is the same as in
        :meth:`Table.read`.

        .. versionchanged:: 3.0
           If the *start* parameter is provided and *stop* is None then the
           table is iterated from *start* to the last line.
           In PyTables < 3.0 only one element was returned.

        """
        index = self._check_sortby_csi(sortby, checkCSI)
        # Adjust the slice to be used.
        (start, stop, step) = self._process_range(
            start, stop, step, warn_negstep=False
        )
        if (start > stop and 0 < step) or (start < stop and 0 > step):
            # Fall-back action is to return an empty iterator
            return iter([])
        row = tableextension.Row(self)
        return row._iter(start, stop, step, coords=index)

    def read_sorted(
        self,
        sortby: Column | str,
        checkCSI: bool = False,  # noqa: N803
        field: str | None = None,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> np.ndarray:
        """Read table data following the order of the index of sortby column.

        The sortby column must have associated a full index.  If you want to
        ensure a fully sorted order, the index must be a CSI one.  You may want
        to use the checkCSI argument in order to explicitly check for the
        existence of a CSI index.

        If field is supplied only the named column will be selected.  If the
        column is not nested, an *array* of the current flavor will be
        returned; if it is, a *structured array* will be used instead.  If no
        field is specified, all the columns will be returned in a structured
        array of the current flavor.

        The meaning of the start, stop and step arguments is the same as in
        :meth:`Table.read`.

        .. versionchanged:: 3.0
           The start, stop and step parameters now behave like in slice.

        """
        self._g_check_open()
        index = self._check_sortby_csi(sortby, checkCSI)
        coords = index[start:stop:step]
        return self.read_coordinates(coords, field)

    def iterrows(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> Iterator[tableextension.Row]:
        """Iterate over the table using a Row instance.

        If a range is not supplied, *all the rows* in the table are iterated
        upon - you can also use the :meth:`Table.__iter__` special method for
        that purpose. If you want to iterate over a given *range of rows* in
        the table, you may use the start, stop and step parameters.

        .. warning::

            When in the middle of a table row iterator, you should not
            use methods that can change the number of rows in the table
            (like :meth:`Table.append` or :meth:`Table.remove_rows`) or
            unexpected errors will happen.

        See Also
        --------
        tableextension.Row : the table row iterator and field accessor

        Examples
        --------
        ::

            result = [ row['var2'] for row in table.iterrows(step=5)
                                                    if row['var1'] <= 20 ]

        .. versionchanged:: 3.0
           If the *start* parameter is provided and *stop* is None then the
           table is iterated from *start* to the last line.
           In PyTables < 3.0 only one element was returned.

        """
        (start, stop, step) = self._process_range(
            start, stop, step, warn_negstep=False
        )
        if (start > stop and 0 < step) or (start < stop and 0 > step):
            # Fall-back action is to return an empty iterator
            return iter([])
        row = tableextension.Row(self)
        return row._iter(start, stop, step)

    def __iter__(self) -> Iterator[tableextension.Row]:
        """Iterate over the table using a Row instance.

        This is equivalent to calling :meth:`Table.iterrows` with default
        arguments, i.e. it iterates over *all the rows* in the table.

        See Also
        --------
        tableextension.Row : the table row iterator and field accessor

        Examples
        --------
        ::

            result = [ row['var2'] for row in table if row['var1'] <= 20 ]

        Which is equivalent to::

            result = [ row['var2'] for row in table.iterrows()
                                                    if row['var1'] <= 20 ]

        """
        return self.iterrows()

    def _read(
        self,
        start: int,
        stop: int,
        step: int,
        field: str | None = None,
        out: np.ndarray | None = None,
    ) -> np.ndarray:
        """Read a range of rows and return an in-memory object."""
        select_field = None
        if field:
            if field not in self.coldtypes:
                if field in self.description._v_names:
                    # Remember to select this field
                    select_field = field
                    field = None
                else:
                    raise KeyError(
                        ("Field {} not found in table " "{}").format(
                            field, self
                        )
                    )
            else:
                # The column hangs directly from the top
                dtype_field = self.coldtypes[field]

        # Return a rank-0 array if start > stop
        if (start >= stop and 0 < step) or (start <= stop and 0 > step):
            if field is None:
                nra = self._get_container(0)
                return nra
            return np.empty(shape=0, dtype=dtype_field)

        nrows = len(range(start, stop, step))

        if out is None:
            # Compute the shape of the resulting column object
            if field:
                # Create a container for the results
                result = np.empty(shape=nrows, dtype=dtype_field)
            else:
                # Recarray case
                result = self._get_container(nrows)
        else:
            # there is no fast way to byteswap, since different columns may
            # have different byteorders
            if not out.dtype.isnative:
                raise ValueError(
                    "output array must be in system's byteorder "
                    "or results will be incorrect"
                )
            if field:
                bytes_required = dtype_field.itemsize * nrows
            else:
                bytes_required = self.rowsize * nrows
            if bytes_required != out.nbytes:
                raise ValueError(
                    f"output array size invalid, got {out.nbytes}"
                    f" bytes, need {bytes_required} bytes"
                )
            if not out.flags["C_CONTIGUOUS"]:
                raise ValueError("output array not C contiguous")
            result = out

        # Call the routine to fill-up the resulting array
        if step == 1 and not field:
            # This optimization works three times faster than
            # the row._fill_col method (up to 170 MB/s on a pentium IV @ 2GHz)
            self._read_records(start, stop - start, result)
        # Warning!: _read_field_name should not be used until
        # H5TBread_fields_name in tableextension will be finished
        # F. Alted 2005/05/26
        # XYX Ho implementem per a PyTables 2.0??
        elif field and step > 15 and 0:
            # For step>15, this seems to work always faster than row._fill_col.
            self._read_field_name(result, start, stop, step, field)
        else:
            self.row._fill_col(result, start, stop, step, field)

        if select_field:
            return result[select_field]
        else:
            return result

    def read(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        field: str | None = None,
        out: np.ndarray | None = None,
    ) -> np.ndarray:
        """Get data in the table as a (record) array.

        The start, stop and step parameters can be used to select only
        a *range of rows* in the table. Their meanings are the same as
        in the built-in Python slices.

        If field is supplied only the named column will be selected.
        If the column is not nested, an *array* of the current flavor
        will be returned; if it is, a *structured array* will be used
        instead.  If no field is specified, all the columns will be
        returned in a structured array of the current flavor.

        Columns under a nested column can be specified in the field
        parameter by using a slash character (/) as a separator (e.g.
        'position/x').

        The out parameter may be used to specify a NumPy array to
        receive the output data.  Note that the array must have the
        same size as the data selected with the other parameters.
        Note that the array's datatype is not checked and no type
        casting is performed, so if it does not match the datatype on
        disk, the output will not be correct.

        When specifying a single nested column with the field parameter,
        and supplying an output buffer with the out parameter, the
        output buffer must contain all columns in the table.
        The data in all columns will be read into the output buffer.
        However, only the specified nested column will be returned from
        the method call.

        When data is read from disk in NumPy format, the output will be
        in the current system's byteorder, regardless of how it is
        stored on disk. If the out parameter is specified, the output
        array also must be in the current system's byteorder.

        .. versionchanged:: 3.0
           Added the *out* parameter.  Also the start, stop and step
           parameters now behave like in slice.

        Examples
        --------
        Reading the entire table::

            t.read()

        Reading record n. 6::

            t.read(6, 7)

        Reading from record n. 6 to the end of the table::

            t.read(6)

        """
        self._g_check_open()

        if field:
            self._check_column(field)

        if out is not None and self.flavor != "numpy":
            msg = (
                f"Optional 'out' argument may only be supplied if array "
                f"flavor is 'numpy', currently is {self.flavor}"
            )
            raise TypeError(msg)

        start, stop, step = self._process_range(
            start, stop, step, warn_negstep=False
        )

        arr = self._read(start, stop, step, field, out)
        return internal_to_flavor(arr, self.flavor)

    def _read_coordinates(
        self, coords: np.ndarray, field: str | None = None
    ) -> np.ndarray:
        """Private part of `read_coordinates()` with no flavor conversion."""
        coords = self._point_selection(coords)

        ncoords = len(coords)
        # Create a read buffer only if needed
        if field is None or ncoords > 0:
            # Doing a copy is faster when ncoords is small (<1000)
            if ncoords < min(1000, self.nrowsinbuf):
                result = self._v_iobuf[:ncoords].copy()
            else:
                result = self._get_container(ncoords)

        # Do the real read
        if ncoords > 0:
            # Turn coords into an array of coordinate indexes, if necessary
            if not (
                isinstance(coords, np.ndarray)
                and coords.dtype.type is _npsizetype
                and coords.flags.contiguous
                and coords.flags.aligned
            ):
                # Get a contiguous and aligned coordinate array
                coords = np.array(coords, dtype=SizeType)
            self._read_elements(coords, result)

        # Do the final conversions, if needed
        if field:
            if ncoords > 0:
                result = get_nested_field(result, field)
            else:
                # Get an empty array from the cache
                result = self._getemptyarray(self.coldtypes[field])
        return result

    def read_coordinates(
        self, coords: np.ndarray, field: str | None = None
    ) -> np.ndarray:
        """Get a set of rows given their indexes as a (record) array.

        This method works much like the :meth:`Table.read` method, but it uses
        a sequence (coords) of row indexes to select the wanted columns,
        instead of a column range.

        The selected rows are returned in an array or structured array of the
        current flavor.

        """
        self._g_check_open()
        result = self._read_coordinates(coords, field)
        return internal_to_flavor(result, self.flavor)

    def get_enum(self, colname: str) -> Enum:
        """Get the enumerated type associated with the named column.

        If the column named colname (a string) exists and is of an enumerated
        type, the corresponding Enum instance (see :ref:`EnumClassDescr`) is
        returned. If it is not of an enumerated type, a TypeError is raised. If
        the column does not exist, a KeyError is raised.

        """
        self._check_column(colname)

        try:
            return self._colenums[colname]
        except KeyError:
            raise TypeError(
                "column ``%s`` of table ``%s`` is not of an enumerated type"
                % (colname, self._v_pathname)
            )

    def col(self, name: str) -> np.ndarray:
        """Get a column from the table.

        If a column called name exists in the table, it is read and returned as
        a NumPy object. If it does not exist, a KeyError is raised.

        Examples
        --------
        ::

            narray = table.col('var2')

        That statement is equivalent to::

            narray = table.read(field='var2')

        Here you can see how this method can be used as a shorthand for the
        :meth:`Table.read` method.

        """
        return self.read(field=name)

    def __getitem__(
        self, key: int | slice | list[int] | list[bool] | np.ndarray
    ) -> np.ndarray:
        """Get a row or a range of rows from the table.

        If key argument is an integer, the corresponding table row is returned
        as a record of the current flavor. If key is a slice, the range of rows
        determined by it is returned as a structured array of the current
        flavor.

        In addition, NumPy-style point selections are supported.  In
        particular, if key is a list of row coordinates, the set of rows
        determined by it is returned.  Furthermore, if key is an array of
        boolean values, only the coordinates where key is True are returned.
        Note that for the latter to work it is necessary that key list would
        contain exactly as many rows as the table has.

        Examples
        --------
        ::

            record = table[4]
            recarray = table[4:1000:2]
            recarray = table[[4,1000]]   # only retrieves rows 4 and 1000
            recarray = table[[True, False, ..., True]]

        Those statements are equivalent to::

            record = table.read(start=4)[0]
            recarray = table.read(start=4, stop=1000, step=2)
            recarray = table.read_coordinates([4,1000])
            recarray = table.read_coordinates([True, False, ..., True])

        Here, you can see how indexing can be used as a shorthand for the
        :meth:`Table.read` and :meth:`Table.read_coordinates` methods.

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
            (start, stop, step) = self._process_range(
                key.start, key.stop, key.step
            )
            return self.read(start, stop, step)
        # Try with a boolean or point selection
        elif type(key) in (list, tuple) or isinstance(key, np.ndarray):
            return self._read_coordinates(key, None)
        else:
            raise IndexError(f"Invalid index or slice: {key!r}")

    def __setitem__(
        self,
        key: int | slice | list[int] | list[bool] | np.ndarray,
        value: Any,
    ) -> int:
        """Set a row or a range of rows in the table.

        It takes different actions depending on the type of the *key*
        parameter: if it is an integer, the corresponding table row is
        set to *value* (a record or sequence capable of being converted
        to the table structure).  If *key* is a slice, the row slice
        determined by it is set to *value* (a record array or sequence
        capable of being converted to the table structure).

        In addition, NumPy-style point selections are supported.  In
        particular, if key is a list of row coordinates, the set of rows
        determined by it is set to value.  Furthermore, if key is an array of
        boolean values, only the coordinates where key is True are set to
        values from value.  Note that for the latter to work it is necessary
        that key list would contain exactly as many rows as the table has.

        Examples
        --------
        ::

            # Modify just one existing row
            table[2] = [456,'db2',1.2]

            # Modify two existing rows
            rows = np.rec.array(
                [[457,'db1',1.2],[6,'de2',1.3]], formats='i4,S3,f8'
            )
            table[1:30:2] = rows             # modify a table slice
            table[[1,3]] = rows              # only modifies rows 1 and 3
            table[[True,False,True]] = rows  # only modifies rows 0 and 2

        Which is equivalent to::

            table.modify_rows(start=2, rows=[456,'db2',1.2])
            rows = np.rec.array(
                [[457,'db1',1.2],[6,'de2',1.3]], formats='i4,S3,f8'
            )
            table.modify_rows(start=1, stop=3, step=2, rows=rows)
            table.modify_coordinates([1,3,2], rows)
            table.modify_coordinates([True, False, True], rows)

        Here, you can see how indexing can be used as a shorthand for the
        :meth:`Table.modify_rows` and :meth:`Table.modify_coordinates`
        methods.

        """
        self._g_check_open()
        self._v_file._check_writable()

        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= self.nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += self.nrows
            return self.modify_rows(key, key + 1, 1, [value])
        elif isinstance(key, slice):
            (start, stop, step) = self._process_range(
                key.start, key.stop, key.step
            )
            return self.modify_rows(start, stop, step, value)
        # Try with a boolean or point selection
        elif type(key) in (list, tuple) or isinstance(key, np.ndarray):
            return self.modify_coordinates(key, value)
        else:
            raise IndexError(f"Invalid index or slice: {key!r}")

    def _save_buffered_rows(self, wbuf_ra: np.ndarray, lenrows: int) -> None:
        """Update the indexes after a flushing of rows."""
        self._open_append(wbuf_ra)
        self._append_records(lenrows)
        self._close_append()
        if self.indexed:
            self._unsaved_indexedrows += lenrows
            # The table caches for indexed queries are dirty now
            self._dirtycache = True
            if self.autoindex:
                # Flush the unindexed rows
                self.flush_rows_to_index(_lastrow=False)
            else:
                # All the columns are dirty now
                self._mark_columns_as_dirty(self.colpathnames)

    def append(self, rows: list | np.ndarray) -> None:
        """Append a sequence of rows to the end of the table.

        The rows argument may be any object which can be converted to
        a structured array compliant with the table structure
        (otherwise, a ValueError is raised).  This includes NumPy
        structured arrays, lists of tuples or array records, and a
        string or Python buffer.

        Examples
        --------
        ::

            import tables as tb

            class Particle(tb.IsDescription):
                name        = tb.StringCol(16, pos=1) # 16-character String
                lati        = tb.IntCol(pos=2)        # integer
                longi       = tb.IntCol(pos=3)        # integer
                pressure    = tb.Float32Col(pos=4)  # float  (single-precision)
                temperature = tb.FloatCol(pos=5)    # double (double-precision)

            fileh = tb.open_file('test4.h5', mode='w')
            table = fileh.create_table(fileh.root, 'table', Particle,
                                       "A table")

            # Append several rows in only one call
            table.append([("Particle:     10", 10, 0, 10 * 10, 10**2),
                          ("Particle:     11", 11, -1, 11 * 11, 11**2),
                          ("Particle:     12", 12, -2, 12 * 12, 12**2)])
            fileh.close()

        """
        self._g_check_open()
        self._v_file._check_writable()

        if not self._chunked:
            raise HDF5ExtError(
                "You cannot append rows to a non-chunked table.", h5bt=False
            )

        if (
            hasattr(rows, "dtype")
            and not self.description._v_is_nested
            and rows.dtype == self.dtype
        ):
            # Shortcut for compliant arrays
            # (for some reason, not valid for nested types)
            wbuf_ra = rows
        else:
            # Try to convert the object into a recarray compliant with table
            try:
                iflavor = flavor_of(rows)
                if iflavor != "python":
                    rows = array_as_internal(rows, iflavor)
                # Works for Python structures and always copies the original,
                # so the resulting object is safe for in-place conversion.
                wbuf_ra = np.rec.array(rows, dtype=self._v_dtype)
            except Exception as exc:  # XXX
                raise ValueError(
                    f"rows parameter cannot be converted into a "
                    f"recarray object compliant with table '{self}'. "
                    f"The error was: <{exc}>"
                )
        lenrows = wbuf_ra.shape[0]
        # If the number of rows to append is zero, don't do anything else
        if lenrows > 0:
            # Save write buffer to disk
            self._save_buffered_rows(wbuf_ra, lenrows)

    def _conv_to_recarr(self, obj: Sequence) -> np.ndarray:
        """Try to convert the object into a recarray."""
        try:
            iflavor = flavor_of(obj)
            if iflavor != "python":
                obj = array_as_internal(obj, iflavor)
            if hasattr(obj, "shape") and obj.shape == ():
                # To allow conversion of scalars (void type) into arrays.
                # See http://projects.scipy.org/scipy/numpy/ticket/315
                # for discussion on how to pass buffers to constructors
                # See also http://projects.scipy.org/scipy/numpy/ticket/348
                recarr = np.array([obj], dtype=self._v_dtype)
            else:
                # Works for Python structures and always copies the original,
                # so the resulting object is safe for in-place conversion.
                recarr = np.rec.array(obj, dtype=self._v_dtype)
        except Exception as exc:  # XXX
            raise ValueError(
                f"Object cannot be converted into a recarray object compliant "
                f"with table format '{self.description._v_nested_descr}'. "
                f"The error was: <{exc}>"
            )

        return recarr

    def modify_coordinates(
        self, coords: list | tuple | np.ndarray, rows: Sequence
    ) -> int:
        """Modify a series of rows in positions specified in coords.

        The values in the selected rows will be modified with the data given in
        rows.  This method returns the number of rows modified.

        The possible values for the rows argument are the same as in
        :meth:`Table.append`.

        """
        if rows is None:  # Nothing to be done
            return SizeType(0)

        # Convert the coordinates to something expected by HDF5
        coords = self._point_selection(coords)

        lcoords = len(coords)
        if len(rows) < lcoords:
            raise ValueError(
                "The value has not enough elements to fill-in "
                "the specified range"
            )

        # Convert rows into a recarray
        recarr = self._conv_to_recarr(rows)

        if len(coords) > 0:
            # Do the actual update of rows
            self._update_elements(lcoords, coords, recarr)

        # Redo the index if needed
        self._reindex(self.colpathnames)

        return SizeType(lcoords)

    def modify_rows(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        rows: Sequence | None = None,
    ) -> int:
        """Modify a series of rows in the slice [start:stop:step].

        The values in the selected rows will be modified with the data given in
        rows.  This method returns the number of rows modified.  Should the
        modification exceed the length of the table, an IndexError is raised
        before changing data.

        The possible values for the rows argument are the same as in
        :meth:`Table.append`.

        """
        if step is None:
            step = 1
        if rows is None:  # Nothing to be done
            return SizeType(0)
        if start is None:
            start = 0

        if start < 0:
            raise ValueError("'start' must have a positive value.")
        if step < 1:
            raise ValueError(
                "'step' must have a value greater or equal than 1."
            )
        if stop is None:
            # compute the stop value. start + len(rows)*step does not work
            stop = start + (len(rows) - 1) * step + 1

        (start, stop, step) = self._process_range(start, stop, step)
        if stop > self.nrows:
            raise IndexError(
                "This modification will exceed the length of "
                "the table. Giving up."
            )
        # Compute the number of rows to read.
        nrows = len(range(start, stop, step))
        if len(rows) != nrows:
            raise ValueError(
                "The value has different elements than the specified range"
            )

        # Convert rows into a recarray
        recarr = self._conv_to_recarr(rows)

        lenrows = len(recarr)
        if start + lenrows > self.nrows:
            raise IndexError(
                "This modification will exceed the length of the "
                "table. Giving up."
            )

        # Do the actual update
        self._update_records(start, stop, step, recarr)

        # Redo the index if needed
        self._reindex(self.colpathnames)

        return SizeType(lenrows)

    def modify_column(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        column: Sequence | None = None,
        colname: str | None = None,
    ):
        """Modify one single column in the row slice [start:stop:step].

        The colname argument specifies the name of the column in the
        table to be modified with the data given in column.  This
        method returns the number of rows modified.  Should the
        modification exceed the length of the table, an IndexError is
        raised before changing data.

        The *column* argument may be any object which can be converted
        to a (record) array compliant with the structure of the column
        to be modified (otherwise, a ValueError is raised).  This
        includes NumPy (record) arrays, lists of scalars, tuples or
        array records, and a string or Python buffer.

        """
        if step is None:
            step = 1
        if not isinstance(colname, str):
            raise TypeError("The 'colname' parameter must be a string.")
        self._v_file._check_writable()

        if column is None:  # Nothing to be done
            return SizeType(0)
        if start is None:
            start = 0

        if start < 0:
            raise ValueError("'start' must have a positive value.")
        if step < 1:
            raise ValueError(
                "'step' must have a value greater or equal than 1."
            )
        # Get the column format to be modified:
        objcol = self._get_column_instance(colname)
        descr = [objcol._v_parent._v_nested_descr[objcol._v_pos]]
        # Try to convert the column object into a NumPy ndarray
        try:
            # If the column is a recarray (or kind of), convert into ndarray
            if hasattr(column, "dtype") and column.dtype.kind == "V":
                column = np.rec.array(column, dtype=descr).field(0)
            else:
                # Make sure the result is always a *copy* of the original,
                # so the resulting object is safe for in-place conversion.
                iflavor = flavor_of(column)
                column = array_as_internal(column, iflavor)
        except Exception as exc:  # XXX
            raise ValueError(
                f"column parameter cannot be converted into a "
                f"ndarray object compliant with specified column "
                f"'{column}'. The error was: <{exc}>"
            )

        # Get rid of single-dimensional dimensions
        column = column.squeeze()
        if column.shape == ():
            # Oops, stripped off too much dimensions
            column.shape = (1,)

        if stop is None:
            # compute the stop value. start + len(rows)*step does not work
            stop = start + (len(column) - 1) * step + 1
        (start, stop, step) = self._process_range(start, stop, step)
        if stop > self.nrows:
            raise IndexError(
                "This modification will exceed the length of "
                "the table. Giving up."
            )
        # Compute the number of rows to read.
        nrows = len(range(start, stop, step))
        if len(column) < nrows:
            raise ValueError(
                "The value has not enough elements to fill-in "
                "the specified range"
            )
        # Now, read the original values:
        mod_recarr = self._read(start, stop, step)
        # Modify the appropriate column in the original recarray
        mod_col = get_nested_field(mod_recarr, colname)
        mod_col[:] = column
        # save this modified rows in table
        self._update_records(start, stop, step, mod_recarr)
        # Redo the index if needed
        self._reindex([colname])

        return SizeType(nrows)

    def modify_columns(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
        columns: Sequence | None = None,
        names: list[str] | tuple[str, ...] | None = None,
    ) -> int:
        """Modify a series of columns in the row slice [start:stop:step].

        The names argument specifies the names of the columns in the
        table to be modified with the data given in columns.  This
        method returns the number of rows modified.  Should the
        modification exceed the length of the table, an IndexError
        is raised before changing data.

        The columns argument may be any object which can be converted
        to a structured array compliant with the structure of the
        columns to be modified (otherwise, a ValueError is raised).
        This includes NumPy structured arrays, lists of tuples or array
        records, and a string or Python buffer.

        """
        if step is None:
            step = 1
        if type(names) not in (list, tuple):
            raise TypeError("The 'names' parameter must be a list of strings.")

        if columns is None:  # Nothing to be done
            return SizeType(0)
        if start is None:
            start = 0
        if start < 0:
            raise ValueError("'start' must have a positive value.")
        if step < 1:
            raise ValueError(
                "'step' must have a value greater or equal than 1."
            )
        descr = []
        for colname in names:
            objcol = self._get_column_instance(colname)
            descr.append(objcol._v_parent._v_nested_descr[objcol._v_pos])
            # descr.append(objcol._v_parent._v_dtype[objcol._v_pos])
        # Try to convert the columns object into a recarray
        try:
            # Make sure the result is always a *copy* of the original,
            # so the resulting object is safe for in-place conversion.
            iflavor = flavor_of(columns)
            if iflavor != "python":
                columns = array_as_internal(columns, iflavor)
                recarray = np.rec.array(columns, dtype=descr)
            else:
                recarray = np.rec.fromarrays(columns, dtype=descr)
        except Exception as exc:  # XXX
            raise ValueError(
                f"columns parameter cannot be converted into a "
                f"recarray object compliant with table '{self}'. "
                f"The error was: <{exc}>"
            )

        if stop is None:
            # compute the stop value. start + len(rows)*step does not work
            stop = start + (len(recarray) - 1) * step + 1
        (start, stop, step) = self._process_range(start, stop, step)
        if stop > self.nrows:
            raise IndexError(
                "This modification will exceed the length of "
                "the table. Giving up."
            )
        # Compute the number of rows to read.
        nrows = len(range(start, stop, step))
        if len(recarray) < nrows:
            raise ValueError(
                "The value has not enough elements to fill-in "
                "the specified range"
            )
        # Now, read the original values:
        mod_recarr = self._read(start, stop, step)
        # Modify the appropriate columns in the original recarray
        for i, name in enumerate(recarray.dtype.names):
            mod_col = get_nested_field(mod_recarr, names[i])
            mod_col[:] = recarray[name].squeeze()
        # save this modified rows in table
        self._update_records(start, stop, step, mod_recarr)
        # Redo the index if needed
        self._reindex(names)

        return SizeType(nrows)

    def flush_rows_to_index(self, _lastrow: bool = True) -> int:
        """Add remaining rows in buffers to non-dirty indexes.

        This can be useful when you have chosen non-automatic indexing
        for the table (see the :attr:`Table.autoindex` property in
        :class:`Table`) and you want to update the indexes on it.

        """
        rowsadded = 0
        if self.indexed:
            # Update the number of unsaved indexed rows
            start = self._indexedrows
            nrows = self._unsaved_indexedrows
            for colname, colindexed in self.colindexed.items():
                if colindexed:
                    col = self.cols._g_col(colname)
                    if nrows > 0 and not col.index.dirty:
                        rowsadded = self._add_rows_to_index(
                            colname, start, nrows, _lastrow, update=True
                        )
            self._unsaved_indexedrows -= rowsadded
            self._indexedrows += rowsadded
        return rowsadded

    def _add_rows_to_index(
        self, colname: str, start: int, nrows: int, lastrow: bool, update: bool
    ) -> int:
        """Add more elements to the existing index."""
        # This method really belongs to Column, but since it makes extensive
        # use of the table, it gets dangerous when closing the file, since the
        # column may be accessing a table which is being destroyed.
        index = self.cols._g_col(colname).index
        slicesize = index.slicesize
        # The next loop does not rely on xrange so that it can
        # deal with long ints (i.e. more than 32-bit integers)
        # This allows to index columns with more than 2**31 rows
        # F. Alted 2005-05-09
        start_lr = index.sorted.nrows * slicesize
        indexedrows = start_lr - start
        stop = start + nrows - slicesize + 1
        while start_lr < stop:
            index.append(
                [self._read(start_lr, start_lr + slicesize, 1, colname)],
                update=update,
            )
            indexedrows += slicesize
            start_lr += slicesize
        # index the remaining rows in last row
        if lastrow and start_lr < self.nrows:
            index.append_last_row(
                [self._read(start_lr, self.nrows, 1, colname)], update=update
            )
            indexedrows += self.nrows - start_lr
        return indexedrows

    def remove_rows(
        self,
        start: int | None = None,
        stop: int | None = None,
        step: int | None = None,
    ) -> int:
        """Remove a range of rows in the table.

        If only start is supplied, that row and all following will be deleted.
        If a range is supplied, i.e. both the start and stop parameters are
        passed, all the rows in the range are removed.

        .. versionchanged:: 3.0
           The start, stop and step parameters now behave like in slice.

        .. seealso:: remove_row()

        Parameters
        ----------
        start : int
            Sets the starting row to be removed. It accepts negative values
            meaning that the count starts from the end.  A value of 0 means the
            first row.
        stop : int
            Sets the last row to be removed to stop-1, i.e. the end point is
            omitted (in the Python range() tradition). Negative values are also
            accepted. If None all rows after start will be removed.
        step : int
            The step size between rows to remove.

            .. versionadded:: 3.0

        Examples
        --------
        Removing rows from 5 to 10 (excluded)::

            t.remove_rows(5, 10)

        Removing all rows starting from the 10th::

            t.remove_rows(10)

        Removing the 6th row::

            t.remove_rows(6, 7)

        .. note::

            removing a single row can be done using the specific
            :meth:`remove_row` method.

        """
        (start, stop, step) = self._process_range(start, stop, step)
        nrows = self._remove_rows(start, stop, step)
        # remove_rows is an invalidating index operation
        self._reindex(self.colpathnames)

        return SizeType(nrows)

    def remove_row(self, n: int) -> None:
        """Remove a row from the table.

        Parameters
        ----------
        n : int
            The index of the row to remove.


        .. versionadded:: 3.0

        Examples
        --------
        Remove row 15::

            table.remove_row(15)

        Which is equivalent to::

            table.remove_rows(15, 16)

        .. warning::

            This is not equivalent to::

                table.remove_rows(15)

        """
        self.remove_rows(start=n, stop=n + 1)

    def _g_update_dependent(self) -> None:
        super()._g_update_dependent()

        # Update the new path in columns
        self.cols._g_update_table_location(self)

        # Update the new path in the Row instance, if cached.  Fixes #224.
        if "row" in self.__dict__:
            self.__dict__["row"] = tableextension.Row(self)

    def _g_move(self, newparent: Group, newname: str) -> None:
        """Move this node in the hierarchy.

        This overloads the Node._g_move() method.

        """
        itgpathname = _index_pathname_of(self)

        # First, move the table to the new location.
        super()._g_move(newparent, newname)

        # Then move the associated index group (if any).
        try:
            itgroup = self._v_file._get_node(itgpathname)
        except NoSuchNodeError:
            pass
        else:
            newigroup = self._v_parent
            newiname = _index_name_of(self)
            itgroup._g_move(newigroup, newiname)

    def _g_remove(self, recursive: bool = False, force: bool = False) -> None:
        # Remove the associated index group (if any).
        itgpathname = _index_pathname_of(self)
        try:
            itgroup = self._v_file._get_node(itgpathname)
        except NoSuchNodeError:
            pass
        else:
            itgroup._f_remove(recursive=True)
            self.indexed = False  # there are indexes no more

        # Remove the leaf itself from the hierarchy.
        super()._g_remove(recursive, force)

    def _set_column_indexing(self, colpathname: str, indexed: bool) -> None:
        """Mark the referred column as indexed or non-indexed."""
        colindexed = self.colindexed
        isindexed, wasindexed = bool(indexed), colindexed[colpathname]
        if isindexed == wasindexed:
            return  # indexing state is unchanged

        # Changing the set of indexed columns invalidates the condition cache
        self._condition_cache.clear()
        colindexed[colpathname] = isindexed
        self.indexed = max(colindexed.values())  # this is an OR :)

    def _mark_columns_as_dirty(self, colnames: list[str]) -> None:
        """Mark column indexes in `colnames` as dirty."""
        assert len(colnames) > 0
        if self.indexed:
            colindexed, cols = self.colindexed, self.cols
            # Mark the proper indexes as dirty
            for colname in colnames:
                if colindexed[colname]:
                    col = cols._g_col(colname)
                    col.index.dirty = True

    def _reindex(self, colnames: list[str]) -> None:
        """Re-index columns in `colnames` if automatic indexing is true."""
        if self.indexed:
            colindexed, cols = self.colindexed, self.cols
            colstoindex = []
            # Mark the proper indexes as dirty
            for colname in colnames:
                if colindexed[colname]:
                    col = cols._g_col(colname)
                    col.index.dirty = True
                    colstoindex.append(colname)
            # Now, re-index the dirty ones
            if self.autoindex and colstoindex:
                self._do_reindex(dirty=True)
            # The table caches for indexed queries are dirty now
            self._dirtycache = True

    def _do_reindex(self, dirty: bool) -> int:
        """Execute common code for `reindex()` and `reindex_dirty()`."""
        indexedrows = 0
        for colname, colindexed in self.colindexed.items():
            if colindexed:
                indexcol = self.cols._g_col(colname)
                indexedrows = indexcol._do_reindex(dirty)
        # Update counters in case some column has been updated
        if indexedrows > 0:
            self._indexedrows = indexedrows
            self._unsaved_indexedrows = self.nrows - indexedrows

        return SizeType(indexedrows)

    def reindex(self) -> None:
        """Recompute all the existing indexes in the table.

        This can be useful when you suspect that, for any reason, the
        index information for columns is no longer valid and want to
        rebuild the indexes on it.

        """
        self._do_reindex(dirty=False)

    def reindex_dirty(self) -> None:
        """Recompute the existing indexes in table, *if* they are dirty.

        This can be useful when you have set :attr:`Table.autoindex`
        (see :class:`Table`) to false for the table and you want to
        update the indexes after an invalidating index operation
        (:meth:`Table.remove_rows`, for example).

        """
        self._do_reindex(dirty=True)

    def _g_copy_rows(
        self,
        obj: Table,
        start: int,
        stop: int,
        step: int,
        sortby: Column | str | None,
        checkCSI: bool,  # noqa: N803
    ) -> None:
        """Copy rows from self to object."""
        if sortby is None:
            self._g_copy_rows_optim(obj, start, stop, step)
            return
        lenbuf = self.nrowsinbuf
        absstep = step
        if step < 0:
            absstep = -step
            start, stop = stop + 1, start + 1
        if sortby is not None:
            index = self._check_sortby_csi(sortby, checkCSI)
        for start2 in range(start, stop, absstep * lenbuf):
            stop2 = start2 + absstep * lenbuf
            if stop2 > stop:
                stop2 = stop
            # The next 'if' is not needed, but it doesn't bother either
            if sortby is None:
                rows = self[start2:stop2:step]
            else:
                coords = index[start2:stop2:step]
                rows = self.read_coordinates(coords)
            # Save the records on disk
            obj.append(rows)
        obj.flush()

    def _g_copy_rows_optim(
        self, obj: Table, start: int, stop: int, step: int
    ) -> None:
        """Copy rows from self to object (optimized version)."""
        nrowsinbuf = self.nrowsinbuf
        obj._open_append(self._v_iobuf)
        nrowsdest = obj.nrows
        for start2 in range(start, stop, step * nrowsinbuf):
            # Save the records on disk
            stop2 = start2 + step * nrowsinbuf
            if stop2 > stop:
                stop2 = stop
            # Optimized version (it saves some conversions)
            nrows = ((stop2 - start2 - 1) // step) + 1
            self.row._fill_col(self._v_iobuf, start2, stop2, step, None)
            # The output buffer is created anew,
            # so the operation is safe to in-place conversion.
            obj._append_records(nrows)
            nrowsdest += nrows
        obj._close_append()

    def _g_prop_indexes(self, other: Table) -> None:
        """Generate index in `other` table for every indexed column here."""
        oldcols, newcols = self.colinstances, other.colinstances
        for colname in newcols:
            if isinstance(oldcols[colname], Column):
                oldcolindexed = oldcols[colname].is_indexed
                if oldcolindexed:
                    oldcolindex = oldcols[colname].index
                    newcol = newcols[colname]
                    newcol.create_index(
                        kind=oldcolindex.kind,
                        optlevel=oldcolindex.optlevel,
                        filters=oldcolindex.filters,
                        tmp_dir=None,
                    )

    def _g_copy_with_stats(
        self,
        group: Group,
        name: str,
        start: int,
        stop: int,
        step: int,
        title: str,
        filters: Filters | None,
        chunkshape: int | tuple[int] | None,
        _log: bool,
        **kwargs,
    ) -> tuple[Table, int]:
        """Private part of Leaf.copy() for each kind of leaf."""
        # Get the private args for the Table flavor of copy()
        sortby = kwargs.pop("sortby", None)
        propindexes = kwargs.pop("propindexes", False)
        check_csi = kwargs.pop("checkCSI", False)
        # Compute the correct indices.
        (start, stop, step) = self._process_range_read(
            start, stop, step, warn_negstep=sortby is None
        )
        # And the number of final rows
        nrows = len(range(start, stop, step))
        # Create the new table and copy the selected data.
        newtable = Table(
            group,
            name,
            self.description,
            title=title,
            filters=filters,
            expectedrows=nrows,
            chunkshape=chunkshape,
            _log=_log,
        )
        self._g_copy_rows(newtable, start, stop, step, sortby, check_csi)
        nbytes = newtable.nrows * newtable.rowsize
        # Generate equivalent indexes in the new table, if required.
        if propindexes and self.indexed:
            self._g_prop_indexes(newtable)
        return (newtable, nbytes)

    # This overloading of copy is needed here in order to document
    # the additional keywords for the Table case.
    def copy(
        self,
        newparent: Group | None = None,
        newname: str | None = None,
        overwrite: bool = False,
        createparents: bool = False,
        **kwargs,
    ) -> Table:
        """Copy this table and return the new one.

        This method has the behavior and keywords described in
        :meth:`Leaf.copy`.  Moreover, it recognises the following additional
        keyword arguments.

        Parameters
        ----------
        sortby
            If specified, and sortby corresponds to a column with an index,
            then the copy will be sorted by this index.  If you want to ensure
            a fully sorted order, the index must be a CSI one.  A reverse
            sorted copy can be achieved by specifying a negative value for the
            step keyword.  If sortby is omitted or None, the original table
            order is used.
        checkCSI
            If true and a CSI index does not exist for the sortby column, an
            error will be raised.  If false (the default), it does nothing.
            You can use this flag in order to explicitly check for the
            existence of a CSI index.
        propindexes
            If true, the existing indexes in the source table are propagated
            (created) to the new one.  If false (the default), the indexes are
            not propagated.

        """
        return super().copy(
            newparent, newname, overwrite, createparents, **kwargs
        )

    def flush(self) -> None:
        """Flush the table buffers."""
        if self._v_file._iswritable():
            # Flush rows that remains to be appended
            if "row" in self.__dict__:
                self.row._flush_buffered_rows()
            if self.indexed and self.autoindex:
                # Flush any unindexed row
                rowsadded = self.flush_rows_to_index(_lastrow=True)
                assert rowsadded <= 0 or self._indexedrows == self.nrows, (
                    "internal error: the number of indexed rows (%d) "
                    "and rows in the table (%d) is not equal; "
                    "please report this to the authors."
                    % (self._indexedrows, self.nrows)
                )
                if self._dirtyindexes:
                    # Finally, re-index any dirty column
                    self.reindex_dirty()

        super().flush()

    def _g_pre_kill_hook(self) -> None:
        """Code to be called before killing the node."""
        # Flush the buffers before to clean-up them
        # self.flush()
        # It seems that flushing during the __del__ phase is a sure receipt for
        # bringing all kind of problems:
        # 1. Illegal Instruction
        # 2. Malloc(): trying to call free() twice
        # 3. Bus Error
        # 4. Segmentation fault
        # So, the best would be doing *nothing* at all in this __del__ phase.
        # As a consequence, the I/O will not be cleaned until a call to
        # Table.flush() would be done. This could lead to a potentially large
        # memory consumption.
        # NOTE: The user should make a call to Table.flush() whenever he has
        #       finished working with his table.
        # I've added a Performance warning in order to compel the user to
        # call self.flush() before the table is being preempted.
        # F. Alted 2006-08-03
        if ("row" in self.__dict__ and self.row._get_unsaved_nrows() > 0) or (
            self.indexed
            and self.autoindex
            and (self._unsaved_indexedrows > 0 or self._dirtyindexes)
        ):
            warnings.warn(
                f"table ``{self._v_pathname}`` is being preempted from "
                f"alive nodes without its buffers being flushed or with some "
                f"index being dirty.  This may lead to very "
                f"ineficient use of resources and even to fatal "
                f"errors in certain situations.  Please do a call "
                f"to the .flush() or .reindex_dirty() methods on "
                f"this table before start using other nodes.",
                PerformanceWarning,
            )
        # Get rid of the IO buffers (if they have been created at all)
        mydict = self.__dict__
        if "_v_iobuf" in mydict:
            del mydict["_v_iobuf"]
        if "_v_wdflts" in mydict:
            del mydict["_v_wdflts"]

    def _f_close(self, flush: bool = True) -> None:
        if not self._v_isopen:
            return  # the node is already closed

        # .. note::
        #
        #   As long as ``Table`` objects access their indices on closing,
        #   ``File.close()`` will need to make *two separate passes*
        #   to first close ``Table`` objects and then ``Index`` hierarchies.
        #

        # Flush right now so the row object does not get in the middle.
        if flush:
            self.flush()

        # Some warnings can be issued after calling `self._g_set_location()`
        # in `self.__init__()`.  If warnings are turned into exceptions,
        # `self._g_post_init_hook` may not be called and `self.cols` not set.
        # One example of this is
        # ``test_create.createTestCase.test05_maxFieldsExceeded()``.
        cols = self.cols
        if cols is not None:
            cols._g_close()

        # Clean address cache
        self._clean_chunk_addrs()

        # Close myself as a leaf.
        super()._f_close(False)

    def __repr__(self) -> str:
        """Return the string representation of `Table` objects.

        This provides column metainfo in addition to standard __str__.
        """
        if self.indexed:
            return f"""\
{self}
  description := {self.description!r}
  byteorder := {self.byteorder!r}
  chunkshape := {self.chunkshape!r}
  autoindex := {self.autoindex!r}
  colindexes := {_ColIndexes(self.colindexes)!r}"""
        else:
            return f"""\
{self}
  description := {self.description!r}
  byteorder := {self.byteorder!r}
  chunkshape := {self.chunkshape!r}"""


class Cols:
    """Container for columns in a table or nested column.

    This class is used as an *accessor* to the columns in a table or nested
    column.  It supports the *natural naming* convention, so that you can
    access the different columns as attributes which lead to Column instances
    (for non-nested columns) or other Cols instances (for nested columns).

    For instance, if table.cols is a Cols instance with a column named col1
    under it, the later can be accessed as table.cols.col1. If col1 is nested
    and contains a col2 column, this can be accessed as table.cols.col1.col2
    and so on. Because of natural naming, the names of members start with
    special prefixes, like in the Group class (see :ref:`GroupClassDescr`).

    Like the Column class (see :ref:`ColumnClassDescr`), Cols supports item
    access to read and write ranges of values in the table or nested column.


    .. rubric:: Cols attributes

    .. attribute:: _v_colnames

        A list of the names of the columns hanging directly
        from the associated table or nested column.  The order of
        the names matches the order of their respective columns in
        the containing table.

    .. attribute:: _v_colpathnames

        A list of the pathnames of all the columns under the
        associated table or nested column (in preorder).  If it does
        not contain nested columns, this is exactly the same as the
        :attr:`Cols._v_colnames` attribute.

    .. attribute:: _v_desc

        The associated Description instance (see
        :ref:`DescriptionClassDescr`).

    """

    @property
    def _v_table(self) -> Table:
        """Return the parent Table instance (see :ref:`TableClassDescr`)."""
        return self._v__tableFile._get_node(self._v__tablePath)

    def __init__(self, table: Table, desc: Description) -> None:
        dict_ = self.__dict__
        dict_["_v__tableFile"] = table._v_file
        dict_["_v__tablePath"] = table._v_pathname
        dict_["_v_desc"] = desc
        dict_["_v_colnames"] = desc._v_names
        dict_["_v_colpathnames"] = table.description._v_pathnames
        # Put the column in the local dictionary
        for name in desc._v_names:
            if name in desc._v_types:
                dict_[name] = Column(table, name, desc)
            else:
                dict_[name] = Cols(table, desc._v_colobjects[name])

    def _g_update_table_location(self, table: Table) -> None:
        """Update the location information about the associated `table`."""
        dict_ = self.__dict__
        dict_["_v__tableFile"] = table._v_file
        dict_["_v__tablePath"] = table._v_pathname

        # Update the locations in individual columns.
        for colname in self._v_colnames:
            dict_[colname]._g_update_table_location(table)

    def __len__(self) -> int:
        """Get the number of top level columns in table."""
        return len(self._v_colnames)

    def _f_col(self, colname: str) -> Cols:
        """Get an accessor to the column colname.

        This method returns a Column instance (see :ref:`ColumnClassDescr`) if
        the requested column is not nested, and a Cols instance (see
        :ref:`ColsClassDescr`) if it is.  You may use full column pathnames in
        colname.

        Calling cols._f_col('col1/col2') is equivalent to using cols.col1.col2.
        However, the first syntax is more intended for programmatic use.  It is
        also better if you want to access columns with names that are not valid
        Python identifiers.

        """
        if not isinstance(colname, str):
            raise TypeError(
                f"Parameter can only be an string. You passed object: "
                f"{colname}"
            )
        if (
            colname.find("/") > -1 and colname not in self._v_colpathnames
        ) and colname not in self._v_colnames:
            raise KeyError(
                f"Cols accessor "
                f"``{self._v__tablePath}.cols{self._v_desc._v_pathname}`` "
                f"does not have a column named ``{colname}``"
            )

        return self._g_col(colname)

    def _g_col(self, colname: str) -> Cols:
        """Like `self._f_col()` but it does not check arguments."""
        # Get the Column or Description object
        inames = colname.split("/")
        cols = self
        for iname in inames:
            cols = cols.__dict__[iname]
        return cols

    def __getitem__(self, key: int | slice) -> Any:
        """Get a row or a range of rows from a table or nested column.

        If key argument is an integer, the corresponding nested type row is
        returned as a record of the current flavor. If key is a slice, the
        range of rows determined by it is returned as a structured array of the
        current flavor.

        Examples
        --------
        ::

            record = table.cols[4]  # equivalent to table[4]
            recarray = table.cols.Info[4:1000:2]

        Those statements are equivalent to::

            nrecord = table.read(start=4)[0]
            nrecarray = table.read(start=4, stop=1000, step=2).field('Info')

        Here you can see how a mix of natural naming, indexing and slicing can
        be used as shorthands for the :meth:`Table.read` method.

        """
        table = self._v_table
        nrows = table.nrows
        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += nrows
            (start, stop, step) = table._process_range(key, key + 1, 1)
            colgroup = self._v_desc._v_pathname
            if colgroup == "":  # The root group
                return table.read(start, stop, step)[0]
            else:
                crecord = table.read(start, stop, step)[0]
                return crecord[colgroup]
        elif isinstance(key, slice):
            (start, stop, step) = table._process_range(
                key.start, key.stop, key.step
            )
            colgroup = self._v_desc._v_pathname
            if colgroup == "":  # The root group
                return table.read(start, stop, step)
            else:
                crecarray = table.read(start, stop, step)
                if hasattr(crecarray, "field"):
                    return crecarray.field(colgroup)  # RecArray case
                else:
                    return get_nested_field(crecarray, colgroup)  # numpy case
        else:
            raise TypeError(f"invalid index or slice: {key!r}")

    def __setitem__(self, key: int | slice, value: Any) -> None:
        """Set a row or a range of rows in a table or nested column.

        If key argument is an integer, the corresponding row is set to
        value. If key is a slice, the range of rows determined by it is set to
        value.

        Examples
        --------
        ::

            table.cols[4] = record
            table.cols.Info[4:1000:2] = recarray

        Those statements are equivalent to::

            table.modify_rows(4, rows=record)
            table.modify_column(4, 1000, 2, colname='Info', column=recarray)

        Here you can see how a mix of natural naming, indexing and slicing
        can be used as shorthands for the :meth:`Table.modify_rows` and
        :meth:`Table.modify_column` methods.

        """
        table = self._v_table
        nrows = table.nrows
        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += nrows
            (start, stop, step) = table._process_range(key, key + 1, 1)
        elif isinstance(key, slice):
            (start, stop, step) = table._process_range(
                key.start, key.stop, key.step
            )
        else:
            raise TypeError(f"invalid index or slice: {key!r}")

        # Actually modify the correct columns
        colgroup = self._v_desc._v_pathname
        if colgroup == "":  # The root group
            table.modify_rows(start, stop, step, rows=value)
        else:
            table.modify_column(
                start, stop, step, colname=colgroup, column=value
            )

    def _g_close(self) -> None:
        # First, close the columns (ie possible indices open)
        for col in self._v_colnames:
            colobj = self._g_col(col)
            if isinstance(colobj, Column):
                colobj.close()
                # Delete the reference to column
                del self.__dict__[col]
            else:
                colobj._g_close()

        self.__dict__.clear()

    def __str__(self) -> str:
        """Return the string representation for this object."""
        # The pathname
        descpathname = self._v_desc._v_pathname
        if descpathname:
            descpathname = "." + descpathname
        return (
            f"{self._v__tablePath}.cols{descpathname} "
            f"({self.__class__.__name__}), "
            f"{len(self._v_colnames)} columns"
        )

    def __repr__(self) -> str:
        """Detailed string representation for this object."""
        lines = [f"{self!s}"]
        for name in self._v_colnames:
            # Get this class name
            classname = getattr(self, name).__class__.__name__
            # The type
            if name in self._v_desc._v_dtypes:
                tcol = self._v_desc._v_dtypes[name]
                # The shape for this column
                shape = (self._v_table.nrows,) + self._v_desc._v_dtypes[
                    name
                ].shape
            else:
                tcol = "Description"
                # Description doesn't have a shape currently
                shape = ()
            lines.append(f"  {name} ({classname}{shape}, {tcol})")
        return "\n".join(lines) + "\n"


class Column:
    """Accessor for a non-nested column in a table.

    Each instance of this class is associated with one *non-nested* column of a
    table. These instances are mainly used to read and write data from the
    table columns using item access (like the Cols class - see
    :ref:`ColsClassDescr`), but there are a few other associated methods to
    deal with indexes.

    .. rubric:: Column attributes

    .. attribute:: descr

        The Description (see :ref:`DescriptionClassDescr`) instance of the
        parent table or nested column.

    .. attribute:: name

        The name of the associated column.

    .. attribute:: pathname

        The complete pathname of the associated column (the same as
        Column.name if the column is not inside a nested column).

    .. attribute:: attrs

        Column attributes (see :ref:`ColClassDescr`).

    Parameters
    ----------
    table
        The parent table instance
    name
        The name of the column that is associated with this object
    descr
        The parent description object

    """

    @lazyattr
    def dtype(self) -> np.dtype:
        """Return the NumPy dtype that most closely matches this column."""
        return self.descr._v_dtypes[self.name].base  # Get rid of shape info

    @lazyattr
    def type(self) -> str:  # noqa: A003
        """Return the PyTables type of the column (a string)."""
        return self.descr._v_types[self.name]

    @property
    def table(self) -> Table:
        """Return the parent Table instance (see :ref:`TableClassDescr`)."""
        return self._table_file._get_node(self._table_path)

    @property
    def index(self) -> Index | None:
        """Return the Index instance associated with this column.

        Return `None` if the column is not indexed.

        See :ref:`IndexClassDescr`.
        """
        index_path = _index_pathname_of_column_(
            self._table_path, self.pathname
        )
        try:
            index = self._table_file._get_node(index_path)
        except NodeError:
            index = None  # The column is not indexed
        return index

    @lazyattr
    def _itemtype(self) -> np.dtype:
        return self.descr._v_dtypes[self.name]

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the shape of this column."""
        return (self.table.nrows,) + self.descr._v_dtypes[self.name].shape

    @property
    def is_indexed(self) -> bool:
        """Return True if the column is indexed, false otherwise."""
        if self.index is None:
            return False
        else:
            return True

    @property
    def maindim(self) -> int:
        """Return the dimension along which iterators work.

        Its value is 0 (i.e. the first dimension).
        """
        return 0

    def __init__(self, table: Table, name: str, descr: Description) -> None:
        self._table_file = table._v_file
        self._table_path = table._v_pathname
        self.name = name
        """The name of the associated column."""
        self.pathname = descr._v_colobjects[name]._v_pathname
        """The complete pathname of the associated column (the same as
        Column.name if the column is not inside a nested column)."""
        self.descr = descr
        """The Description (see :ref:`DescriptionClassDescr`) instance of the
        parent table or nested column."""
        self._v_attrs = ColumnAttributeSet(self)

    def _g_update_table_location(self, table: Table) -> None:
        """Update the location information about the associated `table`."""
        self._table_file = table._v_file
        self._table_path = table._v_pathname

    def __len__(self) -> int:
        """Get the number of elements in the column.

        This matches the length in rows of the parent table.

        """
        return self.table.nrows

    def __getitem__(self, key: int | slice) -> np.ndarray:
        """Get a row or a range of rows from a column.

        If key argument is an integer, the corresponding element in the column
        is returned as an object of the current flavor.  If key is a slice, the
        range of elements determined by it is returned as an array of the
        current flavor.

        Examples
        --------
        ::

            print("Column handlers:")
            for name in table.colnames:
                print(table.cols._f_col(name))
                print("Select table.cols.name[1]-->", table.cols.name[1])
                print("Select table.cols.name[1:2]-->", table.cols.name[1:2])
                print("Select table.cols.name[:]-->", table.cols.name[:])
                print("Select table.cols._f_col('name')[:]-->",
                                                table.cols._f_col('name')[:])

        The output of this for a certain arbitrary table is::

            Column handlers:
            /table.cols.name (Column(), string, idx=None)
            /table.cols.lati (Column(), int32, idx=None)
            /table.cols.longi (Column(), int32, idx=None)
            /table.cols.vector (Column(2,), int32, idx=None)
            /table.cols.matrix2D (Column(2, 2), float64, idx=None)
            Select table.cols.name[1]--> Particle:     11
            Select table.cols.name[1:2]--> ['Particle:     11']
            Select table.cols.name[:]--> ['Particle:     10'
             'Particle:     11' 'Particle:     12'
             'Particle:     13' 'Particle:     14']
            Select table.cols._f_col('name')[:]--> ['Particle:     10'
             'Particle:     11' 'Particle:     12'
             'Particle:     13' 'Particle:     14']

        See the :file:`examples/table2.py` file for a more complete example.

        """
        table = self.table

        # Generalized key support not there yet, but at least allow
        # for a tuple with one single element (the main dimension).
        # (key,) --> key
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]

        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= table.nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += table.nrows
            (start, stop, step) = table._process_range(key, key + 1, 1)
            return table.read(start, stop, step, self.pathname)[0]
        elif isinstance(key, slice):
            (start, stop, step) = table._process_range(
                key.start, key.stop, key.step
            )
            return table.read(start, stop, step, self.pathname)
        else:
            raise TypeError("'%s' key type is not valid in this context" % key)

    def __iter__(self) -> Generator[np.ndarray]:
        """Iterate through all items in the column."""
        table = self.table
        itemsize = self.dtype.itemsize
        nrowsinbuf = table._v_file.params["IO_BUFFER_SIZE"] // itemsize
        buf = np.empty((nrowsinbuf,), self._itemtype)
        max_row = len(self)
        for start_row in range(0, len(self), nrowsinbuf):
            end_row = min(start_row + nrowsinbuf, max_row)
            buf_slice = buf[0 : end_row - start_row]
            table.read(
                start_row, end_row, 1, field=self.pathname, out=buf_slice
            )
            yield from buf_slice

    def __setitem__(self, key: int | slice, value: Any) -> int:
        """Set a row or a range of rows in a column.

        If key argument is an integer, the corresponding element is set to
        value.  If key is a slice, the range of elements determined by it is
        set to value.

        Examples
        --------
        ::

            # Modify row 1
            table.cols.col1[1] = -1

            # Modify rows 1 and 3
            table.cols.col1[1::2] = [2,3]

        Which is equivalent to::

            # Modify row 1
            table.modify_columns(start=1, columns=[[-1]], names=['col1'])

            # Modify rows 1 and 3
            columns = np.rec.fromarrays([[2,3]], formats='i4')
            table.modify_columns(start=1, step=2, columns=columns,
                                 names=['col1'])

        """
        table = self.table
        table._v_file._check_writable()

        # Generalized key support not there yet, but at least allow
        # for a tuple with one single element (the main dimension).
        # (key,) --> key
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]

        if is_idx(key):
            key = operator.index(key)

            # Index out of range protection
            if key >= table.nrows:
                raise IndexError("Index out of range")
            if key < 0:
                # To support negative values
                key += table.nrows
            return table.modify_column(
                key, key + 1, 1, [[value]], self.pathname
            )
        elif isinstance(key, slice):
            (start, stop, step) = table._process_range(
                key.start, key.stop, key.step
            )
            return table.modify_column(start, stop, step, value, self.pathname)
        else:
            raise ValueError("Non-valid index or slice: %s" % key)

    def create_index(
        self,
        optlevel: int = 6,
        kind: str = "medium",
        filters: Filters | None = None,
        tmp_dir: str | None = None,
        _blocksizes: tuple[int, int, int, int] | None = None,
        _testmode: bool = False,
        _verbose: bool = False,
    ) -> int:
        """Create an index for this column.

        .. warning::

            In some situations it is useful to get a completely sorted
            index (CSI).  For those cases, it is best to use the
            :meth:`Column.create_csindex` method instead.

        Parameters
        ----------
        optlevel : int
            The optimization level for building the index.  The levels range
            from 0 (no optimization) up to 9 (maximum optimization).  Higher
            levels of optimization mean better chances for reducing the entropy
            of the index at the price of using more CPU, memory and I/O
            resources for creating the index.
        kind : str
            The kind of the index to be built.  It can take the 'ultralight',
            'light', 'medium' or 'full' values.  Lighter kinds ('ultralight'
            and 'light') mean that the index takes less space on disk, but will
            perform queries slower.  Heavier kinds ('medium' and 'full') mean
            better chances for reducing the entropy of the index (increasing
            the query speed) at the price of using more disk space as well as
            more CPU, memory and I/O resources for creating the index.

            Note that selecting a full kind with an optlevel of 9 (the maximum)
            guarantees the creation of an index with zero entropy, that is, a
            completely sorted index (CSI) - provided that the number of rows in
            the table does not exceed the 2**48 figure (that is more than 100
            trillions of rows).  See :meth:`Column.create_csindex` method for a
            more direct way to create a CSI index.
        filters : Filters
            Specify the Filters instance used to compress the index.  If None,
            default index filters will be used (currently, zlib level 1 with
            shuffling).
        tmp_dir
            When kind is other than 'ultralight', a temporary file is created
            during the index build process.  You can use the tmp_dir argument
            to specify the directory for this temporary file.  The default is
            to create it in the same directory as the file containing the
            original table.

        """
        kinds = ["ultralight", "light", "medium", "full"]
        if kind not in kinds:
            raise ValueError("Kind must have any of these values: %s" % kinds)
        if not isinstance(optlevel, int) or (optlevel < 0 or optlevel > 9):
            raise ValueError(
                "Optimization level must be an integer in the range 0-9"
            )
        if filters is None:
            filters = default_index_filters
        if tmp_dir is None:
            tmp_dir = str(Path(self._table_file.filename).parent)
        else:
            if not Path(tmp_dir).is_dir():
                raise ValueError(
                    f"Temporary directory '{tmp_dir}' does not exist"
                )
        if _blocksizes is not None and (
            not isinstance(_blocksizes, tuple) or len(_blocksizes) != 4
        ):
            raise ValueError(
                "_blocksizes must be a tuple with exactly 4 elements"
            )
        idxrows = _column__create_index(
            self, optlevel, kind, filters, tmp_dir, _blocksizes, _verbose
        )
        return SizeType(idxrows)

    def create_csindex(
        self,
        filters: Filters | None = None,
        tmp_dir: str | None = None,
        _blocksizes: tuple[int, int, int, int] | None = None,
        _testmode: bool = False,
        _verbose: bool = False,
    ) -> int:
        """Create a completely sorted index (CSI) for this column.

        This method guarantees the creation of an index with zero entropy, that
        is, a completely sorted index (CSI) -- provided that the number of rows
        in the table does not exceed the 2**48 figure (that is more than 100
        trillions of rows).  A CSI index is needed for some table methods (like
        :meth:`Table.itersorted` or :meth:`Table.read_sorted`) in order to
        ensure completely sorted results.

        For the meaning of filters and tmp_dir arguments see
        :meth:`Column.create_index`.

        Notes
        -----
        This method is equivalent to
        Column.create_index(optlevel=9, kind='full', ...).

        """
        return self.create_index(
            kind="full",
            optlevel=9,
            filters=filters,
            tmp_dir=tmp_dir,
            _blocksizes=_blocksizes,
            _testmode=_testmode,
            _verbose=_verbose,
        )

    def _do_reindex(self, dirty: bool) -> int:
        """Execute common code for reindex() and reindex_dirty() codes."""
        index = self.index
        dodirty = True
        if dirty and not index.dirty:
            dodirty = False
        if index is not None and dodirty:
            self._table_file._check_writable()
            # Get the old index parameters
            kind = index.kind
            optlevel = index.optlevel
            filters = index.filters
            # We *need* to tell the index that it is going to be undirty.
            # This is needed here so as to unnail() the condition cache.
            index.dirty = False
            # Delete the existing Index
            index._f_remove()
            # Create a new Index with the previous parameters
            return SizeType(
                self.create_index(
                    kind=kind, optlevel=optlevel, filters=filters
                )
            )
        else:
            return SizeType(0)  # The column is not intended for indexing

    def reindex(self) -> None:
        """Recompute the index associated with this column.

        This can be useful when you suspect that, for any reason,
        the index information is no longer valid and you want to rebuild it.

        This method does nothing if the column is not indexed.

        """
        self._do_reindex(dirty=False)

    def reindex_dirty(self) -> None:
        """Recompute the associated index only if it is dirty.

        This can be useful when you have set :attr:`Table.autoindex` to false
        for the table and you want to update the column's index after an
        invalidating index operation (like :meth:`Table.remove_rows`).

        This method does nothing if the column is not indexed.

        """
        self._do_reindex(dirty=True)

    def remove_index(self) -> None:
        """Remove the index associated with this column.

        This method does nothing if the column is not indexed. The removed
        index can be created again by calling the :meth:`Column.create_index`
        method.

        """
        self._table_file._check_writable()

        # Remove the index if existing.
        if self.is_indexed:
            index = self.index
            index._f_remove()
            self.table._set_column_indexing(self.pathname, False)

    def close(self) -> None:
        """Close the column."""
        self.__dict__.clear()

    def __str__(self) -> str:
        """Return the string representation for this object."""
        return (
            f"{self._table_path}.cols.{self.pathname.replace('/', '.')} "
            f"({self.__class__.__name__}{self.shape}, "
            f"{self.descr._v_types[self.name]}, idx={self.index})"
        )

    def __repr__(self) -> str:
        """Detailed string representation for this object."""
        return str(self)

    @lazyattr
    def _v_pos(self) -> int:
        return self.descr._v_colobjects[self.name]._v_pos

    @lazyattr
    def _v_col_attrs(self) -> dict:
        return self.descr._v_colobjects[self.name]._v_col_attrs

    @property
    def attrs(self) -> ColumnAttributeSet:
        """Column attributes."""
        return self._v_attrs


class ColumnAttributeSet:
    """Column attribute set."""

    def __init__(self, column: Column) -> None:

        self.__dict__["_v_tableattrs"] = column.table.attrs
        self.__dict__["_v_fieldindex"] = column._v_pos
        self.__dict__["_v_column_reference"] = weakref.ref(column)

        # Check if this column has _v_col_attrs set and translate them into
        # the table attribute format
        for col_attr_key, col_attr_val in column._v_col_attrs.items():
            self.__setitem__(col_attr_key, col_attr_val)

    def issystemcolumnname(self, key: str) -> bool:
        """Check whether a key is a reserved attribute name."""
        return key in ["_v_tableattrs", "_v_fieldindex", "_v_column_reference"]

    def _prefix(self, string: str) -> str:
        """Prefix keys with a special pattern for storing table attributes."""
        field_index = self.__dict__["_v_fieldindex"]
        return "FIELD_%i_ATTR_%s" % (field_index, string)

    def __getattr__(self, key: str) -> Any:
        """Retrieve a PyTables attribute for this column."""
        if not self.issystemcolumnname(key):
            return getattr(self._v_tableattrs, self._prefix(key))
        else:
            return super().__getattr__(key)

    def __setattr__(self, key: str, val: Any) -> Any:
        """Set a PyTables attribute for this column."""
        if not self.issystemcolumnname(key):
            setattr(self._v_tableattrs, self._prefix(key), val)
        else:
            return super().__setattr__(key, val)

    def __getitem__(self, key: str) -> Any:
        """Dictionary-like interface for __getattr__."""
        if not self.issystemcolumnname(key):
            return self._v_tableattrs[self._prefix(key)]
        else:
            return self[key]

    def __setitem__(self, key: str, value: Any) -> None:
        """Dictionary-like interface for __setattr__."""
        if not self.issystemcolumnname(key):
            self._v_tableattrs[self._prefix(key)] = value
        else:
            self[key] = value

    def __delattr__(self, key: str) -> None:
        """Delete the attribute for this column."""
        if self.issystemcolumnname(key):
            raise TypeError("Deleting system attributes is prohibited")
        else:
            delattr(self._v_tableattrs, self._prefix(key))

    def __delitem__(self, key: str) -> None:
        """Dictionary-like interface for __delattr__."""
        if self.issystemcolumnname(key):
            raise TypeError("Deleting system attributes is prohibited")
        else:
            del self._v_tableattrs[self._prefix(key)]

    def _f_rename(self, oldattrname: str, newattrname: str) -> None:
        """Rename an attribute from oldattrname to newattrname."""
        if oldattrname == newattrname:
            # Do nothing
            return

        if self.issystemcolumnname(oldattrname):
            raise TypeError("Renaming system attributes is prohibited")

        # First, fetch the value of the oldattrname
        attrvalue = getattr(self, oldattrname)

        # Now, create the new attribute
        setattr(self, newattrname, attrvalue)

        # Finally, remove the old attribute
        delattr(self, oldattrname)

    def _f_copy(self, where: Column) -> None:
        """Copy attributes to another column."""
        # Is there a better way to do this?
        if not isinstance(where, Column):
            raise TypeError(f"destination object is not a column: {where!r}")

        for key in self.keys():
            where.attrs[key] = self[key]

    def keys(self) -> list[str]:
        """Return the list of attributes for this column."""
        col_prefix = self._prefix("")
        length = len(col_prefix)
        return [
            key[length:]
            for key in self._v_tableattrs._v_attrnames
            if key.startswith(col_prefix)
        ]

    def contains(self, key: str) -> bool:
        """Return whether a key is in the attribute set."""
        return key in self.keys()

    def __str__(self) -> str:
        """Return the string representation for this object."""
        pathname = self._v_tableattrs._v__nodepath
        classname = self._v_column_reference().__class__.__name__
        # self._v_tableattrs._v_node.__class__.__name__
        attrnumber = sum(1 for _ in self.keys())
        columnname = self._v_column_reference().name

        return (
            f"{pathname}.cols.{columnname}._v_attrs ({classname}), "
            f"{attrnumber} attributes"
        )

    def __repr__(self) -> str:
        """Detailed string representation for this object."""
        # print additional info only if there are attributes to show
        attrnames = self.keys()
        if attrnames:
            rep = [f"{attr} := {getattr(self, attr)!r}" for attr in attrnames]
            return f"{self!s}:\n   [" + ",\n    ".join(rep) + "]"
        else:
            return str(self)
