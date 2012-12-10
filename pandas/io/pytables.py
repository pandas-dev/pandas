"""
High level interface to PyTables for reading and writing pandas data structures
to disk
"""

# pylint: disable-msg=E1101,W0613,W0603

from datetime import datetime, date
import time
import re
import copy
import itertools

import numpy as np
from pandas import (
    Series, TimeSeries, DataFrame, Panel, Index, MultiIndex, Int64Index
)
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparsePanel
from pandas.sparse.array import BlockIndex, IntIndex
from pandas.tseries.api import PeriodIndex, DatetimeIndex
from pandas.core.common import adjoin
from pandas.core.algorithms import match, unique

from pandas.core.categorical import Factor
from pandas.core.common import _asarray_tuplesafe, _try_sort
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.reshape import block2d_to_block3d
import pandas.core.common as com
from pandas.tools.merge import concat

import pandas.lib as lib
import pandas.algos as algos
import pandas.tslib as tslib

from contextlib import contextmanager

# reading and writing the full object in one go
_TYPE_MAP = {
    Series: 'series',
    SparseSeries: 'sparse_series',
    TimeSeries: 'series',
    DataFrame: 'frame',
    SparseDataFrame: 'sparse_frame',
    Panel: 'wide',
    SparsePanel: 'sparse_panel'
}

_NAME_MAP = {
    'series': 'Series',
    'time_series': 'TimeSeries',
    'sparse_series': 'SparseSeries',
    'frame': 'DataFrame',
    'sparse_frame': 'SparseDataFrame',
    'frame_table': 'DataFrame (Table)',
    'wide': 'Panel',
    'sparse_panel': 'SparsePanel',
    'wide_table': 'Panel (Table)',
    'long': 'LongPanel',
    # legacy h5 files
    'Series': 'Series',
    'TimeSeries': 'TimeSeries',
    'DataFrame': 'DataFrame',
    'DataMatrix': 'DataMatrix'
}

# legacy handlers
_LEGACY_MAP = {
    'Series': 'legacy_series',
    'TimeSeries': 'legacy_series',
    'DataFrame': 'legacy_frame',
    'DataMatrix': 'legacy_frame',
    'WidePanel': 'wide_table',
}

# oh the troubles to reduce import time
_table_mod = None
_table_supports_index = False

def _tables():
    global _table_mod
    global _table_supports_index
    if _table_mod is None:
        import tables
        _table_mod = tables

        # version requirements
        major, minor, subv = tables.__version__.split('.')
        if int(major) >= 2 and int(minor[0]) >= 3:
            _table_supports_index = True

    return _table_mod

@contextmanager
def get_store(path, mode='a', complevel=None, complib=None,
              fletcher32=False):
    """
    Creates an HDFStore instance. This function can be used in a with statement

    Parameters
    ----------
    path : string
        File path to HDF5 file
    mode : {'a', 'w', 'r', 'r+'}, default 'a'

        ``'r'``
            Read-only; no data can be modified.
        ``'w'``
            Write; a new file is created (an existing file with the same
            name would be deleted).
        ``'a'``
            Append; an existing file is opened for reading and writing,
            and if the file does not exist it is created.
        ``'r+'``
            It is similar to ``'a'``, but the file must already exist.
    complevel : int, 1-9, default 0
            If a complib is specified compression will be applied
            where possible
    complib : {'zlib', 'bzip2', 'lzo', 'blosc', None}, default None
            If complevel is > 0 apply compression to objects written
            in the store wherever possible
    fletcher32 : bool, default False
            If applying compression use the fletcher32 checksum

    Examples
    --------
    >>> with get_store('test.h5') as store:
    >>>     store['foo'] = bar   # write to HDF5
    >>>     bar = store['foo']   # retrieve
    """
    store = None
    try:
        store = HDFStore(path, mode=mode, complevel=complevel,
                         complib=complib, fletcher32=False)
        yield store
    finally:
        if store is not None:
            store.close()


class HDFStore(object):
    """
    dict-like IO interface for storing pandas objects in PyTables
    format.

    DataFrame and Panel can be stored in Table format, which is slower to
    read and write but can be searched and manipulated more like an SQL
    table. See HDFStore.put for more information

    Parameters
    ----------
    path : string
        File path to HDF5 file
    mode : {'a', 'w', 'r', 'r+'}, default 'a'

        ``'r'``
            Read-only; no data can be modified.
        ``'w'``
            Write; a new file is created (an existing file with the same
            name would be deleted).
        ``'a'``
            Append; an existing file is opened for reading and writing,
            and if the file does not exist it is created.
        ``'r+'``
            It is similar to ``'a'``, but the file must already exist.
    complevel : int, 1-9, default 0
            If a complib is specified compression will be applied
            where possible
    complib : {'zlib', 'bzip2', 'lzo', 'blosc', None}, default None
            If complevel is > 0 apply compression to objects written
            in the store wherever possible
    fletcher32 : bool, default False
            If applying compression use the fletcher32 checksum

    Examples
    --------
    >>> store = HDFStore('test.h5')
    >>> store['foo'] = bar   # write to HDF5
    >>> bar = store['foo']   # retrieve
    >>> store.close()
    """
    _quiet = False

    def __init__(self, path, mode='a', complevel=None, complib=None,
                 fletcher32=False):
        try:
            import tables as _
        except ImportError:  # pragma: no cover
            raise Exception('HDFStore requires PyTables')

        self.path = path
        self.mode = mode
        self.handle = None
        self.complevel = complevel
        self.complib = complib
        self.fletcher32 = fletcher32
        self.filters = None
        self.open(mode=mode, warn=False)

    @property
    def root(self):
        """ return the root node """
        return self.handle.root

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.put(key, value)

    def __delitem__(self, key):
        return self.remove(key)

    def __contains__(self, key):
        """ check for existance of this key
              can match the exact pathname or the pathnm w/o the leading '/'
        """
        node = self.get_node(key)
        if node is not None:
            name = node._v_pathname
            return re.search(key,name) is not None
        return False

    def __len__(self):
        return len(self.groups())

    def __repr__(self):
        output = '%s\nFile path: %s\n' % (type(self), self.path)

        groups = self.groups()
        if len(groups) > 0:
            keys   = []
            values = []
            for n in sorted(groups, key = lambda x: x._v_name):
                kind = getattr(n._v_attrs,'pandas_type',None)

                keys.append(str(n._v_pathname))

                # a table
                if _is_table_type(n):
                    values.append(str(create_table(self, n)))
                
                # a group
                elif kind is None:
                    values.append('unknown type')

                # another type of pandas object
                else:
                    values.append(_NAME_MAP[kind])

            output += adjoin(5, keys, values)
        else:
            output += 'Empty'

        return output

    def keys(self):
        """
        Return a (potentially unordered) list of the keys corresponding to the
        objects stored in the HDFStore. These are ABSOLUTE path-names (e.g. have the leading '/'
        """
        return [ n._v_pathname for n in self.groups() ]

    def open(self, mode='a', warn=True):
        """
        Open the file in the specified mode

        Parameters
        ----------
        mode : {'a', 'w', 'r', 'r+'}, default 'a'
            See HDFStore docstring or tables.openFile for info about modes
        """
        self.mode = mode
        if warn and mode == 'w':  # pragma: no cover
            while True:
                response = raw_input("Re-opening as mode='w' will delete the "
                                     "current file. Continue (y/n)?")
                if response == 'y':
                    break
                elif response == 'n':
                    return
        if self.handle is not None and self.handle.isopen:
            self.handle.close()

        if self.complib is not None:
            if self.complevel is None:
                self.complevel = 9
            self.filters = _tables().Filters(self.complevel,
                                             self.complib,
                                             fletcher32=self.fletcher32)

        try:
            self.handle = _tables().openFile(self.path, self.mode)
        except IOError, e:  # pragma: no cover
            if 'can not be written' in str(e):
                print 'Opening %s in read-only mode' % self.path
                self.handle = _tables().openFile(self.path, 'r')
            else:
                raise

    def close(self):
        """
        Close the PyTables file handle
        """
        self.handle.close()

    def flush(self):
        """
        Force all buffered modifications to be written to disk
        """
        self.handle.flush()

    def get(self, key):
        """
        Retrieve pandas object stored in file

        Parameters
        ----------
        key : object

        Returns
        -------
        obj : type of object stored in file
        """
        group = self.get_node(key)
        if group is None:
            raise KeyError('No object named %s in the file' % key)
        return self._read_group(group)

    def select(self, key, where=None):
        """
        Retrieve pandas object stored in file, optionally based on where
        criteria

        Parameters
        ----------
        key : object
        where : list of Term (or convertable) objects, optional

        """
        group = self.get_node(key)
        if group is None:
            raise KeyError('No object named %s in the file' % key)
        if where is not None and not _is_table_type(group):
            raise Exception('can only select with where on objects written as tables')
        return self._read_group(group, where)

    def put(self, key, value, table=False, append=False,
            compression=None, **kwargs):
        """
        Store object in HDFStore

        Parameters
        ----------
        key : object
        value : {Series, DataFrame, Panel}
        table : boolean, default False
            Write as a PyTables Table structure which may perform worse but
            allow more flexible operations like searching / selecting subsets
            of the data
        append : boolean, default False
            For table data structures, append the input data to the existing
            table
        compression : {None, 'blosc', 'lzo', 'zlib'}, default None
            Use a compression algorithm to compress the data
            If None, the compression settings specified in the ctor will
            be used.
        """
        self._write_to_group(key, value, table=table, append=append,
                             comp=compression, **kwargs)

    def remove(self, key, where=None):
        """
        Remove pandas object partially by specifying the where condition

        Parameters
        ----------
        key : string
            Node to remove or delete rows from
        where : list
            For Table node, delete specified rows. See HDFStore.select for more
            information

        Returns
        -------
        number of rows removed (or None if not a Table)

        """
        group = self.get_node(key)
        if group is not None:

            # remove the node
            if where is None or not len(where):
                group = self.get_node(key)
                group._f_remove(recursive=True)
            
            # delete from the table
            else:
                if not _is_table_type(group):
                    raise Exception('can only remove with where on objects written as tables')
                t = create_table(self, group)
                return t.delete(where)

        return None

    def append(self, key, value, **kwargs):
        """
        Append to Table in file. Node must already exist and be Table
        format.

        Parameters
        ----------
        key : object
        value : {Series, DataFrame, Panel}

        Notes
        -----
        Does *not* check if data being appended overlaps with existing
        data in the table, so be careful
        """
        self._write_to_group(key, value, table=True, append=True, **kwargs)

    def create_table_index(self, key, **kwargs):
        """ Create a pytables index on the table
        Paramaters
        ----------
        key : object (the node to index)

        Exceptions
        ----------
        raises if the node is not a table

        """

        # version requirements
        if not _table_supports_index:
            raise("PyTables >= 2.3 is required for table indexing")

        group = self.get_node(key)
        if group is None: return

        if not _is_table_type(group):
            raise Exception("cannot create table index on a non-table")
        create_table(self, group).create_index(**kwargs)

    def groups(self):
        """ return a list of all the groups (that are not themselves a pandas storage object) """
        return [ g for g in self.handle.walkGroups() if getattr(g._v_attrs,'pandas_type',None) ]

    def get_node(self, key):
        """ return the node with the key or None if it does not exist """
        try:
            if not key.startswith('/'):
                key = '/' + key
            return self.handle.getNode(self.root,key)
        except:
            return None

    ###### private methods ######

    def _get_handler(self, op, kind):
        return getattr(self, '_%s_%s' % (op, kind))

    def _write_to_group(self, key, value, table=False, append=False,
                        comp=None, **kwargs):
        group = self.get_node(key)
        if group is None:
            paths = key.split('/')

            # recursively create the groups
            path = '/'
            for p in paths:
                if not len(p):
                    continue
                new_path = path
                if not path.endswith('/'):
                    new_path += '/'
                new_path += p
                group = self.get_node(new_path)
                if group is None:
                    group = self.handle.createGroup(path, p)
                path  = new_path

        kind = _TYPE_MAP[type(value)]
        if table or (append and _is_table_type(group)):
            kind = '%s_table' % kind
            handler = self._get_handler(op='write', kind=kind)
            wrapper = lambda value: handler(group, value, append=append,
                                            comp=comp, **kwargs)
        else:
            if append:
                raise ValueError('Can only append to Tables')
            if comp:
                raise ValueError('Compression only supported on Tables')

            handler = self._get_handler(op='write', kind=kind)
            wrapper = lambda value: handler(group, value)

        wrapper(value)
        group._v_attrs.pandas_type = kind

    def _write_series(self, group, series):
        self._write_index(group, 'index', series.index)
        self._write_array(group, 'values', series.values)
        group._v_attrs.name = series.name

    def _write_sparse_series(self, group, series):
        self._write_index(group, 'index', series.index)
        self._write_index(group, 'sp_index', series.sp_index)
        self._write_array(group, 'sp_values', series.sp_values)
        group._v_attrs.name = series.name
        group._v_attrs.fill_value = series.fill_value
        group._v_attrs.kind = series.kind

    def _read_sparse_series(self, group, where=None):
        index = self._read_index(group, 'index')
        sp_values = _read_array(group, 'sp_values')
        sp_index = self._read_index(group, 'sp_index')
        name = getattr(group._v_attrs, 'name', None)
        fill_value = getattr(group._v_attrs, 'fill_value', None)
        kind = getattr(group._v_attrs, 'kind', 'block')
        return SparseSeries(sp_values, index=index, sparse_index=sp_index,
                            kind=kind, fill_value=fill_value,
                            name=name)

    def _write_sparse_frame(self, group, sdf):
        for name, ss in sdf.iteritems():
            key = 'sparse_series_%s' % name
            if key not in group._v_children:
                node = self.handle.createGroup(group, key)
            else:
                node = getattr(group, key)
            self._write_sparse_series(node, ss)
        setattr(group._v_attrs, 'default_fill_value',
                sdf.default_fill_value)
        setattr(group._v_attrs, 'default_kind',
                sdf.default_kind)
        self._write_index(group, 'columns', sdf.columns)

    def _read_sparse_frame(self, group, where=None):
        columns = self._read_index(group, 'columns')
        sdict = {}
        for c in columns:
            key = 'sparse_series_%s' % c
            node = getattr(group, key)
            sdict[c] = self._read_sparse_series(node)
        default_kind = getattr(group._v_attrs, 'default_kind')
        default_fill_value = getattr(group._v_attrs, 'default_fill_value')
        return SparseDataFrame(sdict, columns=columns,
                               default_kind=default_kind,
                               default_fill_value=default_fill_value)

    def _write_sparse_panel(self, group, swide):
        setattr(group._v_attrs, 'default_fill_value', swide.default_fill_value)
        setattr(group._v_attrs, 'default_kind', swide.default_kind)
        self._write_index(group, 'items', swide.items)

        for name, sdf in swide.iteritems():
            key = 'sparse_frame_%s' % name
            if key not in group._v_children:
                node = self.handle.createGroup(group, key)
            else:
                node = getattr(group, key)
            self._write_sparse_frame(node, sdf)

    def _read_sparse_panel(self, group, where=None):
        default_fill_value = getattr(group._v_attrs, 'default_fill_value')
        default_kind = getattr(group._v_attrs, 'default_kind')
        items = self._read_index(group, 'items')

        sdict = {}
        for name in items:
            key = 'sparse_frame_%s' % name
            node = getattr(group, key)
            sdict[name] = self._read_sparse_frame(node)
        return SparsePanel(sdict, items=items, default_kind=default_kind,
                           default_fill_value=default_fill_value)

    def _write_frame(self, group, df):
        self._write_block_manager(group, df._data)

    def _read_frame(self, group, where=None):
        return DataFrame(self._read_block_manager(group))

    def _write_block_manager(self, group, data):
        if not data.is_consolidated():
            data = data.consolidate()

        group._v_attrs.ndim = data.ndim
        for i, ax in enumerate(data.axes):
            self._write_index(group, 'axis%d' % i, ax)

        # Supporting mixed-type DataFrame objects...nontrivial
        nblocks = len(data.blocks)
        group._v_attrs.nblocks = nblocks
        for i in range(nblocks):
            blk = data.blocks[i]
            # I have no idea why, but writing values before items fixed #2299
            self._write_array(group, 'block%d_values' % i, blk.values)
            self._write_index(group, 'block%d_items' % i, blk.items)

    def _read_block_manager(self, group):
        ndim = group._v_attrs.ndim

        axes = []
        for i in xrange(ndim):
            ax = self._read_index(group, 'axis%d' % i)
            axes.append(ax)

        items = axes[0]
        blocks = []
        for i in range(group._v_attrs.nblocks):
            blk_items = self._read_index(group, 'block%d_items' % i)
            values = _read_array(group, 'block%d_values' % i)
            blk = make_block(values, blk_items, items)
            blocks.append(blk)

        return BlockManager(blocks, axes)

    def _write_frame_table(self, group, df, append=False, comp=None, **kwargs):
        t = create_table(self, group, typ = 'appendable_frame')
        t.write(axes_to_index=[0], obj=df, append=append, compression=comp, **kwargs)

    def _write_wide(self, group, panel):
        panel._consolidate_inplace()
        self._write_block_manager(group, panel._data)

    def _read_wide(self, group, where=None):
        return Panel(self._read_block_manager(group))

    def _write_wide_table(self, group, panel, append=False, comp=None, **kwargs):
        t = create_table(self, group, typ = 'appendable_panel')
        t.write(axes_to_index=[1,2], obj=panel,
                append=append, compression=comp, **kwargs)

    def _read_wide_table(self, group, where=None):
        t = create_table(self, group)
        return t.read(where)

    def _write_index(self, group, key, index):
        if isinstance(index, MultiIndex):
            setattr(group._v_attrs, '%s_variety' % key, 'multi')
            self._write_multi_index(group, key, index)
        elif isinstance(index, BlockIndex):
            setattr(group._v_attrs, '%s_variety' % key, 'block')
            self._write_block_index(group, key, index)
        elif isinstance(index, IntIndex):
            setattr(group._v_attrs, '%s_variety' % key, 'sparseint')
            self._write_sparse_intindex(group, key, index)
        else:
            setattr(group._v_attrs, '%s_variety' % key, 'regular')
            converted = _convert_index(index).set_name('index')
            self._write_array(group, key, converted.values)
            node = getattr(group, key)
            node._v_attrs.kind = converted.kind
            node._v_attrs.name = index.name

            if isinstance(index, (DatetimeIndex, PeriodIndex)):
                node._v_attrs.index_class = _class_to_alias(type(index))

            if hasattr(index, 'freq'):
                node._v_attrs.freq = index.freq

            if hasattr(index, 'tz') and index.tz is not None:
                zone = tslib.get_timezone(index.tz)
                if zone is None:
                    zone = tslib.tot_seconds(index.tz.utcoffset())
                node._v_attrs.tz = zone

    def _read_index(self, group, key):
        variety = getattr(group._v_attrs, '%s_variety' % key)

        if variety == 'multi':
            return self._read_multi_index(group, key)
        elif variety == 'block':
            return self._read_block_index(group, key)
        elif variety == 'sparseint':
            return self._read_sparse_intindex(group, key)
        elif variety == 'regular':
            _, index = self._read_index_node(getattr(group, key))
            return index
        else:  # pragma: no cover
            raise Exception('unrecognized index variety: %s' % variety)

    def _write_block_index(self, group, key, index):
        self._write_array(group, '%s_blocs' % key, index.blocs)
        self._write_array(group, '%s_blengths' % key, index.blengths)
        setattr(group._v_attrs, '%s_length' % key, index.length)

    def _read_block_index(self, group, key):
        length = getattr(group._v_attrs, '%s_length' % key)
        blocs = _read_array(group, '%s_blocs' % key)
        blengths = _read_array(group, '%s_blengths' % key)
        return BlockIndex(length, blocs, blengths)

    def _write_sparse_intindex(self, group, key, index):
        self._write_array(group, '%s_indices' % key, index.indices)
        setattr(group._v_attrs, '%s_length' % key, index.length)

    def _read_sparse_intindex(self, group, key):
        length = getattr(group._v_attrs, '%s_length' % key)
        indices = _read_array(group, '%s_indices' % key)
        return IntIndex(length, indices)

    def _write_multi_index(self, group, key, index):
        setattr(group._v_attrs, '%s_nlevels' % key, index.nlevels)

        for i, (lev, lab, name) in enumerate(zip(index.levels,
                                                 index.labels,
                                                 index.names)):
            # write the level
            level_key = '%s_level%d' % (key, i)
            conv_level = _convert_index(lev).set_name(level_key)
            self._write_array(group, level_key, conv_level.values)
            node = getattr(group, level_key)
            node._v_attrs.kind = conv_level.kind
            node._v_attrs.name = name

            # write the name
            setattr(node._v_attrs, '%s_name%d' % (key, i), name)

            # write the labels
            label_key = '%s_label%d' % (key, i)
            self._write_array(group, label_key, lab)

    def _read_multi_index(self, group, key):
        nlevels = getattr(group._v_attrs, '%s_nlevels' % key)

        levels = []
        labels = []
        names = []
        for i in range(nlevels):
            level_key = '%s_level%d' % (key, i)
            name, lev = self._read_index_node(getattr(group, level_key))
            levels.append(lev)
            names.append(name)

            label_key = '%s_label%d' % (key, i)
            lab = _read_array(group, label_key)
            labels.append(lab)

        return MultiIndex(levels=levels, labels=labels, names=names)

    def _read_index_node(self, node):
        data = node[:]
        kind = node._v_attrs.kind
        name = None

        if 'name' in node._v_attrs:
            name = node._v_attrs.name

        index_class = _alias_to_class(getattr(node._v_attrs,
                                              'index_class', ''))
        factory = _get_index_factory(index_class)

        kwargs = {}
        if 'freq' in node._v_attrs:
            kwargs['freq'] = node._v_attrs['freq']

        if 'tz' in node._v_attrs:
            kwargs['tz'] = node._v_attrs['tz']

        if kind in ('date', 'datetime'):
            index = factory(_unconvert_index(data, kind), dtype=object,
                            **kwargs)
        else:
            index = factory(_unconvert_index(data, kind), **kwargs)

        index.name = name

        return name, index

    def _write_array(self, group, key, value):
        if key in group:
            self.handle.removeNode(group, key)

        # Transform needed to interface with pytables row/col notation
        empty_array = any(x == 0 for x in value.shape)
        transposed = False

        if not empty_array:
            value = value.T
            transposed = True

        if self.filters is not None:
            atom = None
            try:
                # get the atom for this datatype
                atom = _tables().Atom.from_dtype(value.dtype)
            except ValueError:
                pass

            if atom is not None:
                # create an empty chunked array and fill it from value
                ca = self.handle.createCArray(group, key, atom,
                                              value.shape,
                                              filters=self.filters)
                ca[:] = value
                getattr(group, key)._v_attrs.transposed = transposed
                return

        if value.dtype.type == np.object_:
            vlarr = self.handle.createVLArray(group, key,
                                              _tables().ObjectAtom())
            vlarr.append(value)
        elif value.dtype.type == np.datetime64:
            self.handle.createArray(group, key, value.view('i8'))
            getattr(group, key)._v_attrs.value_type = 'datetime64'
        else:
            if empty_array:
                # ugly hack for length 0 axes
                arr = np.empty((1,) * value.ndim)
                self.handle.createArray(group, key, arr)
                getattr(group, key)._v_attrs.value_type = str(value.dtype)
                getattr(group, key)._v_attrs.shape = value.shape
            else:
                self.handle.createArray(group, key, value)

        getattr(group, key)._v_attrs.transposed = transposed

    def _read_group(self, group, where=None):
        kind = group._v_attrs.pandas_type
        kind = _LEGACY_MAP.get(kind, kind)
        handler = self._get_handler(op='read', kind=kind)
        return handler(group, where)

    def _read_series(self, group, where=None):
        index = self._read_index(group, 'index')
        if len(index) > 0:
            values = _read_array(group, 'values')
        else:
            values = []

        name = getattr(group._v_attrs, 'name', None)
        return Series(values, index=index, name=name)

    def _read_legacy_series(self, group, where=None):
        index = self._read_index_legacy(group, 'index')
        values = _read_array(group, 'values')
        return Series(values, index=index)

    def _read_legacy_frame(self, group, where=None):
        index = self._read_index_legacy(group, 'index')
        columns = self._read_index_legacy(group, 'columns')
        values = _read_array(group, 'values')
        return DataFrame(values, index=index, columns=columns)

    def _read_index_legacy(self, group, key):
        node = getattr(group, key)
        data = node[:]
        kind = node._v_attrs.kind
        return _unconvert_index_legacy(data, kind)

    def _read_frame_table(self, group, where=None):
        t = create_table(self, group)
        return t.read(where)


class IndexCol(object):
    """ an index column description class

        Parameters
        ----------

        axis   : axis which I reference
        values : the ndarray like converted values
        kind   : a string description of this type
        typ    : the pytables type
        pos    : the position in the pytables

        """
    is_indexable = True

    def __init__(self, values = None, kind = None, typ = None, cname = None, itemsize = None, name = None, axis = None, kind_attr = None, pos = None, **kwargs):
        self.values = values
        self.kind   = kind
        self.typ    = typ
        self.itemsize = itemsize
        self.name      = name
        self.cname     = cname
        self.kind_attr = kind_attr
        self.axis      = axis
        self.pos       = pos
        self.table  = None

        if name is not None:
            self.set_name(name, kind_attr)
        if pos is not None:
            self.set_pos(pos)

    def set_name(self, name, kind_attr = None):
        """ set the name of this indexer """
        self.name      = name
        self.kind_attr = kind_attr or "%s_kind" % name
        if self.cname is None:
            self.cname = name

        return self

    def set_axis(self, axis):
        """ set the axis over which I index """
        self.axis      = axis

        return self

    def set_pos(self, pos):
        """ set the position of this column in the Table """
        self.pos       = pos
        if pos is not None and self.typ is not None:
            self.typ._v_pos = pos
        return self

    def set_table(self, table):
        self.table = table
        return self

    def __repr__(self):
        return "name->%s,cname->%s,axis->%s,pos->%s,kind->%s" % (self.name,self.cname,self.axis,self.pos,self.kind)

    __str__ = __repr__

    def copy(self):
        new_self = copy.copy(self)
        return new_self

    def infer(self, table):
        """ infer this column from the table: create and return a new object """
        new_self = self.copy()
        new_self.set_table(table)
        new_self.get_attr()
        return new_self

    def convert(self, sel):
        """ set the values from this selection """
        self.values = _maybe_convert(sel.values[self.cname], self.kind)

    @property
    def attrs(self):
        return self.table._v_attrs

    @property
    def description(self):
        return self.table.description

    @property
    def col(self):
        """ return my current col description """
        return getattr(self.description,self.cname,None)

    @property
    def cvalues(self):
        """ return my cython values """
        return self.values

    def __iter__(self):
        return iter(self.values)

    def maybe_set_size(self, min_itemsize = None, **kwargs):
        """ maybe set a string col itemsize:
               min_itemsize can be an interger or a dict with this columns name with an integer size """
        if self.kind == 'string':

            if isinstance(min_itemsize, dict):
                min_itemsize = min_itemsize.get(self.name)

            if min_itemsize is not None and self.typ.itemsize < min_itemsize:
                self.typ = _tables().StringCol(itemsize = min_itemsize, pos = self.pos)

    def validate_and_set(self, table, append, **kwargs):
        self.set_table(table)
        self.validate_col()
        self.validate_attr(append)
        self.set_attr()

    def validate_col(self):
        """ validate this column & set table data for it """

        # validate this column for string truncation (or reset to the max size)
        if self.kind == 'string':

            c = self.col
            if c is not None:
                if c.itemsize < self.itemsize:
                    raise Exception("[%s] column has a min_itemsize of [%s] but itemsize [%s] is required!" % (self.cname,self.itemsize,c.itemsize))


    def validate_attr(self, append):
        # check for backwards incompatibility
        if append:
            existing_kind = getattr(self.attrs,self.kind_attr,None)
            if existing_kind is not None and existing_kind != self.kind:
                raise TypeError("incompatible kind in col [%s - %s]" %
                                (existing_kind, self.kind))

    def get_attr(self):
        """ set the kind for this colummn """
        self.kind = getattr(self.attrs,self.kind_attr,None)

    def set_attr(self):
        """ set the kind for this colummn """
        setattr(self.attrs,self.kind_attr,self.kind)

class DataCol(IndexCol):
    """ a data holding column, by definition this is not indexable

        Parameters
        ----------

        data   : the actual data
        cname  : the column name in the table to hold the data (typeically values)
        """
    is_indexable = False

    @classmethod
    def create_for_block(cls, i, **kwargs):
        """ return a new datacol with the block i """
        return cls(name = 'values_%d' % i, cname = 'values_block_%d' % i, **kwargs)

    def __init__(self, values = None, kind = None, typ = None, cname = None, data = None, **kwargs):
        super(DataCol, self).__init__(values = values, kind = kind, typ = typ, cname = cname, **kwargs)
        self.dtype = None
        self.dtype_attr = "%s_dtype" % self.name
        self.set_data(data)

    def __repr__(self):
        return "name->%s,cname->%s,dtype->%s,shape->%s" % (self.name,self.cname,self.dtype,self.shape)

    def set_data(self, data):
        self.data = data
        if data is not None:
            if self.dtype is None:
                self.dtype = data.dtype.name

    @property
    def shape(self):
        return getattr(self.data,'shape',None)

    @property
    def cvalues(self):
        """ return my cython values """
        return self.data

    def validate_attr(self, append):
        """ validate that we have the same order as the existing & same dtype """
        if append:
            existing_fields = getattr(self.attrs, self.kind_attr, None)
            if (existing_fields is not None and
                existing_fields != list(self.values)):
                raise Exception("appended items do not match existing items"
                                " in table!")

            existing_dtype = getattr(self.attrs, self.dtype_attr, None)
            if (existing_dtype is not None and
                existing_dtype != self.dtype):
                raise Exception("appended items dtype do not match existing items dtype"
                                " in table!")

    def convert(self, sel):
        """ set the data from this selection (and convert to the correct dtype if we can) """
        self.set_data(sel.values[self.cname])

        # convert to the correct dtype
        if self.dtype is not None:
            try:
                self.data = self.data.astype(self.dtype)
            except:
                self.data = self.data.astype('O')

    def get_attr(self):
        """ get the data for this colummn """
        self.values = getattr(self.attrs,self.kind_attr,None)
        self.dtype  = getattr(self.attrs,self.dtype_attr,None)

    def set_attr(self):
        """ set the data for this colummn """
        setattr(self.attrs,self.kind_attr,self.values)
        if self.dtype is not None:
            setattr(self.attrs,self.dtype_attr,self.dtype)

class Table(object):
    """ represent a table:
          facilitate read/write of various types of tables
          this is an abstract base class

        Parameters
        ----------

        parent : my parent HDFStore
        group  : the group node where the table resides

        Attrs in Table Node
        -------------------
        These are attributes that are store in the main table node, they are necessary
        to recreate these tables when read back in.

        index_axes: a list of tuples of the (original indexing axis and index column)
        non_index_axes: a list of tuples of the (original index axis and columns on a non-indexing axis)
        values_axes : a list of the columns which comprise the data of this table

        """
    table_type = None
    ndim       = None

    def __init__(self, parent, group):
        self.parent      = parent
        self.group       = group
        self.index_axes     = []
        self.non_index_axes = []
        self.values_axes    = []
        self.selection      = None

    @property
    def pandas_type(self):
        return getattr(self.group._v_attrs,'pandas_type',None)

    def __repr__(self):
        """ return a pretty representatgion of myself """
        return "%s (typ->%s,nrows->%s)" % (self.pandas_type,self.table_type,self.nrows)

    __str__ = __repr__

    @property
    def nrows(self):
        return getattr(self.table,'nrows',None)

    @property
    def table(self):
        """ return the table group """
        return getattr(self.group, 'table', None)

    @property
    def handle(self):
        return self.parent.handle

    @property
    def _quiet(self):
        return self.parent._quiet

    @property
    def filters(self):
        return self.parent.filters

    @property
    def complevel(self):
        return self.parent.complevel

    @property
    def fletcher32(self):
        return self.parent.fletcher32

    @property
    def complib(self):
        return self.parent.complib

    @property
    def attrs(self):
        return self.group._v_attrs

    @property
    def description(self):
        return self.table.description

    @property
    def axes(self):
        return itertools.chain(self.index_axes, self.values_axes)

    def kinds_map(self):
        """ return a list of the kinds for each columns """
        return [ (a.cname,a.kind) for a in self.index_axes ]

    def index_cols(self):
        """ return a list of my index cols """
        return [ (i.axis,i.cname) for i in self.index_axes ]

    def values_cols(self):
        """ return a list of my values cols """
        return [ i.cname for i in self.values_axes ]

    def set_attrs(self):
        """ set our table type & indexables """
        self.attrs.table_type  = self.table_type
        self.attrs.index_cols  = self.index_cols()
        self.attrs.values_cols = self.values_cols()
        self.attrs.non_index_axes = self.non_index_axes

    def validate(self):
        """ raise if we have an incompitable table type with the current """
        et = getattr(self.attrs,'table_type',None)
        if et is not None and et != self.table_type:
                raise TypeError("incompatible table_type with existing [%s - %s]" %
                                (et, self.table_type))
        ic  = getattr(self.attrs,'index_cols',None)
        if ic is not None and ic != self.index_cols():
            raise TypeError("incompatible index cols with existing [%s - %s]" %
                            (ic, self.index_cols()))

    @property
    def indexables(self):
        """ create/cache the indexables if they don't exist """
        if self._indexables is None:

            d               = self.description
            self._indexables = []

            # index columns
            self._indexables.extend([ IndexCol(name = name, axis = axis, pos = i) for i, (axis, name) in enumerate(self.attrs.index_cols) ])

            # data columns
            base_pos = len(self._indexables)
            self._indexables.extend([ DataCol.create_for_block(i = i, pos = base_pos + i ) for i, c in enumerate(self.attrs.values_cols) ])

        return self._indexables

    def create_index(self, columns = None, optlevel = None, kind = None):
        """
        Create a pytables index on the specified columns
          note: cannot index Time64Col() currently; PyTables must be >= 2.3.1


        Paramaters
        ----------
        columns : None or list_like (the indexers to index)
        optlevel: optimization level (defaults to 6)
        kind    : kind of index (defaults to 'medium')

        Exceptions
        ----------
        raises if the node is not a table

        """

        table = self.table
        if table is None: return

        self.infer_axes()

        if columns is None:
            columns = [ self.index_axes[0].name ]
        if not isinstance(columns, (tuple,list)):
            columns = [ columns ]

        kw = dict()
        if optlevel is not None:
            kw['optlevel'] = optlevel
        if kind is not None:
            kw['kind']     = kind

        for c in columns:
            v = getattr(table.cols,c,None)
            if v is not None and not v.is_indexed:
                v.createIndex(**kw)

    def read_axes(self, where):
        """ create and return the axes sniffed from the table """

        # infer the data kind
        self.infer_axes()

        # create the selection
        self.selection = Selection(self, where)
        self.selection.select()

        # convert the data
        for a in self.axes:
            a.convert(self.selection)

    def infer_axes(self):
        """ infer the axes from the indexables """
        self.index_axes, self.values_axes = [ a.infer(self.table) for a in self.indexables if a.is_indexable ], [ a.infer(self.table) for a in self.indexables if not a.is_indexable ]
        self.non_index_axes = getattr(self.attrs,'non_index_axes',None) or []

    def create_axes(self, axes_to_index, obj, validate = True, min_itemsize = None):
        """ create and return the axes
              leagcy tables create an indexable column, indexable index, non-indexable fields

        """

        self.index_axes     = []
        self.non_index_axes = []

        # create axes to index and non_index
        j = 0
        for i, a in enumerate(obj.axes):

            if i in axes_to_index:
                name = obj._AXIS_NAMES[i]
                self.index_axes.append(_convert_index(a).set_name(name).set_axis(i).set_pos(j))
                j += 1

            else:
                self.non_index_axes.append((i,list(a)))

        # check for column conflicts
        if validate:
            for a in self.axes:
                a.maybe_set_size(min_itemsize = min_itemsize)

        # add my values
        self.values_axes = []
        for i, b in enumerate(obj._data.blocks):
            values = b.values

            # a string column
            if b.dtype.name == 'object':
                atom  = _tables().StringCol(itemsize = values.dtype.itemsize, shape = b.shape[0])
                utype = 'S8'
            else:
                atom  = getattr(_tables(),"%sCol" % b.dtype.name.capitalize())(shape = b.shape[0])
                utype = atom._deftype

            # coerce data to this type
            try:
                values = values.astype(utype)
            except (Exception), detail:
                raise Exception("cannot coerce data type -> [dtype->%s]" % b.dtype.name)

            dc = DataCol.create_for_block(i = i, values = list(b.items), kind = b.dtype.name, typ = atom, data = values, pos = j)
            j += 1
            self.values_axes.append(dc)

    def create_description(self, compression = None, complevel = None):
        """ create the description of the table from the axes & values """

        d = { 'name' : 'table' }

        # description from the axes & values
        d['description'] = dict([ (a.cname,a.typ) for a in self.axes ])

        if compression:
            complevel = self.complevel
            if complevel is None:
                complevel = 9
            filters = _tables().Filters(complevel=complevel,
                                        complib=compression,
                                        fletcher32=self.fletcher32)
            d['filters'] = filters
        elif self.filters is not None:
            d['filters'] = self.filters

        return d

    def read(self, **kwargs):
        raise NotImplementedError("cannot read on an abstract table: subclasses should implement")

    def write(self, **kwargs):
        raise NotImplementedError("cannot write on an abstract table")

    def delete(self, where = None, **kwargs):
        """ support fully deleting the node in its entirety (only) - where specification must be None """
        if where is None:
            self.handle.removeNode(self.group, recursive=True)
            return None

        raise NotImplementedError("cannot delete on an abstract table")

class WORMTable(Table):
    """ a write-once read-many table: this format DOES NOT ALLOW appending to a
         table. writing is a one-time operation the data are stored in a format
         that allows for searching the data on disk
         """
    table_type = 'worm'

    def read(self, **kwargs):
        """ read the indicies and the indexing array, calculate offset rows and
        return """
        raise NotImplementedError("WORMTable needs to implement read")

    def write(self, **kwargs):
        """ write in a format that we can search later on (but cannot append
               to): write out the indicies and the values using _write_array
               (e.g. a CArray) create an indexing table so that we can search"""
        raise NotImplementedError("WORKTable needs to implement write")

class LegacyTable(Table):
    """ an appendable table: allow append/query/delete operations to a
          (possibily) already existing appendable table this table ALLOWS
          append (but doesn't require them), and stores the data in a format
          that can be easily searched

        """
    _indexables = [IndexCol(name = 'index',  axis = 0, pos = 0),
                   IndexCol(name = 'column', axis = 1, pos = 1, index_kind = 'columns_kind'), 
                   DataCol( name = 'fields', cname = 'values', kind_attr = 'fields', pos = 2) ]
    table_type = 'legacy'

    def write(self, **kwargs):
        raise Exception("write operations are not allowed on legacy tables!")

    def read(self, where=None):
        """ we have 2 indexable columns, with an arbitrary number of data axes """

        self.read_axes(where)

        index  = self.index_axes[0].values
        column = self.index_axes[1].values

        major = Factor.from_array(index)
        minor = Factor.from_array(column)

        J, K = len(major.levels), len(minor.levels)
        key = major.labels * K + minor.labels

        panels = []
        if len(unique(key)) == len(key):
            sorter, _ = algos.groupsort_indexer(com._ensure_int64(key), J * K)
            sorter = com._ensure_platform_int(sorter)

            # create the panels
            for c in self.values_axes:

                # the data need to be sorted
                sorted_values = c.data.take(sorter, axis=0)
                major_labels = major.labels.take(sorter)
                minor_labels = minor.labels.take(sorter)
                items        = Index(c.values)

                block = block2d_to_block3d(sorted_values, items, (J, K),
                                           major_labels, minor_labels)

                mgr = BlockManager([block], [items, major.levels, minor.levels])
                panels.append(Panel(mgr))

        else:
            if not self._quiet:  # pragma: no cover
                print ('Duplicate entries in table, taking most recently '
                       'appended')

            # reconstruct
            long_index = MultiIndex.from_arrays([index, column])

            panels = []
            for c in self.values_axes:
                lp = DataFrame(c.data, index=long_index, columns=c.values)

                # need a better algorithm
                tuple_index = long_index._tuple_index

                unique_tuples = lib.fast_unique(tuple_index)
                unique_tuples = _asarray_tuplesafe(unique_tuples)

                indexer = match(unique_tuples, tuple_index)
                indexer = com._ensure_platform_int(indexer)

                new_index = long_index.take(indexer)
                new_values = lp.values.take(indexer, axis=0)

                lp = DataFrame(new_values, index=new_index, columns=lp.columns)
                panels.append(lp.to_panel())

        # append the panels
        wp = concat(panels, axis = 0, verify_integrity = True)

        # reorder by any non_index_axes
        for axis,labels in self.non_index_axes:
            wp = wp.reindex_axis(labels,axis=axis,copy=False)

        if self.selection.filter:
            new_minor = sorted(set(wp.minor_axis) & self.selection.filter)
            wp = wp.reindex(minor=new_minor, copy = False)

        return wp

class LegacyFrameTable(LegacyTable):
    """ support the legacy frame table """
    table_type = 'legacy_frame'
    def read(self, *args, **kwargs):
        return super(LegacyFrameTable, self).read(*args, **kwargs)['value']

class LegacyPanelTable(LegacyTable):
    """ support the legacy panel table """
    table_type = 'legacy_panel'

class AppendableTable(LegacyTable):
    """ suppor the new appendable table formats """
    _indexables = None
    table_type = 'appendable'

    def write(self, axes_to_index, obj, append=False, compression=None,
              complevel=None, min_itemsize = None, **kwargs):

        # create the table if it doesn't exist (or get it if it does)
        if not append:
            if 'table' in self.group:
                self.handle.removeNode(self.group, 'table')

        # create the axes
        self.create_axes(axes_to_index = axes_to_index, obj = obj, validate = append, min_itemsize = min_itemsize)

        if 'table' not in self.group:

            # create the table
            options = self.create_description(compression = compression, complevel = complevel)

            # set the table attributes
            self.set_attrs()

            # create the table
            table = self.handle.createTable(self.group, **options)

        else:

            # the table must already exist
            table = self.table

            # validate the table
            self.validate()

        # validate the axes and set the kinds
        for a in self.axes:
            a.validate_and_set(table, append)

        # add the rows
        self.write_data()
        self.handle.flush()

    def write_data(self):
        """ fast writing of data: requires specific cython routines each axis shape """

        masks  = []

        # create the masks
        for a in self.values_axes:

            # figure the mask: only do if we can successfully process this column, otherwise ignore the mask
            try:
                mask = np.isnan(a.data).all(axis=0)
                masks.append(mask.astype('u1'))
            except:

                # need to check for Nan in a non-numeric type column!!!
                masks.append(np.zeros((a.data.shape[1:]), dtype = 'u1'))

        # consolidate masks
        mask = masks[0]
        for m in masks[1:]:
            m = mask & m

        # the arguments & values
        args   = [ a.cvalues for a in self.index_axes ]
        values = [ a.data for a in self.values_axes ]

        # get our function
        try:
            func = getattr(lib,"create_hdf_rows_%sd" % self.ndim)
            args.append(mask)
            args.append(values)
            rows = func(*args)
            if len(rows):
                self.table.append(rows)
        except (Exception), detail:
            raise Exception("tables cannot write this data -> %s" % str(detail))

    def delete(self, where = None):
        if where is None:
            return super(LegacyTable, self).delete()

        # infer the data kind
        table = self.table
        self.infer_axes()

        # create the selection
        self.selection = Selection(self, where)
        self.selection.select_coords()

        # delete the rows in reverse order
        l  = list(self.selection.values)
        ln = len(l)

        if ln:

            # if we can do a consecutive removal - do it!
            if l[0]+ln-1 == l[-1]:
                table.removeRows(start = l[0], stop = l[-1]+1)

            # one by one
            else:
                l.reverse()
                for c in l:
                    table.removeRows(c)

                    self.handle.flush()

        # return the number of rows removed
        return ln

class AppendableFrameTable(AppendableTable):
    """ suppor the new appendable table formats """
    table_type = 'appendable_frame'
    ndim       = 2

    def read(self, where=None):

        self.read_axes(where)

        index   = Index(self.index_axes[0].values)
        frames  = []
        for a in self.values_axes:
            columns = Index(a.values)
            block   = make_block(a.cvalues.T, columns, columns)
            mgr     = BlockManager([ block ], [ columns, index ])
            frames.append(DataFrame(mgr))
        df = concat(frames, axis = 1, verify_integrity = True)

        # sort the indicies & reorder the columns
        for axis,labels in self.non_index_axes:
            df = df.reindex_axis(labels,axis=axis,copy=False)
        columns_ordered = df.columns

        # apply the column filter (but keep columns in the same order)
        if self.selection.filter:
            columns = Index(set(columns_ordered) & self.selection.filter)
            columns = sorted(columns_ordered.get_indexer(columns))
            df      = df.reindex(columns = columns_ordered.take(columns), copy = False)

        else:
            df      = df.reindex(columns = columns_ordered, copy = False)

        return df

class AppendablePanelTable(AppendableTable):
    """ suppor the new appendable table formats """
    table_type = 'appendable_panel'
    ndim       = 3

# table maps
_TABLE_MAP = {
    'appendable_frame' : AppendableFrameTable,
    'appendable_panel' : AppendablePanelTable,
    'worm'             : WORMTable,
    'legacy_frame'     : LegacyFrameTable,
    'legacy_panel'     : LegacyPanelTable,
    'default'          : AppendablePanelTable,
}

def create_table(parent, group, typ = None, **kwargs):
    """ return a suitable Table class to operate """

    pt = getattr(group._v_attrs,'pandas_type',None)
    tt = getattr(group._v_attrs,'table_type',None)

    # a new node
    if pt is None:

        return (_TABLE_MAP.get(typ) or _TABLE_MAP.get('default'))(parent, group, **kwargs)

    # existing node (legacy)
    if tt is None:

        # distiguish between a frame/table
        tt = 'legacy_panel'
        try:
            if group.table.description.values.shape[0] == 1:
                tt = 'legacy_frame'
        except:
            pass

    return _TABLE_MAP.get(tt)(parent, group, **kwargs)


def _convert_index(index):
    if isinstance(index, DatetimeIndex):
        converted = index.asi8
        return IndexCol(converted, 'datetime64', _tables().Int64Col())
    elif isinstance(index, (Int64Index, PeriodIndex)):
        atom = _tables().Int64Col()
        return IndexCol(index.values, 'integer', atom)

    if isinstance(index, MultiIndex):
        raise Exception('MultiIndex not supported here!')

    inferred_type = lib.infer_dtype(index)

    values = np.asarray(index)

    if inferred_type == 'datetime64':
        converted = values.view('i8')
        return IndexCol(converted, 'datetime64', _tables().Int64Col())
    elif inferred_type == 'datetime':
        converted = np.array([(time.mktime(v.timetuple()) +
                            v.microsecond / 1E6) for v in values],
                            dtype=np.float64)
        return IndexCol(converted, 'datetime', _tables().Time64Col())
    elif inferred_type == 'date':
        converted = np.array([time.mktime(v.timetuple()) for v in values],
                            dtype=np.int32)
        return IndexCol(converted, 'date', _tables().Time32Col())
    elif inferred_type == 'string':
        # atom = _tables().ObjectAtom()
        # return np.asarray(values, dtype='O'), 'object', atom

        converted = np.array(list(values), dtype=np.str_)
        itemsize = converted.dtype.itemsize
        return IndexCol(converted, 'string', _tables().StringCol(itemsize), itemsize = itemsize)
    elif inferred_type == 'unicode':
        atom = _tables().ObjectAtom()
        return IndexCol(np.asarray(values, dtype='O'), 'object', atom)
    elif inferred_type == 'integer':
        # take a guess for now, hope the values fit
        atom = _tables().Int64Col()
        return IndexCol(np.asarray(values, dtype=np.int64), 'integer', atom)
    elif inferred_type == 'floating':
        atom = _tables().Float64Col()
        return IndexCol(np.asarray(values, dtype=np.float64), 'float', atom)
    else:  # pragma: no cover
        atom = _tables().ObjectAtom()
        return IndexCol(np.asarray(values, dtype='O'), 'object', atom)


def _read_array(group, key):
    import tables
    node = getattr(group, key)
    data = node[:]
    attrs = node._v_attrs

    transposed = getattr(attrs, 'transposed', False)

    if isinstance(node, tables.VLArray):
        ret = data[0]
    else:
        dtype = getattr(attrs, 'value_type', None)
        shape = getattr(attrs, 'shape', None)

        if shape is not None:
            # length 0 axis
            ret = np.empty(shape, dtype=dtype)
        else:
            ret = data

        if dtype == 'datetime64':
            ret = np.array(ret, dtype='M8[ns]')

    if transposed:
        return ret.T
    else:
        return ret


def _unconvert_index(data, kind):
    if kind == 'datetime64':
        index = DatetimeIndex(data)
    elif kind == 'datetime':
        index = np.array([datetime.fromtimestamp(v) for v in data],
                         dtype=object)
    elif kind == 'date':
        index = np.array([date.fromtimestamp(v) for v in data], dtype=object)
    elif kind in ('string', 'integer', 'float'):
        index = np.array(data)
    elif kind == 'object':
        index = np.array(data[0])
    else:  # pragma: no cover
        raise ValueError('unrecognized index type %s' % kind)
    return index


def _unconvert_index_legacy(data, kind, legacy=False):
    if kind == 'datetime':
        index = lib.time64_to_datetime(data)
    elif kind in ('string', 'integer'):
        index = np.array(data, dtype=object)
    else:  # pragma: no cover
        raise ValueError('unrecognized index type %s' % kind)
    return index


def _maybe_convert(values, val_kind):
    if _need_convert(val_kind):
        conv = _get_converter(val_kind)
        # conv = np.frompyfunc(conv, 1, 1)
        values = conv(values)
    return values


def _get_converter(kind):
    if kind == 'datetime64':
        return lambda x: np.array(x, dtype='M8[ns]')
    if kind == 'datetime':
        return lib.convert_timestamps
    else:  # pragma: no cover
        raise ValueError('invalid kind %s' % kind)


def _need_convert(kind):
    if kind in ('datetime', 'datetime64'):
        return True
    return False


def _is_table_type(group):
    try:
        return 'table' in group._v_attrs.pandas_type
    except AttributeError:
        # new node, e.g.
        return False

_index_type_map = {DatetimeIndex: 'datetime',
                   PeriodIndex: 'period'}

_reverse_index_map = {}
for k, v in _index_type_map.iteritems():
    _reverse_index_map[v] = k


def _class_to_alias(cls):
    return _index_type_map.get(cls, '')


def _alias_to_class(alias):
    if isinstance(alias, type):  # pragma: no cover
        return alias  # compat: for a short period of time master stored types
    return _reverse_index_map.get(alias, Index)


class Term(object):
    """ create a term object that holds a field, op, and value

        Parameters
        ----------
        field : dict, string term expression, or the field to operate (must be a valid index/column type of DataFrame/Panel)
        op    : a valid op (defaults to '=') (optional)
                >, >=, <, <=, =, != (not equal) are allowed
        value : a value or list of values (required)
        kinds : the kinds map (dict of column name -> kind)

        Returns
        -------
        a Term object

        Examples
        --------
        Term(dict(field = 'index', op = '>', value = '20121114'))
        Term('index', '20121114')
        Term('index', '>', '20121114')
        Term('index', ['20121114','20121114'])
        Term('index', datetime(2012,11,14))
        Term('major>20121114')
        Term('minor', ['A','B'])

    """

    _ops     = ['<=','<','>=','>','!=','=']
    _search  = re.compile("^(?P<field>\w+)(?P<op>%s)(?P<value>.+)$" % '|'.join(_ops))
    _index   = ['index','major_axis','major']
    _column  = ['column','columns','minor_axis','minor']

    def __init__(self, field, op = None, value = None, kinds = None):
        self.field = None
        self.op    = None
        self.value = None

        if kinds is None:
            kinds = []
        self.kinds = dict(kinds)
        self.filter     = None
        self.condition  = None

        # unpack lists/tuples in field
        while(isinstance(field,(tuple,list))):
            f = field
            field = f[0]
            if len(f) > 1:
                op = f[1]
            if len(f) > 2:
                value = f[2]

        # backwards compatible
        if isinstance(field, dict):
            self.field = field.get('field')
            self.op    = field.get('op') or '='
            self.value = field.get('value')

        # passed a term
        elif isinstance(field,Term):
            self.field = field.field
            self.op    = field.op
            self.value = field.value

        # a string expression (or just the field)
        elif isinstance(field,basestring):

            # is a term is passed
            s = self._search.match(field)
            if s is not None:
                self.field = s.group('field')
                self.op    = s.group('op')
                self.value = s.group('value')

            else:
                self.field = field

                # is an op passed?
                if isinstance(op, basestring) and op in self._ops:
                    self.op    = op
                    self.value = value
                else:
                    self.op    = '='
                    self.value = op

        else:
            raise Exception("Term does not understand the supplied field [%s]" % field)

        # we have valid fields
        if self.field is None or self.op is None or self.value is None:
            raise Exception("Could not create this term [%s]" % str(self))

        # map alias for field names
        if self.field in self._index and len(kinds) > 0:
            self.field = kinds[0][0]
        elif self.field in self._column and len(kinds) > 1:
            self.field = kinds[1][0]

        # we have valid conditions
        if self.op in ['>','>=','<','<=']:
            if hasattr(self.value,'__iter__') and len(self.value) > 1:
                raise Exception("an inequality condition cannot have multiple values [%s]" % str(self))

        if not hasattr(self.value,'__iter__'):
            self.value = [ self.value ]

        self.eval()

    def __str__(self):
        return "field->%s,op->%s,value->%s" % (self.field,self.op,self.value)

    __repr__ = __str__

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an actual column in the table) """
        return self.field in self.kinds

    @property
    def kind(self):
        """ the kind of my field """
        return self.kinds.get(self.field)

    def eval(self):
        """ set the numexpr expression for this term """

        # convert values
        values = [ self.convert_value(v) for v in self.value ]

        # equality conditions
        if self.op in ['=','!=']:

            if self.is_in_table:

                # too many values to create the expression?
                if len(values) <= 61:
                    self.condition = "(%s)" % ' | '.join([ "(%s == %s)" % (self.field,v[0]) for v in values])

                # use a filter after reading
                else:
                    self.filter = set([ v[1] for v in values ])

            else:

                self.filter = set([ v[1] for v in values ])

        else:

            if self.is_in_table:

                self.condition = '(%s %s %s)' % (self.field, self.op, values[0][0])

    def convert_value(self, v):

        #### a little hacky here, need to really figure out what we should convert ####x
        if self.field == 'index' or self.field == 'major_axis':
            if self.kind == 'datetime64' :
                return [lib.Timestamp(v).value, None]
            elif isinstance(v, datetime):
                return [time.mktime(v.timetuple()), None]
        elif not isinstance(v, basestring):
            return [str(v), None]

        # string quoting
        return ["'" + v + "'", v]

class Selection(object):
    """
    Carries out a selection operation on a tables.Table object.

    Parameters
    ----------
    table : a Table object
    where : list of Terms (or convertable to)

    """
    def __init__(self, table, where=None):
        self.table      = table
        self.where      = where
        self.values     = None
        self.condition  = None
        self.filter     = None
        self.terms      = self.generate(where)

        # create the numexpr & the filter
        if self.terms:
            conds = [ t.condition for t in self.terms if t.condition is not None ]
            if len(conds):
                self.condition = "(%s)" % ' & '.join(conds)
            self.filter    = set()
            for t in self.terms:
                if t.filter is not None:
                    self.filter |= t.filter

    def generate(self, where):
        """ where can be a : dict,list,tuple,string """
        if where is None: return None

        if not isinstance(where, (list,tuple)):
            where = [ where ]
        else:
            # do we have all list/tuple
            if not any([ isinstance(w, (list,tuple,Term)) for w in where ]):
                where = [ where ]

        return [ Term(c, kinds = self.table.kinds_map()) for c in where ]

    def select(self):
        """
        generate the selection
        """
        if self.condition is not None:
            self.values = self.table.table.readWhere(self.condition)
        else:
            self.values = self.table.table.read()

    def select_coords(self):
        """
        generate the selection
        """
        self.values = self.table.table.getWhereList(self.condition)


def _get_index_factory(klass):
    if klass == DatetimeIndex:
        def f(values, freq=None, tz=None):
            return DatetimeIndex._simple_new(values, None, freq=freq,
                                             tz=tz)
        return f
    return klass
