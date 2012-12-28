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
import warnings

import numpy as np
from pandas import (
    Series, TimeSeries, DataFrame, Panel, Panel4D, Index, MultiIndex, Int64Index
)
from pandas.sparse.api import SparseSeries, SparseDataFrame, SparsePanel
from pandas.sparse.array import BlockIndex, IntIndex
from pandas.tseries.api import PeriodIndex, DatetimeIndex
from pandas.core.common import adjoin
from pandas.core.algorithms import match, unique, factorize
from pandas.core.strings import str_len
from pandas.core.categorical import Categorical
from pandas.core.common import _asarray_tuplesafe, _try_sort
from pandas.core.internals import BlockManager, make_block, form_blocks
from pandas.core.reshape import block2d_to_block3d, block2d_to_blocknd, factor_indexer
import pandas.core.common as com
from pandas.tools.merge import concat

import pandas.lib as lib
import pandas.algos as algos
import pandas.tslib as tslib

from contextlib import contextmanager

# versioning attribute
_version = '0.10'

class IncompatibilityWarning(Warning): pass

# reading and writing the full object in one go
_TYPE_MAP = {
    Series: 'series',
    SparseSeries: 'sparse_series',
    TimeSeries: 'series',
    DataFrame: 'frame',
    SparseDataFrame: 'sparse_frame',
    Panel: 'wide',
    Panel4D : 'ndim',
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
        ver = tables.__version__.split('.')
        try:
            if int(ver[0]) >= 2 and int(ver[1][0]) >= 3:
                _table_supports_index = True
        except:
            pass

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

    def select(self, key, where=None, **kwargs):
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
        return self._read_group(group, where, **kwargs)

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
            if where is None:
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
        _tables()
        if not _table_supports_index:
            raise Exception("PyTables >= 2.3 is required for table indexing")

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
        group._v_attrs.pandas_version = _version
        #group._v_attrs.meta = getattr(value,'meta',None)

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

    def _write_wide(self, group, panel):
        panel._consolidate_inplace()
        self._write_block_manager(group, panel._data)

    def _read_wide(self, group, where=None):
        return Panel(self._read_block_manager(group))

    def _write_ndim_table(self, group, obj, append=False, comp=None, axes=None, **kwargs):
        if axes is None:
            axes = [1,2,3]
        t = create_table(self, group, typ = 'appendable_ndim')
        t.write(axes=axes, obj=obj,
                append=append, compression=comp, **kwargs)

    def _read_ndim_table(self, group, where=None, **kwargs):
        t = create_table(self, group, **kwargs)
        return t.read(where)

    def _write_frame_table(self, group, df, append=False, comp=None, axes=None, **kwargs):
        if axes is None:
            axes = [0]
        t = create_table(self, group, typ = 'appendable_frame')
        t.write(axes=axes, obj=df, append=append, compression=comp, **kwargs)

    _read_frame_table = _read_ndim_table

    def _write_wide_table(self, group, panel, append=False, comp=None, axes=None, **kwargs):
        if axes is None:
            axes = [1,2]
        t = create_table(self, group, typ = 'appendable_panel')
        t.write(axes=axes, obj=panel,
                append=append, compression=comp, **kwargs)

    _read_wide_table = _read_ndim_table

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

    def _read_group(self, group, where=None, **kwargs):
        kind = group._v_attrs.pandas_type
        kind = _LEGACY_MAP.get(kind, kind)
        handler = self._get_handler(op='read', kind=kind)
        v = handler(group, where, **kwargs)
        #if v is not None:
        #    meta = getattr(group._v_attrs,'meta',None)
        #    if meta is not None:
        #        v.meta = meta
        return v

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

    def __eq__(self, other):
        """ compare 2 col items """
        return all([ getattr(self,a,None) == getattr(other,a,None) for a in ['name','cname','axis','pos'] ])

    def copy(self):
        new_self = copy.copy(self)
        return new_self

    def infer(self, table):
        """ infer this column from the table: create and return a new object """
        new_self = self.copy()
        new_self.set_table(table)
        new_self.get_attr()
        return new_self

    def convert(self, values):
        """ set the values from this selection """
        self.values = Index(_maybe_convert(values[self.cname], self.kind))

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

    def validate_col(self, itemsize = None):
        """ validate this column: return the compared against itemsize """

        # validate this column for string truncation (or reset to the max size)
        dtype = getattr(self,'dtype',None)
        if self.kind == 'string' or (dtype is not None and dtype.startswith('string')):

            c = self.col
            if c is not None:
                if itemsize is None:
                    itemsize = self.itemsize
                if c.itemsize < itemsize:
                    raise Exception("[%s] column has a min_itemsize of [%s] but itemsize [%s] is required!" % (self.cname,itemsize,c.itemsize))
                return c.itemsize

        return None


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

    def __eq__(self, other):
        """ compare 2 col items """
        return all([ getattr(self,a,None) == getattr(other,a,None) for a in ['name','cname','dtype','pos'] ])

    def set_data(self, data):
        self.data = data
        if data is not None:
            if self.dtype is None:
                self.dtype = data.dtype.name

    def take_data(self):
        """ return the data & release the memory """
        self.data, data = None, self.data
        return data

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

    def convert(self, values):
        """ set the data from this selection (and convert to the correct dtype if we can) """
        self.set_data(values[self.cname])

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
    obj_type   = None
    ndim       = None

    def __init__(self, parent, group, **kwargs):
        self.parent      = parent
        self.group       = group
        self.version     = getattr(group._v_attrs,'pandas_version',None)
        self.index_axes     = []
        self.non_index_axes = []
        self.values_axes    = []
        self.selection      = None

    @property
    def table_type_short(self):
        return self.table_type.split('_')[0]

    @property
    def pandas_type(self):
        return getattr(self.group._v_attrs,'pandas_type',None)

    def __repr__(self):
        """ return a pretty representatgion of myself """
        self.infer_axes()
        return "%s (typ->%s,nrows->%s,indexers->[%s])" % (self.pandas_type,self.table_type_short,self.nrows,','.join([ a.name for a in self.index_axes ]))

    __str__ = __repr__

    def copy(self):
        new_self = copy.copy(self)
        return new_self

    def validate(self, other):
        """ validate against an existing table """
        if other is None: return

        if other.table_type != self.table_type:
            raise TypeError("incompatible table_type with existing [%s - %s]" %
                            (other.table_type, self.table_type))

        for c in ['index_axes','non_index_axes','values_axes']:
            if getattr(self,c,None) != getattr(other,c,None):
                raise Exception("invalid combinate of [%s] on appending data [%s] vs current table [%s]" % (c,getattr(self,c,None),getattr(other,c,None)))

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

    @property
    def is_transposed(self):
        return False

    @property
    def data_orientation(self):
        """ return a tuple of my permutated axes, non_indexable at the front """
        return tuple(itertools.chain([ a[0] for a in self.non_index_axes ], [ a.axis for a in self.index_axes ]))

    def queryables(self):
        """ return a dict of the kinds allowable columns for this object """
        return dict([ (a.cname,a.kind) for a in self.index_axes ] + [ (self.obj_type._AXIS_NAMES[axis],None) for axis, values in self.non_index_axes ])

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

    def validate_version(self, where = None):
        """ are we trying to operate on an old version? """
        if where is not None:
            if self.version is None or float(self.version) < 0.1:
                warnings.warn("where criteria is being ignored as we this version is too old (or not-defined) [%s]" % self.version, IncompatibilityWarning)

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

        if not self.infer_axes(): return

        if columns is None:
            columns = [ self.index_axes[0].name ]
        if not isinstance(columns, (tuple,list)):
            columns = [ columns ]

        kw = dict()
        if optlevel is not None:
            kw['optlevel'] = optlevel
        if kind is not None:
            kw['kind']     = kind

        table = self.table
        for c in columns:
            v = getattr(table.cols,c,None)
            if v is not None:

                # remove the index if the kind/optlevel have changed
                if v.is_indexed:
                    index = v.index
                    cur_optlevel = index.optlevel
                    cur_kind     = index.kind

                    if kind is not None and cur_kind != kind:
                        v.removeIndex()
                    else:
                        kw['kind'] = cur_kind

                    if optlevel is not None and cur_optlevel != optlevel:
                        v.removeIndex()
                    else:
                        kw['optlevel'] = cur_optlevel

                # create the index
                if not v.is_indexed:
                    v.createIndex(**kw)

    def read_axes(self, where):
        """ create and return the axes sniffed from the table: return boolean for success """

        # validate the version
        self.validate_version(where)

        # infer the data kind
        if not self.infer_axes(): return False

        # create the selection
        self.selection = Selection(self, where)
        values = self.selection.select()

        # convert the data
        for a in self.axes:
            a.convert(values)

        return True

    def infer_axes(self):
        """ infer the axes from the indexables:
              return a boolean indicating if we have a valid table or not """

        table = self.table
        if table is None:
            return False

        self.index_axes, self.values_axes = [ a.infer(self.table) for a in self.indexables if a.is_indexable ], [ a.infer(self.table) for a in self.indexables if not a.is_indexable ]
        self.non_index_axes   = getattr(self.attrs,'non_index_axes',None) or []

        return True

    def get_data_blocks(self, obj):
        """ return the data blocks for this obj """
        return obj._data.blocks

    def create_axes(self, axes, obj, validate = True, min_itemsize = None):
        """ create and return the axes
              leagcy tables create an indexable column, indexable index, non-indexable fields

        """

        # map axes to numbers
        axes = [ obj._get_axis_number(a) for a in axes ]

        # do we have an existing table (if so, use its axes)?
        if self.infer_axes():
            existing_table = self.copy()
            axes = [ a.axis for a in existing_table.index_axes]
        else:
            existing_table = None

        # currently support on ndim-1 axes
        if len(axes) != self.ndim-1:
            raise Exception("currenctly only support ndim-1 indexers in an AppendableTable")

        # create according to the new data
        self.non_index_axes   = []

        # create axes to index and non_index
        index_axes_map = dict()
        for i, a in enumerate(obj.axes):

            if i in axes:
                name = obj._AXIS_NAMES[i]
                index_axes_map[i] = _convert_index(a).set_name(name).set_axis(i)
            else:

                # we might be able to change the axes on the appending data if necessary
                append_axis = list(a)
                if existing_table is not None:
                    indexer = len(self.non_index_axes)
                    exist_axis = existing_table.non_index_axes[indexer][1]
                    if append_axis != exist_axis:

                        # ahah! -> reindex
                        if sorted(append_axis) == sorted(exist_axis):
                            append_axis = exist_axis

                self.non_index_axes.append((i,append_axis))

        # set axis positions (based on the axes)
        self.index_axes = [ index_axes_map[a].set_pos(j) for j, a in enumerate(axes) ]
        j = len(self.index_axes)

        # check for column conflicts
        if validate:
            for a in self.axes:
                a.maybe_set_size(min_itemsize = min_itemsize)

        # reindex by our non_index_axes
        for a in self.non_index_axes:
            obj = obj.reindex_axis(a[1], axis = a[0], copy = False)

        blocks = self.get_data_blocks(obj)

        # add my values
        self.values_axes = []
        for i, b in enumerate(blocks):

            # shape of the data column are the indexable axes
            shape  = b.shape[0]
            values = b.values

            # a string column
            if b.dtype.name == 'object':

                # itemsize is the maximum length of a string (along any dimension)
                itemsize = _itemsize_string_array(values)

                # specified min_itemsize?
                if isinstance(min_itemsize, dict):
                    itemsize = max(int(min_itemsize.get('values')),itemsize)

                # check for column in the values conflicts
                if existing_table is not None and validate:
                    eci = existing_table.values_axes[i].validate_col(itemsize)
                    if eci > itemsize:
                        itemsize = eci

                atom  = _tables().StringCol(itemsize = itemsize, shape = shape)
                utype = 'S%s' % itemsize
                kind  = 'string'

            else:
                atom  = getattr(_tables(),"%sCol" % b.dtype.name.capitalize())(shape = shape)
                utype = atom._deftype
                kind  = b.dtype.name

            # coerce data to this type
            try:
                values = values.astype(utype)
            except (Exception), detail:
                raise Exception("cannot coerce data type -> [dtype->%s]" % b.dtype.name)

            dc = DataCol.create_for_block(i = i, values = list(b.items), kind = kind, typ = atom, data = values, pos = j)
            j += 1
            self.values_axes.append(dc)

        # validate the axes if we have an existing table
        if validate:
            self.validate(existing_table)

    def process_axes(self, obj):
        """ process axes filters """

        def reindex(obj, axis, filt, ordered):
            axis_name = obj._get_axis_name(axis)
            ordd = ordered & filt
            ordd = sorted(ordered.get_indexer(ordd))
            return obj.reindex_axis(ordered.take(ordd), axis = obj._get_axis_number(axis_name), copy = False)

        # apply the selection filters (but keep in the same order)
        if self.selection.filter:
            for axis, filt in self.selection.filter:
                obj = reindex(obj, axis, filt, getattr(obj,obj._get_axis_name(axis)))

        return obj

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
    _indexables = [IndexCol(name = 'index',  axis = 1, pos = 0),
                   IndexCol(name = 'column', axis = 2, pos = 1, index_kind = 'columns_kind'),
                   DataCol( name = 'fields', cname = 'values', kind_attr = 'fields', pos = 2) ]
    table_type = 'legacy'
    ndim       = 3

    def write(self, **kwargs):
        raise Exception("write operations are not allowed on legacy tables!")

    def read(self, where=None):
        """ we have n indexable columns, with an arbitrary number of data axes """


        if not self.read_axes(where): return None

        factors  = [ Categorical.from_array(a.values) for a in self.index_axes ]
        levels   = [ f.levels for f in factors ]
        N        = [ len(f.levels) for f in factors ]
        labels   = [ f.labels for f in factors ]

        # compute the key
        key      = factor_indexer(N[1:], labels)

        objs = []
        if len(unique(key)) == len(key):

            sorter, _ = algos.groupsort_indexer(com._ensure_int64(key), np.prod(N))
            sorter = com._ensure_platform_int(sorter)

            # create the objs
            for c in self.values_axes:

                # the data need to be sorted
                sorted_values = c.take_data().take(sorter, axis=0)

                take_labels   = [ l.take(sorter) for l in labels ]
                items         = Index(c.values)
                block         = block2d_to_blocknd(sorted_values, items, tuple(N), take_labels)

                # create the object
                mgr = BlockManager([block], [items] + levels)
                obj = self.obj_type(mgr)

                # permute if needed
                if self.is_transposed:
                    obj = obj.transpose(*tuple(Series(self.data_orientation).argsort()))

                objs.append(obj)

        else:
            if not self._quiet:  # pragma: no cover
                print ('Duplicate entries in table, taking most recently '
                       'appended')

            # reconstruct
            long_index = MultiIndex.from_arrays([ i.values for i in self.index_axes ])


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
                objs.append(lp.to_panel())

        # create the composite object
        if len(objs) == 1:
            wp = objs[0]
        else:
            wp = concat(objs, axis = 0, verify_integrity = True)

        # reorder by any non_index_axes
        for axis,labels in self.non_index_axes:
            wp = wp.reindex_axis(labels,axis=axis,copy=False)

        # apply the selection filters & axis orderings
        wp = self.process_axes(wp)

        return wp

class LegacyFrameTable(LegacyTable):
    """ support the legacy frame table """
    table_type = 'legacy_frame'
    obj_type   = Panel
    def read(self, *args, **kwargs):
        return super(LegacyFrameTable, self).read(*args, **kwargs)['value']

class LegacyPanelTable(LegacyTable):
    """ support the legacy panel table """
    table_type = 'legacy_panel'
    obj_type   = Panel

class AppendableTable(LegacyTable):
    """ suppor the new appendable table formats """
    _indexables = None
    table_type = 'appendable'

    def write(self, axes, obj, append=False, compression=None,
              complevel=None, min_itemsize = None, **kwargs):

        # create the table if it doesn't exist (or get it if it does)
        if not append:
            if 'table' in self.group:
                self.handle.removeNode(self.group, 'table')

        # create the axes
        self.create_axes(axes = axes, obj = obj, validate = append, min_itemsize = min_itemsize)

        if 'table' not in self.group:

            # create the table
            options = self.create_description(compression = compression, complevel = complevel)

            # set the table attributes
            self.set_attrs()

            # create the table
            table = self.handle.createTable(self.group, **options)

        else:
            table = self.table

        # validate the axes and set the kinds
        for a in self.axes:
            a.validate_and_set(table, append)

        # add the rows
        self.write_data()
        self.handle.flush()

    def write_data(self):
        """ fast writing of data: requires specific cython routines each axis shape """

        # create the masks & values
        masks  = []
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

        # the arguments
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

        # delete all rows (and return the nrows)
        if where is None or not len(where):
            nrows = self.nrows
            self.handle.removeNode(self.group, recursive=True)
            return nrows

        # infer the data kind
        if not self.infer_axes(): return None

        # create the selection
        table = self.table
        self.selection = Selection(self, where)
        values = self.selection.select_coords()

        # delete the rows in reverse order
        l  = Series(values).order()
        ln = len(l)

        if ln:

            # construct groups of consecutive rows
            diff   = l.diff()
            groups = list(diff[diff>1].index)

            # 1 group
            if not len(groups):
                groups = [0]

            # final element
            if groups[-1] != ln:
                groups.append(ln)

            # initial element
            if groups[0] != 0:
                groups.insert(0,0)

            # we must remove in reverse order!
            pg = groups.pop()
            for g in reversed(groups):
                rows = l.take(range(g,pg))
                table.removeRows(start = rows[rows.index[0]], stop = rows[rows.index[-1]]+1)
                pg = g

            self.handle.flush()

        # return the number of rows removed
        return ln

class AppendableFrameTable(AppendableTable):
    """ suppor the new appendable table formats """
    table_type = 'appendable_frame'
    ndim       = 2
    obj_type   = DataFrame

    @property
    def is_transposed(self):
        return self.index_axes[0].axis == 1

    def get_data_blocks(self, obj):
        """ these are written transposed """
        if self.is_transposed:
            obj = obj.T
        return obj._data.blocks

    def read(self, where=None):

        if not self.read_axes(where): return None

        index   = self.index_axes[0].values
        frames  = []
        for a in self.values_axes:
            columns = Index(a.values)

            if self.is_transposed:
                values   = a.cvalues
                index_   = columns
                columns_ = index
            else:
                values   = a.cvalues.T
                index_   = index
                columns_ = columns

            block   = make_block(values, columns_, columns_)
            mgr     = BlockManager([ block ], [ columns_, index_ ])
            frames.append(DataFrame(mgr))
        df = concat(frames, axis = 1, verify_integrity = True)

        # sort the indicies & reorder the columns
        for axis,labels in self.non_index_axes:
            df = df.reindex_axis(labels,axis=axis,copy=False)

        # apply the selection filters & axis orderings
        df = self.process_axes(df)

        return df

class AppendablePanelTable(AppendableTable):
    """ suppor the new appendable table formats """
    table_type = 'appendable_panel'
    ndim       = 3
    obj_type   = Panel

    def get_data_blocks(self, obj):
        """ these are written transposed """
        if self.is_transposed:
            obj = obj.transpose(*self.data_orientation)
        return obj._data.blocks

    @property
    def is_transposed(self):
        return self.data_orientation != tuple(range(self.ndim))

class AppendableNDimTable(AppendablePanelTable):
    """ suppor the new appendable table formats """
    table_type = 'appendable_ndim'
    ndim       = 4
    obj_type   = Panel4D

# table maps
_TABLE_MAP = {
    'appendable_frame' : AppendableFrameTable,
    'appendable_panel' : AppendablePanelTable,
    'appendable_ndim'  : AppendableNDimTable,
    'worm'             : WORMTable,
    'legacy_frame'     : LegacyFrameTable,
    'legacy_panel'     : LegacyPanelTable,
    'default'          : AppendablePanelTable,
}

def create_table(parent, group, typ = None, **kwargs):
    """ return a suitable Table class to operate """

    pt = getattr(group._v_attrs,'pandas_type',None)
    tt = getattr(group._v_attrs,'table_type',None) or typ

    # a new node
    if pt is None:

        return (_TABLE_MAP.get(typ) or _TABLE_MAP.get('default'))(parent, group, **kwargs)

    # existing node (legacy)
    if tt is None:

        # distiguish between a frame/table
        tt = 'legacy_panel'
        try:
            fields = group.table._v_attrs.fields
            if len(fields) == 1 and fields[0] == 'value':
                tt = 'legacy_frame'
        except:
            pass

    return _TABLE_MAP.get(tt)(parent, group, **kwargs)


def _itemsize_string_array(arr):
    """ return the maximum size of elements in a strnig array """
    return max([ str_len(arr[v].ravel()).max() for v in range(arr.shape[0]) ])

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
        queryables : a kinds map (dict of column name -> kind), or None i column is non-indexable

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
        Term('major_axis>20121114')
        Term('minor_axis', ['A','B'])

    """

    _ops     = ['<=','<','>=','>','!=','=']
    _search  = re.compile("^(?P<field>\w+)(?P<op>%s)(?P<value>.+)$" % '|'.join(_ops))

    def __init__(self, field, op = None, value = None, queryables = None):
        self.field = None
        self.op    = None
        self.value = None
        self.q     = queryables or dict()
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

        # we have valid conditions
        if self.op in ['>','>=','<','<=']:
            if hasattr(self.value,'__iter__') and len(self.value) > 1:
                raise Exception("an inequality condition cannot have multiple values [%s]" % str(self))

        if not hasattr(self.value,'__iter__'):
            self.value = [ self.value ]

        if len(self.q):
            self.eval()

    def __str__(self):
        return "field->%s,op->%s,value->%s" % (self.field,self.op,self.value)

    __repr__ = __str__

    @property
    def is_valid(self):
        """ return True if this is a valid field """
        return self.field in self.q

    @property
    def is_in_table(self):
        """ return True if this is a valid column name for generation (e.g. an actual column in the table) """
        return self.q.get(self.field) is not None

    @property
    def kind(self):
        """ the kind of my field """
        return self.q.get(self.field)

    def eval(self):
        """ set the numexpr expression for this term """

        if not self.is_valid:
            raise Exception("query term is not valid [%s]" % str(self))

        # convert values if we are in the table
        if self.is_in_table:
            values = [ self.convert_value(v) for v in self.value ]
        else:
            values = [ [v, v] for v in self.value ]

        # equality conditions
        if self.op in ['=','!=']:

            if self.is_in_table:

                # too many values to create the expression?
                if len(values) <= 61:
                    self.condition = "(%s)" % ' | '.join([ "(%s == %s)" % (self.field,v[0]) for v in values])

                # use a filter after reading
                else:
                    self.filter = (self.field,Index([ v[1] for v in values ]))

            else:

                self.filter = (self.field,Index([ v[1] for v in values ]))

        else:

            if self.is_in_table:

                self.condition = '(%s %s %s)' % (self.field, self.op, values[0][0])

            else:

                raise Exception("passing a filterable condition to a non-table indexer [%s]" % str(self))

    def convert_value(self, v):

        #### a little hacky here, need to really figure out what we should convert ####x
        if self.field == 'index' or self.field == 'major_axis':
            if self.kind == 'datetime64' :
                return [lib.Timestamp(v).value, None]
            elif isinstance(v, datetime) or hasattr(v,'timetuple') or self.kind == 'date':
                return [time.mktime(v.timetuple()), None]
            elif self.kind == 'integer':
                v = int(float(v))
                return [v, v]
            elif self.kind == 'float':
                v = float(v)
                return [v, v]
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
        self.condition  = None
        self.filter     = None
        self.terms      = self.generate(where)

        # create the numexpr & the filter
        if self.terms:
            conds = [ t.condition for t in self.terms if t.condition is not None ]
            if len(conds):
                self.condition = "(%s)" % ' & '.join(conds)
            self.filter = []
            for t in self.terms:
                if t.filter is not None:
                    self.filter.append(t.filter)

    def generate(self, where):
        """ where can be a : dict,list,tuple,string """
        if where is None: return None

        if not isinstance(where, (list,tuple)):
            where = [ where ]
        else:
            # do we have all list/tuple
            if not any([ isinstance(w, (list,tuple,Term)) for w in where ]):
                where = [ where ]

        queryables = self.table.queryables()
        return [ Term(c, queryables = queryables) for c in where ]

    def select(self):
        """
        generate the selection
        """
        if self.condition is not None:
            return self.table.table.readWhere(self.condition)
        else:
            return self.table.table.read()

    def select_coords(self):
        """
        generate the selection
        """
        return self.table.table.getWhereList(self.condition, sort = True)


def _get_index_factory(klass):
    if klass == DatetimeIndex:
        def f(values, freq=None, tz=None):
            return DatetimeIndex._simple_new(values, None, freq=freq,
                                             tz=tz)
        return f
    return klass
